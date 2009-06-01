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
     &     ,sname_strlen,units_strlen,lname_strlen
      IMPLICIT NONE
      SAVE
      INTEGER, PARAMETER :: KOIJ=21,KOIJL=26,KOL=6,KOLNST=8
!@var OIJ   lat-lon ocean diagnostics (on ocean grid)
!@var OIJL  3-dimensional ocean diagnostics
!@var OL    vertical ocean diagnostics
!@var OLNST strait diagnostics
!ny?  logical :: allocated_odiag_glob = .false.
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:)  :: OIJ_loc   !ny? ,OIJ
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: OIJL_loc !ny? ,OIJL
      REAL*8, DIMENSION(:,:,:), allocatable      :: OIJ
      REAL*8, DIMENSION(:,:,:,:), allocatable :: OIJL
      REAL*8, DIMENSION(LMO,KOL)   :: OL
      REAL*8, DIMENSION(LMO,NMST,KOLNST):: OLNST

!@var OIJL_out like OIJL_loc, but rescaled for postprocessing
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: OIJL_out

!@var IJ_xxx Names for OIJ diagnostics
      INTEGER IJ_HBL,IJ_BO,IJ_BOSOL,IJ_USTAR,IJ_SSH,IJ_PB,IJ_SF
!@var lname_oij Long names for OIJ diagnostics
      CHARACTER(len=lname_strlen), DIMENSION(KOIJ) :: LNAME_OIJ
!@var sname_oij Short names for OIJ diagnostics
      CHARACTER(len=sname_strlen), DIMENSION(KOIJ) :: SNAME_OIJ
!@var units_oij Units for OIJ diagnostics
      CHARACTER(len=units_strlen), DIMENSION(KOIJ) :: UNITS_OIJ
!@var ia_oij IDACC numbers for OIJ diagnostics
      INTEGER, DIMENSION(KOIJ) :: IA_OIJ
!@var scale_oij scales for OIJ diagnostics
      REAL*8, DIMENSION(KOIJ) :: SCALE_OIJ
!@var [ij]grid_oij Grid descriptor for OIJ diagnostics
      INTEGER, DIMENSION(KOIJ) :: IGRID_OIJ,JGRID_OIJ
#ifdef TRACERS_OceanBiology
!@var ij_pCO2 surface ocean partial CO2 pressure
       INTEGER :: IJ_dic,IJ_pCO2,IJ_nitr,IJ_diat,ij_herb
     .           ,ij_amm,ij_sil,ij_iron,ij_chlo,ij_cyan
     .           ,ij_cocc,ij_doc,IJ_alk
     .           ,ij_flux
#endif
!@var CDL_OIJ consolidated metadata for OIJ output fields in CDL notation
      character(len=100), dimension(koij*6) :: cdl_oij

!@var IJL_xxx Names for OIJL diagnostics
      INTEGER IJL_MO,IJL_G0M,IJL_S0M,IJL_GFLX,IJL_SFLX,IJL_MFU,IJL_MFV
     *     ,IJL_MFW,IJL_GGMFL,IJL_SGMFL,IJL_KVM,IJL_KVG,IJL_WGFL
     *     ,IJL_WSFL,IJL_PTM,IJL_PDM,IJL_MOU,IJL_MOV
!@var lname_oijl Long names for OIJL diagnostics
      CHARACTER(len=lname_strlen), DIMENSION(KOIJL) :: LNAME_OIJL
!@var sname_oijl Short names for OIJL diagnostics
      CHARACTER(len=sname_strlen), DIMENSION(KOIJL) :: SNAME_OIJL
!@var units_oijl Units for OIJL diagnostics
      CHARACTER(len=units_strlen), DIMENSION(KOIJL) :: UNITS_OIJL
!@var ia_oijl IDACC numbers for OIJL diagnostics
      INTEGER, DIMENSION(KOIJL) :: IA_OIJL
!@var denom_oijl denominators for OIJL diagnostics
      INTEGER, DIMENSION(KOIJL) :: DENOM_OIJL
!@var scale_oijl scales for OIJL diagnostics
      REAL*8, DIMENSION(KOIJL) :: SCALE_OIJL
!@var [ijl]grid_oijl Grid descriptors for OIJL diagnostics
      INTEGER, DIMENSION(KOIJL) :: IGRID_OIJL,JGRID_OIJL,LGRID_OIJL
!@var CDL_OIJL consolidated metadata for OIJL output fields in CDL notation
      character(len=100), dimension(koijl*6) :: cdl_oijl

!@var LN_xxx Names for OLNST diagnostics
      INTEGER LN_KVM,LN_KVG,LN_WGFL,LN_WSFL,LN_MFLX,LN_GFLX,LN_SFLX
     *     ,LN_ICFL
!@var lname_olnst Long names for OLNST diagnostics
      CHARACTER(len=lname_strlen), DIMENSION(KOLNST) :: LNAME_OLNST
!@var sname_olnst Short names for OLNST diagnostics
      CHARACTER(len=sname_strlen), DIMENSION(KOLNST) :: SNAME_OLNST
!@var units_olnst Units for OLNST diagnostics
      CHARACTER(len=units_strlen), DIMENSION(KOLNST) :: UNITS_OLNST
!@var ia_olnst IDACC numbers for OLNST diagnostics
      INTEGER, DIMENSION(KOLNST) :: IA_OLNST
!@var scale_olnst scales for OLNST diagnostics
      REAL*8, DIMENSION(KOLNST) :: SCALE_OLNST
!@var lgrid_olnst Grid descriptors for OLNST diagnostics
      INTEGER, DIMENSION(KOLNST) :: LGRID_OLNST
!@var CDL_OLNST consolidated metadata for OLNST output fields in CDL notation
      character(len=100), dimension(kolnst*6) :: cdl_olnst

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
!@var NBAS number of ocean basins 
!@var nqty number of output qtys zonally averaged over basins
      INTEGER, PARAMETER :: NBAS=4,nqty=3,NCIRC=3
!@var BASIN names of ocean basins for diag output
      CHARACTER*16, DIMENSION(NBAS) :: BASIN=
     *     (/"Atlantic","Pacific ","Indian  ","Global  "/)
      character(len=4), dimension(nqty), parameter :: qtyname=(/
     &     'Mass','Heat','Salt' /)
      character(len=9), dimension(nqty), parameter :: qtyflxunit=(/
     &     '10^9 kg/s','10^15 W  ','10^6 kg/s' /)
      character(len=3), dimension(ncirc), parameter :: circstr=(/
     &     'moc','gmf','gyr' /)
      character(len=7), dimension(ncirc), parameter :: circname=(/
     &     'Overtrn', 'GM flx ', 'Hor gyr'/)

!@var KOJL,OJL (number of qtys having) zonal sums/means over basins
      INTEGER, PARAMETER :: KOJL=4*NBAS
      REAL*8, DIMENSION(JM,LMO,NBAS,KOJL/NBAS) :: OJL
      REAL*8, DIMENSION(JM,LMO,KOJL) :: OJL_out

!@var JL_xxx indices for qtys in OJL
      INTEGER :: JL_M,JL_PT,JL_S,JL_SF

!@var lname_ojl Long names for OJL diagnostics
      CHARACTER(len=lname_strlen), DIMENSION(KOJL) :: LNAME_OJL
!@var sname_ojl Short names for OJL diagnostics
      CHARACTER(len=sname_strlen), DIMENSION(KOJL) :: SNAME_OJL
!@var units_ojl Units for OJL diagnostics
      CHARACTER(len=units_strlen), DIMENSION(KOJL) :: UNITS_OJL
!@var ia_ojl IDACC numbers for OJL diagnostics
      INTEGER, DIMENSION(KOJL) :: IA_OJL
!@var denom_ojl denominators for OJL diagnostics
      INTEGER, DIMENSION(KOJL) :: DENOM_OJL
!@var scale_ojl scales for OJL diagnostics
      REAL*8, DIMENSION(KOJL) :: SCALE_OJL
!@var [jl]grid_ojl Grid descriptors for OJL diagnostics
      INTEGER, DIMENSION(KOJL) :: JGRID_OJL,LGRID_OJL
!@var CDL_OJL consolidated metadata for OJL output fields in CDL notation
      character(len=100), dimension(KOJL*6) :: cdl_ojl

!@var NSEC number of lat/lon sections for diags
      INTEGER, PARAMETER :: NSEC=3
!@var SEC_LAT, SEC_LON lat/lon for sectional tracer profiles
      REAL*8, PARAMETER :: SEC_LAT(NSEC) = (/-64.,0.,48./),
     *     SEC_LON(NSEC) = (/-165.0,-30.,65./)
C**** OTJ is integrated flux array OTJ(LATITUDE,BASIN,KQ)
C****   KQ 1   Mass (kg)
C****      2   Heat (J)
C****      3   Salt (kg)
C****  OTJCOMP OTJ(LATITUDE,BASIN,COMP,KQ) (KQ = 2,3)
C**** COMP 1   advected by overturning
C****      2   flux from GM
C****      3   advected by horizontal gyres (residual)
C****
      REAL*8 OTJ(0:JM,4,3),OTJCOMP(0:JM,4,3,3)
      integer, parameter :: kotj=4*4*3

!@var OTJ_out reshaped combination of OTJ and OTJCOMP
      REAL*8 OTJ_out(JM,kotj)

!@var ia_otj IDACC numbers for OTJ diagnostics
      integer, dimension(kotj) :: ia_otj
!@var scale_otj scales for OTJ diagnostics
      real*8, dimension(kotj) :: scale_otj
!@var sname_otj short names for OTJ diagnostics
      character(len=sname_strlen), dimension(kotj) :: sname_otj
!@var lname_otj Long names for OTJ diagnostics
      character(len=lname_strlen), dimension(kotj) :: lname_otj
!@var units_otj units for OTJ diagnostics
      character(len=units_strlen), dimension(kotj) :: units_otj
!@var CDL_OTJ consolidated metadata for OTJ output fields in CDL notation
      character(len=100), dimension(kotj*6) :: cdl_otj

!@var SFM meridional overturning stream function for each basin
      REAL*8, DIMENSION(JM,0:LMO,4) :: SFM!,SFS SFS is for salt

C****

#ifdef TRACERS_OCEAN
!@var KTOIJL number of 3-dimensional ocean tracer diagnostics
      INTEGER, PARAMETER :: KTOIJL=10
!@var TOIJL  3-dimensional ocean tracer diagnostics
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:,:) :: TOIJL_loc !ny? ,TOIJL
      REAL*8, DIMENSION(:,:,:,:,:), ALLOCATABLE    :: TOIJL
!@var toijl_xxx indices for TOIJL diags
      INTEGER, PARAMETER :: toijl_conc=1,toijl_tflx=2,toijl_gmfl=6
     *     ,toijl_wtfl=10
!@var TLNST strait diagnostics
      REAL*8, DIMENSION(LMO,NMST,KOLNST,NTM):: TLNST
#endif

#ifdef NEW_IO
c declarations that facilitate summation of acc-files when using
c the new i/o system
      target :: oij_loc,oijl_loc,ol,olnst,oijl_out
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
      REAL*4, DIMENSION(:,:,:), ALLOCATABLE  :: OIJ4
      REAL*4, DIMENSION(:,:,:,:), ALLOCATABLE :: OIJL4
      REAL*4, DIMENSION(LMO,KOL)   :: OL4
      REAL*4, DIMENSION(LMO,NMST,KOLNST):: OLNST4
#ifdef TRACERS_OCEAN
!@var TOIJL4 work array for read/wrote operations
      REAL*4, DIMENSION(:,:,:,:,:), ALLOCATABLE :: TOIJL4
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
            allocate(OIJ4(IM,JM,KOIJ),OIJL4(IM,JM,LMO,KOIJL))
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
            deallocate(OIJ4,OIJL4)
            IF (HEADER(1:LHEAD).NE.MODULE_HEADER(1:LHEAD)) THEN
              PRINT*,"Discrepancy in module version ",HEADER
     *           ,MODULE_HEADER
              GO TO 10
            END IF
#ifdef TRACERS_OCEAN
            allocate(TOIJL4(IM,JM,LMO,KTOIJL,NTM))
            READ (kunit,err=10) TR_HEADER,TOIJL4,TLNST4,it
C**** accumulate diagnostics
            TOIJL=TOIJL+TOIJL4
            TLNST=TLNST+TLNST4
            deallocate(TOIJL4)
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
      use odiag, only : ol,olnst,oij=>oij_ioptr,oijl=>oijl_ioptr
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

      subroutine def_meta_ocdiag(fid)
!@sum  def_meta_ocdiag defines metadata in ocean acc files
!@auth M. Kelley
!@ver  beta
      use odiag
      use pario, only : defvar,write_attr
      USE OCEANR_DIM, only : grid=>ogrid
      use diag_com, only : zoc,zoc1
      use ocean, only : olat_dg,olon_dg
      use straits, only : nmst,lmst,name_st
      implicit none
      integer :: fid         !@var fid file id

      call defvar(grid,fid,name_st,'strait_name(strait_strlen,nmst)')
      call defvar(grid,fid,lmst,'lmst(nmst)')

      call defvar(grid,fid,olat_dg(:,1),'lato(jmo)')
      call defvar(grid,fid,olat_dg(:,2),'lato2(jmo)')
      call defvar(grid,fid,olon_dg(:,1),'lono(imo)')
      call defvar(grid,fid,olon_dg(:,2),'lono2(imo)')
      call defvar(grid,fid,zoc(1:lmo),'zoc(lmo)')
      call defvar(grid,fid,zoc1(1:lmo),'zoce(lmo)')

      call write_attr(grid,fid,'oij','reduction','sum')
      call write_attr(grid,fid,'oij','split_dim',3)
      call defvar(grid,fid,ia_oij,'ia_oij(koij)')
      call defvar(grid,fid,scale_oij,'scale_oij(koij)')
      call defvar(grid,fid,sname_oij,'sname_oij(sname_strlen,koij)')
      call defvar(grid,fid,cdl_oij,'cdl_oij(cdl_strlen,kcdl_oij)')

      call write_attr(grid,fid,'oijl','reduction','sum')
      call write_attr(grid,fid,'oijl','split_dim',4)
      call defvar(grid,fid,ia_oijl,'ia_oijl(koijl)')
      call defvar(grid,fid,denom_oijl,'denom_oijl(koijl)')
      call defvar(grid,fid,scale_oijl,'scale_oijl(koijl)')
      call defvar(grid,fid,sname_oijl,'sname_oijl(sname_strlen,koijl)')
      call defvar(grid,fid,cdl_oijl,'cdl_oijl(cdl_strlen,kcdl_oijl)')

      call write_attr(grid,fid,'ol','reduction','sum')
      call write_attr(grid,fid,'ol','split_dim',2)

      call write_attr(grid,fid,'olnst','reduction','sum')
      call write_attr(grid,fid,'olnst','split_dim',3)
      call defvar(grid,fid,ia_olnst,'ia_olnst(kolnst)')
      call defvar(grid,fid,scale_olnst,'scale_olnst(kolnst)')
      call defvar(grid,fid,sname_olnst,
     &     'sname_olnst(sname_strlen,kolnst)')
      call defvar(grid,fid,cdl_olnst,'cdl_olnst(cdl_strlen,kcdl_olnst)')

      call defvar(grid,fid,ojl_out,'ojl(jmo,lmo,kojl)',
     &     r4_on_disk=.true.)
      call write_attr(grid,fid,'ojl','reduction','sum')
      call write_attr(grid,fid,'ojl','split_dim',3)
      call defvar(grid,fid,ia_ojl,'ia_ojl(kojl)')
      call defvar(grid,fid,denom_ojl,'denom_ojl(kojl)')
      call defvar(grid,fid,scale_ojl,'scale_ojl(kojl)')
      call defvar(grid,fid,sname_ojl,'sname_ojl(sname_strlen,kojl)')
      call defvar(grid,fid,cdl_ojl,'cdl_ojl(cdl_strlen,kcdl_ojl)')

      call defvar(grid,fid,otj_out,'otj(jmo,kotj)',
     &     r4_on_disk=.true.)
      call write_attr(grid,fid,'otj','reduction','sum')
      call write_attr(grid,fid,'otj','split_dim',2)
      call defvar(grid,fid,ia_otj,'ia_otj(kotj)')
      call defvar(grid,fid,scale_otj,'scale_otj(kotj)')
      call defvar(grid,fid,sname_otj,'sname_otj(sname_strlen,kotj)')
      call defvar(grid,fid,cdl_otj,'cdl_otj(cdl_strlen,kcdl_otj)')

      return
      end subroutine def_meta_ocdiag

      subroutine write_meta_ocdiag(fid)
!@sum  write_meta_ocdiag write ocean accumulation metadata to file
!@auth M. Kelley
      use odiag
      use pario, only : write_dist_data,write_data
      USE OCEANR_DIM, only : grid=>ogrid
      use diag_com, only : zoc,zoc1
      use ocean, only : olat_dg,olon_dg
      use straits, only : nmst,lmst,name_st
      implicit none
      integer :: fid         !@var fid file id

      call write_data(grid,fid,'strait_name',name_st)
      call write_data(grid,fid,'lmst',lmst)

      call write_data(grid,fid,'lato',olat_dg(:,1))
      call write_data(grid,fid,'lato2',olat_dg(:,2))
      call write_data(grid,fid,'lono',olon_dg(:,1))
      call write_data(grid,fid,'lono2',olon_dg(:,2))
      call write_data(grid,fid,'zoc',zoc(1:lmo))
      call write_data(grid,fid,'zoce',zoc1(2:lmo+1))

      call write_data(grid,fid,'ia_oij',ia_oij)
      call write_data(grid,fid,'scale_oij',scale_oij)
      call write_data(grid,fid,'sname_oij',sname_oij)
      call write_data(grid,fid,'cdl_oij',cdl_oij)

      call write_data(grid,fid,'ia_oijl',ia_oijl)
      call write_data(grid,fid,'denom_oijl',denom_oijl)
      call write_data(grid,fid,'scale_oijl',scale_oijl)
      call write_data(grid,fid,'sname_oijl',sname_oijl)
      call write_data(grid,fid,'cdl_oijl',cdl_oijl)

      call write_data(grid,fid,'ojl',ojl_out)
      call write_data(grid,fid,'ia_ojl',ia_ojl)
      call write_data(grid,fid,'denom_ojl',denom_ojl)
      call write_data(grid,fid,'scale_ojl',scale_ojl)
      call write_data(grid,fid,'sname_ojl',sname_ojl)
      call write_data(grid,fid,'cdl_ojl',cdl_ojl)

      call write_data(grid,fid,'otj',otj_out)
      call write_data(grid,fid,'ia_otj',ia_otj)
      call write_data(grid,fid,'scale_otj',scale_otj)
      call write_data(grid,fid,'sname_otj',sname_otj)
      call write_data(grid,fid,'cdl_otj',cdl_otj)

      call write_data(grid,fid,'ia_olnst',ia_olnst)
      call write_data(grid,fid,'scale_olnst',scale_olnst)
      call write_data(grid,fid,'sname_olnst',sname_olnst)
      call write_data(grid,fid,'cdl_olnst',cdl_olnst)

      return
      end subroutine write_meta_ocdiag

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

      subroutine set_ioptrs_ocnacc_extended
c point i/o pointers for diagnostic accumlations to the
c instances of the arrays containing derived quantities
      use odiag
      implicit none
      oijl_ioptr   => oijl_out
      return
      end subroutine set_ioptrs_ocnacc_extended

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
#ifndef CUBE_GRID
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
#endif
C****
      RETURN
      END SUBROUTINE DIAGCO


      SUBROUTINE conserv_ODIAG (M,CONSFN,ICON)
!@sum  conserv_ODIAG generic routine keeps track of conserved properties
!@+    uses OJ_BUDG mapping from ocean sub-domain to budget grid
!@auth Gary Russell/Gavin Schmidt/Denis Gueyffier
!@ver  1.0
      USE OCEAN, only : oJ_BUDG, oWTBUDG, oJ_0B, oJ_1B,imaxj
      USE DIAG_COM, only : consrv=>consrv_loc,nofm,jm_budg
      USE DOMAIN_DECOMP_1D, only : GET
      USE OCEANR_DIM, only : oGRID
      IMPLICIT NONE
!@var M index denoting from where routine is called
      INTEGER, INTENT(IN) :: M
!@var ICON index for the quantity concerned
      INTEGER, INTENT(IN) :: ICON
!@var CONSFN external routine that calculates total conserved quantity
      EXTERNAL CONSFN
!@var TOTAL amount of conserved quantity at this time
      REAL*8, DIMENSION(oGRID%I_STRT_HALO:oGRID%I_STOP_HALO,
     &                  oGRID%J_STRT_HALO:oGRID%J_STOP_HALO) :: TOTAL
      REAL*8, DIMENSION(JM_BUDG) :: TOTALJ
      INTEGER :: I,J,NM,NI
      INTEGER :: J_0,J_1

      CALL GET(ogrid, J_STRT=J_0, J_STOP=J_1)

C**** NOFM contains the indexes of the CONSRV array where each
C**** change is to be stored for each quantity. If NOFM(M,ICON)=0,
C**** no calculation is done.
C**** NOFM(1,ICON) is the index for the instantaneous value.
      IF (NOFM(M,ICON).gt.0) THEN
C**** Calculate current value TOTAL
        CALL CONSFN(TOTAL)
        NM=NOFM(M,ICON)
        NI=NOFM(1,ICON)
C**** Calculate zonal sums
        TOTALJ(oJ_0B:oJ_1B)=0.
        DO J=J_0,J_1
          DO I=1,IMAXJ(J) 
            TOTALJ(oJ_BUDG(I,J)) = TOTALJ(oJ_BUDG(I,J)) + TOTAL(I,J)
     &           *oWTBUDG(I,J)
          END DO
        END DO
C**** Accumulate difference from last time in CONSRV(NM)
        IF (M.GT.1) THEN
          DO J=oJ_0B,oJ_1B
            CONSRV(J,NM)=CONSRV(J,NM)+(TOTALJ(J)-CONSRV(J,NI))
          END DO
        END IF
C**** Save current value in CONSRV(NI)
        DO J=oJ_0B,oJ_1B
          CONSRV(J,NI)=TOTALJ(J)
        END DO
      END IF
      RETURN
C****
      END SUBROUTINE conserv_ODIAG

 
      SUBROUTINE init_ODIAG
!@sum  init_ODIAG initialises ocean diagnostics
!@auth Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only : rhows
      USE CONSTANT, only : bygrav
      USE MODEL_COM, only : dtsrc
      USE OCEAN, only : ze,dts,ndyno
      USE DIAG_COM, only : ia_src,conpt0,zoc,zoc1
      USE ODIAG
      use straits, only : nmst,name_st
      IMPLICIT NONE
      LOGICAL :: QCON(NPTS), T = .TRUE. , F = .FALSE.
      CHARACTER CONPT(NPTS)*10
      INTEGER k,kb,kq,kc,kk,n
      character(len=10) :: xstr,ystr,zstr

C**** Set names for OLNST diagnostics
c      LN_KVM=1  ; LN_KVG=2  ; LN_WGFL=3 ; LN_WSFL=4
c      LN_MFLX=5 ; LN_GFLX=6 ; LN_SFLX=7 ; LN_ICFL=8

C**** set properties for OLNST diagnostics
      do k=1,kolnst
        sname_olnst(k) = 'unused'
        lname_olnst(k) = 'no output'
        units_olnst(k) = 'no output'
        ia_olnst(k) = ia_src
        scale_olnst(k) = 1.
        lgrid_olnst(k) = 1
      enddo
      k=0
c
      k=k+1
      LN_KVM = k
      sname_olnst(k) = 'kvm'
      lname_olnst(k) = 'Vertical Diffusion Coefficient'
      units_olnst(k) = 'm^2/s'
      scale_olnst(k) = .5
      lgrid_olnst(k) = 2
c
      k=k+1
      LN_KVG = k
c
      k=k+1
      LN_WGFL = k
c
      k=k+1
      LN_WSFL = k
c
      k=k+1
      LN_MFLX = k
      sname_olnst(k) = 'mflx'
      lname_olnst(k) = 'Strait Transport of Mass'
      units_olnst(k) = '10^6 kg/s'
      scale_olnst(k) = 1.D-6 / DTS
c
      k=k+1
      LN_GFLX = k
      sname_olnst(k) = 'hflx'
      lname_olnst(k) = 'Strait Trans of Potential Enthalpy'
      units_olnst(k) = '10^11 W'
      scale_olnst(k) = .5D-11 / DTS
c
      k=k+1
      LN_SFLX = k
      sname_olnst(k) = 'sflx'
      lname_olnst(k) = 'Strait Transport of Salt'
      units_olnst(k) = '10^5 kg/s'
      scale_olnst(k) = .5D-5 / DTS
c
      k=k+1
      LN_ICFL = k

C**** Set names for OL diagnostics
      L_RHO=1   ; L_TEMP=2  ; L_SALT=3

C**** set properties for OIJL diagnostics
      do k=1,koijl
        sname_oijl(k) = 'unused'
        ia_oijl(k) = ia_src
        denom_oijl(k) = 0
        scale_oijl(k) = 1.
        lname_oijl(k) = 'no output'
        units_oijl(k) = 'no output'
        igrid_oijl(k) = 1
        jgrid_oijl(k) = 1
        lgrid_oijl(k) = 1
      enddo
      k=0
      
c
      k=k+1
      IJL_MO = k
      sname_oijl(k) = 'mo' ! gridbox mass; serves as a denominator
c
      k=k+1
      IJL_MOU = k
      sname_oijl(k) = 'mou' ! denominator for east-west velocity
      igrid_oijl(k) = 2
c
      k=k+1
      IJL_MOV = k
      sname_oijl(k) = 'mov' ! denominator for north-south velocity
      jgrid_oijl(k) = 2
c
      k=k+1
      IJL_G0M = k
      denom_oijl(k) = IJL_MO
      sname_oijl(k) = 'heat'
      units_oijl(k) = '10^6 J/m^2'
      lname_oijl(k) = 'OCEAN HEAT CONTENT'
      scale_oijl(k) = 1d-6
c
      k=k+1
      IJL_S0M = k
      denom_oijl(k) = IJL_MO
      sname_oijl(k) = 'salt'
      units_oijl(k) = 'psu'
      lname_oijl(k) = 'OCEAN SALINITY'
      scale_oijl(k) = 1d3
c
      k=k+1
      IJL_MFU = k
      denom_oijl(k) = IJL_MOU
      sname_oijl(k) = 'u'
      units_oijl(k) = 'cm/s'
      lname_oijl(k) = 'EAST-WEST VELOCITY'
      scale_oijl(k) = 1d2*(2./ndyno)
      igrid_oijl(k) = 2
c
      k=k+1
      IJL_MFV = k
      denom_oijl(k) = IJL_MOV
      sname_oijl(k) = 'v'
      units_oijl(k) = 'cm/s'
      lname_oijl(k) = 'NORTH-SOUTH VELOCITY'
      scale_oijl(k) = 1d2*(2./ndyno)
      jgrid_oijl(k) = 2
c
      k=k+1
      IJL_MFW = k
      sname_oijl(k) = 'w'
      units_oijl(k) = 'cm/s'
      lname_oijl(k) = 'DOWNWARD VERTICAL VELOCITY'
      scale_oijl(k) = (1d2/RHOWS)*(2./ndyno)
      lgrid_oijl(k) = 2
c
      k=k+1
      IJL_GFLX = k
      sname_oijl(k) = 'gflx_x'
      units_oijl(k) = '10^15 W'
      lname_oijl(k) = "EAST-WEST HEAT FLUX"
      scale_oijl(k) = 1d-15/dts
      igrid_oijl(k) = 2
c
      k=k+1
c      IJL_GFLX_ns = k
      sname_oijl(k) = 'gflx_y'
      units_oijl(k) = '10^15 W'
      lname_oijl(k) = 'NORTH-SOUTH HEAT FLUX'
      scale_oijl(k) = 1d-15/dts
      jgrid_oijl(k) = 2
c
      k=k+1
c      IJL_GFLX_vert = k
c
      k=k+1
      IJL_SFLX = k
      units_oijl(k) = '10^6 kg/s'
      lname_oijl(k) = 'EAST-WEST SALT FLUX'
      scale_oijl(k) = 1d-6/dts
      igrid_oijl(k) = 2
c
      k=k+1
c      IJL_SFLX_ns = k
      units_oijl(k) = '10^6 kg/s'
      lname_oijl(k) = 'NORTH-SOUTH SALT FLUX'
      scale_oijl(k) = 1d-6/dts
      jgrid_oijl(k) = 2
c
      k=k+1
c      IJL_SFLX_vert = k
c
      k=k+1
      IJL_KVM = k
      sname_oijl(k) = 'kvm'
      units_oijl(k) = 'cm^2/s'
      lname_oijl(k) = 'VERT. MOM. DIFF.'
      scale_oijl(k) = 1d4*.00097d0**2/dts
      lgrid_oijl(k) = 2
c
      k=k+1
      IJL_KVG = k
      sname_oijl(k) = 'kvg'
      units_oijl(k) = 'cm^2/s'
      lname_oijl(k) = 'Vertical heat diffusivity'
      lgrid_oijl(k) = 2
c
      k=k+1
      IJL_WGFL = k
      sname_oijl(k) = 'wgfl'
      units_oijl(k) = 'W/m^2'
      lname_oijl(k) = 'VERT. HEAT DIFF.'
      scale_oijl(k) = 1./dts
      lgrid_oijl(k) = 2
c
      k=k+1
      IJL_WSFL = k
      sname_oijl(k) = 'wsfl'
      units_oijl(k) = '10^-6 kg/m^2'
      lname_oijl(k) = 'VERT. SALT DIFF.'
      scale_oijl(k) = 1d6/dts
      lgrid_oijl(k) = 2
c
      k=k+1
      IJL_GGMFL = k
      sname_oijl(k) = 'ggmflx_x'
      units_oijl(k) = '10^9 W'
      lname_oijl(k) = 'GM/EDDY E-W HEAT FLUX'
      scale_oijl(k) = 1d-9/dts
      igrid_oijl(k) = 2
c
      k=k+1
c      IJL_GGMFL_ns = k
      sname_oijl(k) = 'ggmflx_y'
      units_oijl(k) = '10^9 W'
      lname_oijl(k) = 'GM/EDDY N-S HEAT FLUX'
      scale_oijl(k) = 1d-9/dts
      jgrid_oijl(k) = 2
c
      k=k+1
c      IJL_GGMFL_vert = k
      sname_oijl(k) = 'ggmflx_z'
      units_oijl(k) = 'W/m^2'
      lname_oijl(k) = 'GM/EDDY VERT. HEAT FLUX'
      scale_oijl(k) = 1./dts
      lgrid_oijl(k) = 2
c
      k=k+1
      IJL_SGMFL = k
      sname_oijl(k) = 'sgmflx_x'
      units_oijl(k) = 'kg/s'
      lname_oijl(k) = 'GM/EDDY E-W SALT FLUX'
      scale_oijl(k) = 1./dts
      igrid_oijl(k) = 2
c
      k=k+1
c      IJL_SGMFL_ns = k
      sname_oijl(k) = 'sgmflx_y'
      units_oijl(k) = 'kg/s'
      lname_oijl(k) = 'GM/EDDY N-S SALT FLUX'
      scale_oijl(k) = 1./dts
      jgrid_oijl(k) = 2
c
      k=k+1
c      IJL_SGMFL_vert = k
      sname_oijl(k) = 'sgmflx_z'
      units_oijl(k) = 'kg/s'
      lname_oijl(k) = 'GM/EDDY VERT. SALT FLUX'
      scale_oijl(k) = 1./dts
      lgrid_oijl(k) = 2
c
      k=k+1
      IJL_PTM = k
      denom_oijl(k) = IJL_MO
      sname_oijl(k) = 'pot_temp'
      units_oijl(k) = 'C'
      lname_oijl(k) = 'OCEAN POTENTIAL TEMPERATURE'
c
      k=k+1
      IJL_PDM = k
      denom_oijl(k) = IJL_MO
      sname_oijl(k) = 'pot_dens'
      units_oijl(k) = 'KG/M^3'
      lname_oijl(k) = 'OCEAN POTENTIAL DENSITY (SIGMA_0)'
c

C**** set properties for OIJ diagnostics
      do k=1,koij
        sname_oij(k) = 'unused'
        lname_oij(k) = 'no output'
        units_oij(k) = 'no output'
        igrid_oij(k) = 1
        jgrid_oij(k) = 1
      enddo
      k=0
c
      k=k+1
      IJ_HBL=k
      lname_oij(k)="Ocean Boundary layer depth (KPP)"
      sname_oij(k)="oij_hbl"
      units_oij(k)="m"
      ia_oij(k)=ia_src
      scale_oij(k)=0.25

      k=k+1
      IJ_BO=k
      lname_oij(k)="Surface buoyancy forcing (KPP)"
      sname_oij(k)="oij_bo"
      units_oij(k)="10^-7 m^2/s^3"
      ia_oij(k)=ia_src
      scale_oij(k)=0.25*1d7

      k=k+1
      IJ_BOSOL=k
      lname_oij(k)="Surface solar buoyancy flux"
      sname_oij(k)="oij_bosol"
      units_oij(k)="10^-7 m^2/s^3"
      ia_oij(k)=ia_src
      scale_oij(k)=0.25*1d7

      k=k+1
      IJ_USTAR=k
      lname_oij(k)="Surface friction speed"
      sname_oij(k)="oij_ustar"
      units_oij(k)="m/s"
      ia_oij(k)=ia_src
      scale_oij(k)=0.25

#ifdef TRACERS_OceanBiology
      k=k+1
      IJ_nitr=k
      lname_oij(k)="Surface ocean Nitrates"
      sname_oij(k)="oij_nitr"
      units_oij(k)="uM"
      ia_oij(k)=ia_src
      scale_oij(k)=1

      k=k+1
      IJ_amm=k
      lname_oij(k)="Surface ocean Ammonium"
      sname_oij(k)="oij_amm"
      units_oij(k)="uM"
      ia_oij(k)=ia_src
      scale_oij(k)=1

      k=k+1
      IJ_sil=k
      lname_oij(k)="Surface ocean Silicate"
      sname_oij(k)="oij_sil"
      units_oij(k)="uM"
      ia_oij(k)=ia_src
      scale_oij(k)=1

      k=k+1
      IJ_iron=k
      lname_oij(k)="Surface ocean Iron"
      sname_oij(k)="oij_iron"
      units_oij(k)="nM"
      ia_oij(k)=ia_src
      scale_oij(k)=1

      k=k+1
      IJ_diat=k
      lname_oij(k)="Surface ocean Diatoms"
      sname_oij(k)="oij_diat"
      units_oij(k)="mg/m3"
      ia_oij(k)=ia_src
      scale_oij(k)=1

      k=k+1
      IJ_chlo=k
      lname_oij(k)="Surface ocean Chlorophytes"
      sname_oij(k)="oij_chlo"
      units_oij(k)="mg/m3"
      ia_oij(k)=ia_src
      scale_oij(k)=1

      k=k+1
      IJ_cyan=k
      lname_oij(k)="Surface ocean Cyanobacteria"
      sname_oij(k)="oij_cyan"
      units_oij(k)="mg/m3"
      ia_oij(k)=ia_src
      scale_oij(k)=1

      k=k+1
      IJ_cocc=k
      lname_oij(k)="Surface ocean Coccolithophores"
      sname_oij(k)="oij_cocc"
      units_oij(k)="mg/m3"
      ia_oij(k)=ia_src
      scale_oij(k)=1

      k=k+1
      IJ_herb=k
      lname_oij(k)="Surface ocean Herbivores"
      sname_oij(k)="oij_herb"
      units_oij(k)="mg/m3"
      ia_oij(k)=ia_src
      scale_oij(k)=1

      k=k+1
      IJ_doc=k
      lname_oij(k)="Surface ocean DOC"
      sname_oij(k)="oij_doc"
      units_oij(k)="uM"
      ia_oij(k)=ia_src
      scale_oij(k)=1

      k=k+1
      IJ_dic=k
      lname_oij(k)="Surface ocean DIC"
      sname_oij(k)="oij_dic"
      units_oij(k)="uM"
      ia_oij(k)=ia_src
      scale_oij(k)=1

      k=k+1
      IJ_pCO2=k
      lname_oij(k)="Surface ocean partial CO2 pressure"
      sname_oij(k)="oij_pCO2"
      units_oij(k)="uatm"
      ia_oij(k)=ia_src
      scale_oij(k)=1

      k=k+1
      IJ_alk=k
      lname_oij(k)="Surface ocean alkalinity"
      sname_oij(k)="oij_alk"
      units_oij(k)="umol/kg"
      ia_oij(k)=ia_src
      scale_oij(k)=1
    
      k=k+1
      IJ_flux=k
      lname_oij(k)="Air_Sea Flux CO2 (watson)"
      sname_oij(k)="oij_flux"
      units_oij(k)="mole/m2/s"
      ia_oij(k)=ia_src
      scale_oij(k)=1

#endif

      k=k+1
      IJ_SSH=k
      lname_oij(k)="Ocean surface height"
      sname_oij(k)="oij_ssh"
      units_oij(k)="m"
      ia_oij(k)=ia_src
      scale_oij(k)=bygrav

      k=k+1
      IJ_PB=k
      lname_oij(k)="Ocean bottom pressure anomaly"
      sname_oij(k)="oij_pb"
      units_oij(k)="Pa"
      ia_oij(k)=ia_src
      scale_oij(k)=1.

      k=k+1
      IJ_SF=k
      lname_oij(k)='HORIZONTAL MASS TRANSPORT STREAMFUNCTION'
      sname_oij(k)='osfij'
      units_oij(k)='Sv'
      ia_oij(k)=ia_src
      scale_oij(k)=1.
      igrid_oij(k) = 2
      jgrid_oij(k) = 2

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
      CALL SET_CON(QCON,CONPT,"OCN MASS","(10**2 KG/M^2)  ",
     *     "(10**-8 KG/SM^2)",1d-2,1d8,icon_OMS)
C**** Oceanic angular momentum
      QCON=(/ F, F, F, T, F, F, F, T, T, T, T/)
      CALL SET_CON(QCON,CONPT,"OCN AM  ","(10**12 JS/M^2) ",
     *     "(10**2 J/M^2)   ",1d-12,1d-2,icon_OAM)
C**** Oceanic kinetic energy
      QCON=(/ F, F, F, T, F, F, F, T, T, T, T/)
      CALL SET_CON(QCON,CONPT,"OCEAN KE","(J/M^2)         ",
     *     "(10**-6 W/M^2)  ",1d0,1d6,icon_OKE)
C**** Oceanic potential enthalpy (heat)
      QCON=(/ F, F, F, T, F, F, F, T, T, T, T/)
      CALL SET_CON(QCON,CONPT,"OCN HEAT","(10**6 J/M^2)   ",
     *     "(10**-2 W/M^2)  ",1d-6,1d2,icon_OCE)
C**** Oceanic salt mass
      QCON=(/ F, F, F, T, F, F, F, T, T, T, T/)
      CALL SET_CON(QCON,CONPT,"OCN SALT","(10 KG/M^2)     ",
     *     "(10**-9 KG/SM^2)",1d-1,1d9,icon_OSL)
C**** Initialise ocean basins
      CALL OBASIN

C**** Define ocean depths for diagnostic output
      ZOC1(1:LMO+1) = ZE(0:LMO)
      ZOC(1:LMO) = 0.5*(ZE(1:LMO)+ZE(0:LMO-1))

C**** Metadata for northward transports
      scale_otj(:) = 1.
      ia_otj(:) = ia_src
      kk = 0
      do kq=1,3
      do kb=1,4
        kk = kk + 1
        sname_otj(kk)='nt_'//trim(qtyname(kq))//'_'//basin(kb)(1:3)
        lname_otj(kk)='North. Trans. of '//trim(qtyname(kq))//
     &       ' in the '//trim(basin(kb))//' basin'
        units_otj(kk)=qtyflxunit(kq)
        do kc=1,3
          kk = kk + 1
          sname_otj(kk)='nt_'//trim(qtyname(kq))//'_'//basin(kb)(1:3)//
     &       '_'//trim(circstr(kc))
          lname_otj(kk)='North. Trans. of '//trim(qtyname(kq))//
     &       ' in the '//trim(basin(kb))//' basin by '//
     &       trim(circname(kc))
          units_otj(kk)=qtyflxunit(kq)
        enddo
      enddo
      enddo

C**** set properties for OJL diagnostics
      do k=1,kojl
        sname_ojl(k) = 'unused'
        ia_ojl(k) = ia_src
        denom_ojl(k) = 0
        scale_ojl(k) = 1.
        lname_ojl(k) = 'no output'
        units_ojl(k) = 'no output'
        jgrid_ojl(k) = 1
        lgrid_ojl(k) = 1
      enddo
c
      k=0
      kk=0
c
      k=k+1
      JL_M = k
      do n=1,nbas
        kk = kk + 1
        sname_ojl(kk) = 'mo_'//basin(n)(1:3)
      enddo
c
      k=k+1
      JL_PT = k
      do n=1,nbas
        kk = kk + 1
        denom_ojl(kk) = JL_M +n-1
        sname_ojl(kk) = 'temp_'//basin(n)(1:3)
        units_ojl(kk) = 'C'
        lname_ojl(kk) = 'Temperature, '//trim(basin(n))//' Basin'
      enddo
c
      k=k+1
      JL_S = k
      do n=1,nbas
        kk = kk + 1
        denom_ojl(kk) = JL_M +n-1
        sname_ojl(kk) = 'salt_'//basin(n)(1:3)
        units_ojl(kk) = 'psu'
        lname_ojl(kk) = 'Salinity, '//trim(basin(n))//' Basin'
      enddo
c
      k=k+1
      JL_SF = k
      do n=1,nbas
        kk = kk + 1
        sname_ojl(kk) = 'sf_'//basin(n)(1:3)
        units_ojl(kk) = 'Sv'
        lname_ojl(kk) = 'Stream Function, '//trim(basin(n))//' Basin'
        jgrid_ojl(kk) = 2
        lgrid_ojl(kk) = 2
      enddo

#ifdef NEW_IO
c
c Declare the dimensions and metadata of output fields using
c netcdf CDL notation.  The C convention for dimension ordering
c must be used (reversed wrt Fortran).
c

      cdl_oij = ''
      cdl_oij(1:2)(:) = (/
     &     'netcdf xxx { ', 'dimensions:  ' /)
      write(cdl_oij(3),'(a,i3,a)') '   lono = ',im,' ;'
      write(cdl_oij(4),'(a,i3,a)') '   lato = ',jm,' ;'
      write(cdl_oij(5),'(a,i3,a)') '   lono2 = ',im,' ;'
      write(cdl_oij(6),'(a,i3,a)') '   lato2 = ',jm,' ;'
      cdl_oij(7:15)(:) = (/
     &     'variables:                       ',
     &     'float lono(lono) ;               ',
     &     '   lono:units = "degrees_east" ; ',
     &     'float lato(lato) ;               ',
     &     '   lato:units = "degrees_north" ;',
     &     'float lono2(lono2) ;             ',
     &     '   lono2:units = "degrees_east" ;',
     &     'float lato2(lato2) ;             ',
     &     '   lato2:units = "degrees_north" ;'
     &     /)
      kk = count(len_trim(cdl_oij).gt.0)
      do k=1,koij
        if(trim(sname_oij(k)).eq.'unused') cycle
        xstr='lono) ;'
        if(igrid_oij(k).eq.2) xstr='lono2) ;'
        ystr='(lato,'
        if(jgrid_oij(k).eq.2) ystr='(lato2,'
        kk = kk + 1
        cdl_oij(kk) = 'float '//trim(sname_oij(k))//
     &       trim(ystr)//trim(xstr)
        kk = kk + 1
        cdl_oij(kk) = '   '//trim(sname_oij(k))//':long_name = "'//
     &       trim(lname_oij(k))//'" ;'
        kk = kk + 1
        cdl_oij(kk) = '   '//trim(sname_oij(k))//':units = "'//
     &       trim(units_oij(k))//'" ;'
      enddo
      kk = kk + 1
      cdl_oij(kk) = '}'

      cdl_oijl = ''
      cdl_oijl(1:2)(:) = (/
     &     'netcdf xxx { ', 'dimensions:  ' /)
      write(cdl_oijl(3),'(a,i3,a)') '   lono = ',im,' ;'
      write(cdl_oijl(4),'(a,i3,a)') '   lato = ',jm,' ;'
      write(cdl_oijl(5),'(a,i3,a)') '   lono2 = ',im,' ;'
      write(cdl_oijl(6),'(a,i3,a)') '   lato2 = ',jm,' ;'
      write(cdl_oijl(7),'(a,i3,a)') '   zoc = ',lmo,' ;'
      write(cdl_oijl(8),'(a,i3,a)') '   zoce = ',lmo,' ;'
      cdl_oijl(9:21)(:) = (/
     &     'variables:                        ',
     &     'float lono(lono) ;                ',
     &     '   lono:units = "degrees_east" ;  ',
     &     'float lato(lato) ;                ',
     &     '   lato:units = "degrees_north" ; ',
     &     'float lono2(lono2) ;              ',
     &     '   lono2:units = "degrees_east" ; ',
     &     'float lato2(lato2) ;              ',
     &     '   lato2:units = "degrees_north" ;',
     &     'float zoc(zoc) ;                  ',
     &     '   zoc:units = "m" ;              ',
     &     'float zoce(zoce) ;                ',
     &     '   zoce:units = "m" ;             '
     &     /)
      kk = count(len_trim(cdl_oijl).gt.0)
      do k=1,koijl
        if(trim(sname_oijl(k)).eq.'unused') cycle
        if(trim(lname_oijl(k)).eq.'no output') cycle
        xstr='lono) ;'
        if(igrid_oijl(k).eq.2) xstr='lono2) ;'
        ystr='lato,'
        if(jgrid_oijl(k).eq.2) ystr='lato2,'
        zstr='(zoc,'
        if(lgrid_oijl(k).eq.2) zstr='(zoce,'
        kk = kk + 1
        cdl_oijl(kk) = 'float '//trim(sname_oijl(k))//
     &       trim(zstr)//trim(ystr)//trim(xstr)
        kk = kk + 1
        cdl_oijl(kk) = '   '//trim(sname_oijl(k))//':long_name = "'//
     &       trim(lname_oijl(k))//'" ;'
        kk = kk + 1
        cdl_oijl(kk) = '   '//trim(sname_oijl(k))//':units = "'//
     &       trim(units_oijl(k))//'" ;'
      enddo
      kk = kk + 1
      cdl_oijl(kk) = '}'

      cdl_ojl = ''
      cdl_ojl(1:2)(:) = (/
     &     'netcdf xxx { ', 'dimensions:  ' /)
      write(cdl_ojl(3),'(a,i3,a)') '   lato = ',jm,' ;'
      write(cdl_ojl(4),'(a,i3,a)') '   lato2 = ',jm,' ;'
      write(cdl_ojl(5),'(a,i3,a)') '   zoc = ',lmo,' ;'
      write(cdl_ojl(6),'(a,i3,a)') '   zoce = ',lmo,' ;'
      cdl_ojl(7:15)(:) = (/
     &     'variables:                          ',
     &     'float lato(lato) ;                  ',
     &     '   lato:units = "degrees_north" ;   ',
     &     'float lato2(lato2) ;                ',
     &     '   lato2:units = "degrees_north" ;  ',
     &     'float zoc(zoc) ;                    ',
     &     '   zoc:units = "m" ;                ',
     &     'float zoce(zoce) ;                  ',
     &     '   zoce:units = "m" ;               '
     &     /)
      kk = count(len_trim(cdl_ojl).gt.0)
      do k=1,kojl
        if(trim(sname_ojl(k)).eq.'unused') cycle
        if(trim(lname_ojl(k)).eq.'no output') cycle
        ystr='lato) ;'
        if(jgrid_ojl(k).eq.2) ystr='lato2) ;'
        zstr='(zoc,'
        if(lgrid_ojl(k).eq.2) zstr='(zoce,'
        kk = kk + 1
        cdl_ojl(kk) = 'float '//trim(sname_ojl(k))//
     &       trim(zstr)//trim(ystr)
        kk = kk + 1
        cdl_ojl(kk) = '   '//trim(sname_ojl(k))//':long_name = "'//
     &       trim(lname_ojl(k))//'" ;'
        kk = kk + 1
        cdl_ojl(kk) = '   '//trim(sname_ojl(k))//':units = "'//
     &       trim(units_ojl(k))//'" ;'
      enddo
      kk = kk + 1
      cdl_ojl(kk) = '}'

      cdl_olnst = ''
      cdl_olnst(1:2)(:) = (/
     &     'netcdf xxx { ', 'dimensions:  ' /)
      write(cdl_olnst(3),'(a,i3,a)') '   zoc = ',lmo,' ;'
      write(cdl_olnst(4),'(a,i3,a)') '   zoce = ',lmo,' ;'
      write(cdl_olnst(5),'(a,i3,a)') '   nmst = ',nmst,' ;'
      write(cdl_olnst(6),'(a,i3,a)') '   strait_strlen = ',
     &     len(name_st(1)),' ;'
      cdl_olnst(7:13)(:) = (/
     &     'variables:                            ',
     &     'float zoc(zoc) ;                      ',
     &     '   zoc:units = "m" ;                  ',
     &     'float zoce(zoce) ;                    ',
     &     '   zoce:units = "m" ;                 ',
     &     'char strait_name(nmst,strait_strlen) ;',
     &     'int lmst(nmst) ;                      '
     &     /)
      kk = count(len_trim(cdl_olnst).gt.0)
      do k=1,kolnst
        if(trim(sname_olnst(k)).eq.'unused') cycle
        if(trim(lname_olnst(k)).eq.'no output') cycle
        zstr='zoc) ;'
        if(lgrid_olnst(k).eq.2) zstr='zoce) ;'
        kk = kk + 1
        cdl_olnst(kk) = 'float '//trim(sname_olnst(k))//
     &       '(nmst,'//trim(zstr)
        kk = kk + 1
        cdl_olnst(kk) = '   '//trim(sname_olnst(k))//':long_name = "'//
     &       trim(lname_olnst(k))//'" ;'
        kk = kk + 1
        cdl_olnst(kk) = '   '//trim(sname_olnst(k))//':units = "'//
     &       trim(units_olnst(k))//'" ;'
      enddo
      kk = kk + 1
      cdl_olnst(kk) = '}'

      cdl_otj = ''
      cdl_otj(1:2)(:) = (/
     &     'netcdf xxx { ', 'dimensions:  ' /)
      write(cdl_otj(3),'(a,i3,a)') '   lato2 = ',jm,' ;'
      cdl_otj(4:6)(:) = (/
     &     'variables:                        ',
     &     'float lato2(lato2) ;              ',
     &     '   lato2:units = "degrees_north" ;'
     &     /)
      kk = count(len_trim(cdl_otj).gt.0)
      do k=1,kotj
        if(trim(sname_otj(k)).eq.'unused') cycle
        kk = kk + 1
        cdl_otj(kk) = 'float '//trim(sname_otj(k))//'(lato2) ;'
        kk = kk + 1
        cdl_otj(kk) = '   '//trim(sname_otj(k))//':long_name = "'//
     &       trim(lname_otj(k))//'" ;'
        kk = kk + 1
        cdl_otj(kk) = '   '//trim(sname_otj(k))//':units = "'//
     &       trim(units_otj(k))//'" ;'
      enddo
      kk = kk + 1
      cdl_otj(kk) = '}'
#endif

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

      USE DOMAIN_DECOMP_1D, only : dist_grid,get,am_i_root

      USE ODIAG

      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(IN) :: grid

      INTEGER :: J_1H, J_0H
      INTEGER :: IER

      CALL GET(grid, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)

      ALLOCATE(        OIJ_loc (IM,J_0H:J_1H,KOIJ), STAT=IER )
      ALLOCATE(       OIJL_loc (IM,J_0H:J_1H,LMO,KOIJL), STAT=IER )
      ALLOCATE(       OIJL_out (IM,J_0H:J_1H,LMO,KOIJL), STAT=IER )
#ifdef TRACERS_OCEAN
      ALLOCATE(      TOIJL_loc (IM,J_0H:J_1H,LMO,KTOIJL,NTM), STAT=IER )
#endif

      if(am_i_root()) then
        ALLOCATE( OIJ (IM,JM,KOIJ), STAT=IER )
        ALLOCATE(OIJL (IM,JM,LMO,KOIJL), STAT=IER )
#ifdef TRACERS_OCEAN
        ALLOCATE(TOIJL (IM,JM,LMO,KTOIJL,NTM), STAT=IER )
#endif
      endif

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


      subroutine set_owtbudg()
!@sum Precomputes area weights for zonal means on budget grid
!auth M. Kelley
      USE DIAG_COM, only : jm_budg, area_of_zone=>dxyp_budg
      USE OCEAN, only : im,jm, owtbudg, imaxj, dxypo, oJ_BUDG
      USE DOMAIN_DECOMP_1D, only :GET
      USE OCEANR_DIM, only : oGRID
      IMPLICIT NONE
      INTEGER :: I,J,J_0,J_1
      
      CALL GET(ogrid, J_STRT=J_0,J_STOP=J_1)

c Note that when ocean quantities are to be averaged into "budget"
c zones whose latitudinal boundaries do not coincide with those of
c the ocean grid, the calculation here makes the weights at a given
c budget index k sum up to a number slightly different
c than unity:  sum(owtbudg, mask=(oj_budg==k)) != 1.
c This ensures that global means of "zonal averages" calculated using
c budget grid areas are identical to those calculated using ocean
c native-grid areas.
      do J=J_0,J_1
        owtbudg(:,j)=dxypo(j)/area_of_zone(oj_budg(1,j))
      enddo

c compensate for polar loops in conserv_ODIAG not going from 1 to im
      if(ogrid%have_south_pole) owtbudg(:,1) = owtbudg(:,1)*im
      if(ogrid%have_north_pole) owtbudg(:,jm) = owtbudg(:,jm)*im
        
      END SUBROUTINE set_owtbudg


      SUBROUTINE SET_OJ_BUDG
!@sum set mapping from ocean domain to budget grid 
!@auth Gavin Schmidt/Denis Gueyffier
      USE OCEAN, only : oJ_BUDG,oJ_0B,oJ_1B,oLAT2D_DG,IMO=>IM
      USE DIAG_COM, only : jm_budg
      USE DOMAIN_DECOMP_1D, only :GET
      USE OCEANR_DIM, only : oGRID
      IMPLICIT NONE
!@var I,J are atm grid point values for the accumulation
      INTEGER :: I,J,J_0,J_1,J_0H,J_1H
 
C**** define atmospheric grid
      CALL GET(ogrid,
     &     J_STRT=J_0,J_STOP=J_1,
     &     J_STRT_HALO=J_0H,J_STOP_HALO=J_1H)

      DO J=J_0H,J_1H
        DO I=1,IMO
           oJ_BUDG(I,J)=NINT(1+(olat2d_dg(I,J)+90)*(JM_BUDG-1)/180.)
        END DO
      END DO

      oJ_0B=MINVAL( oJ_BUDG(1:IMO,J_0:J_1) )
      oJ_1B=MAXVAL( oJ_BUDG(1:IMO,J_0:J_1) )

      END SUBROUTINE SET_OJ_BUDG

