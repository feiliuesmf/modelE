#include "rundeck_opts.h"
#ifdef TRACERS_ATM_ONLY
#undef TRACERS_ON
#undef TRACERS_WATER
#endif

      MODULE LANDICE_COM
!@sum  LANDICE_COM contains the model arrays for land ice
!@auth Gavin Schmidt
!@ver  2010/10/13
!@cont io_landice
      USE RESOLUTION, only : im,jm
#ifdef TRACERS_WATER
      USE TRACER_COM, only : NTM
#endif
#ifdef COUPLE_GLIMMER
      USE glint_main, only : glint_params
#endif
      use cdl_mod, only : cdl_type
      use mdiag_com, only : sname_strlen,units_strlen,lname_strlen
      IMPLICIT NONE
      SAVE

!@var nhc number of height classes
      integer :: nhc=1
!@fhc fraction of landice area in each height class (static for testing purposes)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: fhc

!@var ELEVHC: surface elevation, per height class (m)
! ZATMO should be kept consistent with this.
! The value of this ONLY MATTERS for grid cells with landice
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:)   :: ELEVHC

!@var SNOWLI snow amount on land ice (kg/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: SNOWLI
!@var TLANDI temperature of each land ice layer (C)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: TLANDI
!@var MDWNIMP downward implicit ice amount accumulator (kg)
!@var EDWNIMP downward implicit energy amount accumulator (J)
!@var FSHGLM = fraction of SH GMELT water; Sum[FSHGLM(:,:)] = 1
!@var FNHGLM = fraction of NH GMELT water; Sum[FNHGLM(:,:)] = 1
      Real*8,Allocatable,Dimension(:,:) ::
     *  MDWNIMP,EDWNIMP, FSHGLM,FNHGLM

#ifdef TRACERS_WATER
!@var TRSNOWLI tracer amount in land ice snow (kg/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: TRSNOWLI
!@var TRLNDI tracer amount in land ice (kg/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: TRLNDI
!@var TDWNIMP downward implicit tracer amount accumulator (kg)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: TRDWNIMP
#endif

#ifdef COUPLE_GLIMMER
!@var glint_greenland Glimmer model of the Greenland ice sheet
      type(glint_params) :: glint_greenland;
#endif

!@param kijhc number of ijhc accumulations
      integer, parameter :: kijhc=11
!@var ijhc accumulations for glacial ice height-classified diagnostics
      real*8, dimension(:,:,:,:), allocatable :: ijhc
!@var scale_ijhc scale factor for ijhc diagnostics
      real*8, dimension(kijhc) :: scale_ijhc
!@var ia_ijhc,denom_ijhc  idacc-numbers,weights for ijhc diagnostics
      integer, dimension(kijhc) :: ia_ijhc,denom_ijhc,lgrid_ijhc
!@var sname_ijhc short names of ijhc diagnostics
      character(len=sname_strlen), dimension(kijhc) :: sname_ijhc
!@var lname_ijhc,units_ijhc descriptions/units of ijhc diagnostics
      character(len=lname_strlen), dimension(kijhc) :: lname_ijhc
      character(len=units_strlen), dimension(kijhc) :: units_ijhc
!@var cdl_ijhc consolidated metadata for ijhc output fields in cdl notation
      type(cdl_type) :: cdl_ijhc,cdl_ijhc_latlon
!@var ijhc_xxx indices for accumulations
      integer ::
     &     ijhc_frac,ijhc_tsurf,
     &     IJHC_PRECLI,IJHC_RUNLI,IJHC_EVAPLI,IJHC_F0LI,IJHC_TSLI,
     &     IJHC_SHDTLI,IJHC_EVHDT,IJHC_TRHDT,IJHC_IMPMLI,IJHC_IMPHLI

      END MODULE LANDICE_COM

      SUBROUTINE ALLOC_LANDICE_COM(grid)
!@sum  To allocate arrays whose sizes now need to be determined at
!@+    run time
!@auth NCCS (Goddard) Development Team
      USE DOMAIN_DECOMP_ATM, ONLY : DIST_GRID
      USE RESOLUTION, ONLY : IM,LM
      Use LANDICE_COM, Only: NHC,FHC,SNOWLI,TLANDI, MDWNIMP,EDWNIMP,
     *                       FSHGLM,FNHGLM, ELEVHC
#ifdef TRACERS_WATER
      USE LANDICE_COM, ONLY : TRSNOWLI, TRLNDI, TRDWNIMP
      USE TRACER_COM, only : NTM
#ifdef TRACERS_OCEAN
      ! landice_com should really be the owner of *ACC[PDA,PDG]
      USE LANDICE, only : TRACCPDA, TRACCPDG
#endif
#endif
      use Dictionary_mod, only : sync_param, get_param
      USE LANDICE_COM, only : KIJHC,IJHC

      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(IN) :: grid

      INTEGER :: I_1H, I_0H, J_1H, J_0H
      INTEGER :: IER

      call sync_param("NHC",NHC)

      I_0H = grid%I_STRT_HALO
      I_1H = grid%I_STOP_HALO
      J_0H = grid%J_STRT_HALO
      J_1H = grid%J_STOP_HALO

      ALLOCATE( FHC(I_0H:I_1H,J_0H:J_1H,NHC),
     *          ELEVHC(I_0H:I_1H,J_0H:J_1H,NHC),
     *          SNOWLI(I_0H:I_1H,J_0H:J_1H,NHC),
     *          TLANDI(2,I_0H:I_1H,J_0H:J_1H,NHC),
     *          MDWNIMP(I_0H:I_1H,J_0H:J_1H),
     *          EDWNIMP(I_0H:I_1H,J_0H:J_1H),
     *          FSHGLM(I_0H:I_1H,J_0H:J_1H),
     *          FNHGLM(I_0H:I_1H,J_0H:J_1H),
     *          STAT=IER)
      fhc(:,:,:) = 1d0/nhc
      elevhc(:,:,:) = 0
#ifdef TRACERS_WATER
      ALLOCATE( TRSNOWLI(NTM,I_0H:I_1H,J_0H:J_1H,NHC),
     *          TRLNDI  (NTM,I_0H:I_1H,J_0H:J_1H,NHC),
     *          TRDWNIMP(NTM,I_0H:I_1H,J_0H:J_1H),
     *          STAT=IER)
#ifdef TRACERS_OCEAN
      ALLOCATE(TRACCPDA(NTM), TRACCPDG(NTM))
      TRACCPDA = 0.; TRACCPDG = 0.
#endif
#endif

      ALLOCATE(IJHC(I_0H:I_1H,J_0H:J_1H,NHC,KIJHC))

      RETURN
      END SUBROUTINE ALLOC_LANDICE_COM

c      SUBROUTINE io_landice(kunit,iaction,ioerr)
c!@sum  io_landice reads and writes landice variables to file
c!@auth Gavin Schmidt
c      USE MODEL_COM, only : ioread,iowrite,lhead,irsfic,irsficno,irerun
c      USE DOMAIN_DECOMP_ATM, only : grid
c      USE DOMAIN_DECOMP_1D, only : GET, AM_I_ROOT
c      USE DOMAIN_DECOMP_1D, only : PACK_DATA, UNPACK_DATA, PACK_COLUMN
c      USE DOMAIN_DECOMP_1D, only : UNPACK_COLUMN, broadcast,
c     *     BACKSPACE_PARALLEL
c      USE LANDICE_COM
c      USE LANDICE
c      IMPLICIT NONE
c
c      INTEGER kunit   !@var kunit unit number of read/write
c      INTEGER iaction !@var iaction flag for reading or writing to file
c!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
c      INTEGER, INTENT(INOUT) :: IOERR
c!@var HEADER Character string label for individual records
c      CHARACTER*80 :: HEADER, MODULE_HEADER = "GLAIC03"
c!@var SNOWLI_GLOB dummy global array used for esmf i/o
c      REAL*8, DIMENSION(IM,JM) :: SNOWLI_GLOB
c!@var TLANDI_GLOB dummy global array used for esmf i/o
c      REAL*8, DIMENSION(2,IM,JM) :: TLANDI_GLOB
c      REAL*8, DIMENSION(IM,JM) :: MDWNIMP_GLOB, EDWNIMP_GLOB
c      INTEGER :: J_0H,J_1H
c#ifdef TRACERS_WATER
c!@var TRHEADER Character string label for individual records
c      CHARACTER*80 :: TRHEADER, TRMODULE_HEADER = "TRGLAC02"
c!@var TRSNOWLI_GLOB  dummy global arrays used for i/o
c      REAL*8, DIMENSION(NTM,IM,JM) :: TRSNOWLI_GLOB
c!@var TRLNDI_GLOB  dummy global arrays used for i/o
c      REAL*8, DIMENSION(NTM,IM,JM) :: TRLNDI_GLOB
c      REAL*8, DIMENSION(NTM,IM,JM) :: TRDWNIMP_GLOB
c      write (TRMODULE_HEADER(lhead+1:80)
c     *     ,'(a7,i3,a)')'R8 dim(',NTM,',im,jm):TRSNOWLI,TRLNDI,TRDWN'
c#ifdef TRACERS_OCEAN
c     *     //',TRACC,TRIMP*2'
c#endif
c#endif
c
c      MODULE_HEADER(lhead+1:80) = 'R8 SNOW(im,jm),T(2,im,jm),MDWN,EDWN'
c     *     //',MACC,EACC,IMP*8'
c
c      SELECT CASE (IACTION)
c
c      CASE (:IOWRITE)            ! output to standard restart file
cC**** Gather into global arrays
c        CALL PACK_DATA(grid,SNOWLI,SNOWLI_GLOB)
c        CALL PACK_COLUMN( grid,TLANDI,TLANDI_GLOB)
c        CALL PACK_COLUMN( grid,MDWNIMP,MDWNIMP_GLOB)
c        CALL PACK_COLUMN( grid,EDWNIMP,EDWNIMP_GLOB)
c        IF (AM_I_ROOT())
c     &       WRITE (kunit,err=10) MODULE_HEADER,SNOWLI_GLOB,TLANDI_GLOB
c     *       ,MDWNIMP_GLOB,EDWNIMP_GLOB,ACCPDA,ACCPDG,EACCPDA,EACCPDG
c     *       ,MICBIMP,EICBIMP
c#ifdef TRACERS_WATER
cC**** Gather into global arrays
c          CALL PACK_COLUMN( grid,TRSNOWLI,TRSNOWLI_GLOB )
c          CALL PACK_COLUMN( grid,TRLNDI,  TRLNDI_GLOB  )
c          CALL PACK_COLUMN( grid,TRDWNIMP,TRDWNIMP_GLOB )
c        IF (AM_I_ROOT())
c     &         WRITE (kunit,err=10) TRMODULE_HEADER,TRSNOWLI_GLOB
c     *         ,TRLNDI_GLOB,TRDWNIMP_GLOB
c#ifdef TRACERS_OCEAN
c     *         ,TRACCPDA,TRACCPDG  !  ,TRICBIMP
c#endif
c#endif
c
c      CASE (IOREAD:)            ! input from restart file
c        if ( AM_I_ROOT() ) then
c          READ (kunit,err=10) HEADER
c          CALL BACKSPACE_PARALLEL(kunit)
c          IF (HEADER(1:LHEAD).EQ.MODULE_HEADER(1:LHEAD)) THEN
c            READ (kunit,err=10) HEADER,SNOWLI_GLOB,TLANDI_GLOB
c     *           ,MDWNIMP_GLOB,EDWNIMP_GLOB,ACCPDA,ACCPDG,EACCPDA
c     *           ,EACCPDG,MICBIMP,EICBIMP
c          ELSEIF (HEADER(1:LHEAD).EQ."GLAIC01" .or. HEADER(1:LHEAD).EQ
c     *           ."GLAIC02") THEN
c            READ (kunit,err=10) HEADER,SNOWLI_GLOB,TLANDI_GLOB
c            MDWNIMP_GLOB=0 ; EDWNIMP_GLOB=0
c            ACCPDA=0 ; ACCPDG=0 ; EACCPDA=0 ; EACCPDG=0
c            MICBIMP(:) = 0  ;  EICBIMP(:) = 0
c          ELSE
c            PRINT*,"Discrepancy in module version ",HEADER,MODULE_HEADER
c            GO TO 10
c          END IF
c        end if
cC****** Get useful ESMF parameters
c        CALL GET( GRID, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H )
cC****** Load data into distributed arrays
c        CALL UNPACK_DATA( GRID, SNOWLI_GLOB, SNOWLI)
c        CALL UNPACK_COLUMN( GRID, TLANDI_GLOB, TLANDI)
c        CALL UNPACK_COLUMN( GRID, MDWNIMP_GLOB, MDWNIMP)
c        CALL UNPACK_COLUMN( GRID, EDWNIMP_GLOB, EDWNIMP)
c        call broadcast(grid,  ACCPDA)
c        call broadcast(grid,  ACCPDG)
c        call broadcast(grid, EACCPDA)
c        call broadcast(grid, EACCPDG)
c        call broadcast(grid, MICBIMP)
c        call broadcast(grid, EICBIMP)
c
c#ifdef TRACERS_WATER
c        SELECT CASE (IACTION)
c        CASE (IRERUN,IOREAD,IRSFIC,IRSFICNO)    ! from reruns/restarts
c          if ( AM_I_ROOT() ) then
c            READ (kunit,err=10) TRHEADER
c            CALL BACKSPACE_PARALLEL(kunit)
c            IF (TRHEADER(1:LHEAD).EQ.TRMODULE_HEADER(1:LHEAD)) THEN
c              READ (kunit,err=10) TRHEADER,TRSNOWLI_GLOB,TRLNDI_GLOB
c     *             ,TRDWNIMP_GLOB
c#ifdef TRACERS_OCEAN
c     *             ,TRACCPDA,TRACCPDG  !  ,TRICBIMP
c#endif
c            ELSEIF (TRHEADER(1:LHEAD).EQ."TRGLAC01") THEN
c              READ (kunit,err=10) TRHEADER,TRSNOWLI_GLOB,TRLNDI_GLOB
c            ELSE
c              PRINT*,"Discrepancy in module version ",TRHEADER
c     *             ,TRMODULE_HEADER
c              GO TO 10
c            END IF
c          end if                !..am_i_root
cC********* Load data into distributed arrays
c          CALL UNPACK_COLUMN(GRID, TRSNOWLI_GLOB, TRSNOWLI)
c          CALL UNPACK_COLUMN(GRID, TRLNDI_GLOB,   TRLNDI)
c          CALL UNPACK_COLUMN(GRID, TRDWNIMP_GLOB, TRDWNIMP)
c#ifdef TRACERS_OCEAN
c          call broadcast(grid, TRACCPDA)
c          call broadcast(grid, TRACCPDG)
cc          call broadcast(grid, TRICBIMP)
c#endif
c        END SELECT
c#endif
c
c      END SELECT
c
c      RETURN
c 10   IOERR=1
c      RETURN
c      END SUBROUTINE io_landice


#ifdef NEW_IO
      subroutine read_landice_ic
!@sum   read_landice_ic read land ice initial conditions file.
      use model_com, only : ioread
      use domain_decomp_atm, only : grid
      use pario, only : par_open,par_close
      implicit none
      integer :: fid
      fid = par_open(grid,'GIC','read')
      call new_io_landice(fid,ioread)
      call par_close(grid,fid)
      return
      end subroutine read_landice_ic

      subroutine def_rsf_landice(fid)
!@sum  def_rsf_landice defines landice array structure in restart files
!@auth M. Kelley
!@ver  beta
      use landice_com
      use landice
      use domain_decomp_atm, only : grid
      use pario, only : defvar
      use conserv_diags
      implicit none
      integer fid   !@var fid file id
      call defvar(grid,fid,fhc,'fhc(dist_im,dist_jm,nhc)')
      call defvar(grid,fid,elevhc,'elevhc(dist_im,dist_jm,nhc)')
      call defvar(grid,fid,snowli,'snowli(dist_im,dist_jm,nhc)')
      call defvar(grid,fid,tlandi,'tlandi(d2,dist_im,dist_jm,nhc)')
      call defvar(grid,fid,mdwnimp,'mdwnimp(dist_im,dist_jm)')
      call defvar(grid,fid,edwnimp,'edwnimp(dist_im,dist_jm)')
      call defvar(grid,fid,accpda,'accpda')
      call defvar(grid,fid,eaccpda,'eaccpda')
      call defvar(grid,fid,accpdg,'accpdg')
      call defvar(grid,fid,eaccpdg,'eaccpdg')
      call defvar(grid,fid,micbimp,'micbimp(two)')
      call defvar(grid,fid,eicbimp,'eicbimp(two)')
#ifdef TRACERS_WATER
      call defvar(grid,fid,trsnowli,
     &     'trsnowli(ntm,dist_im,dist_jm,nhc)')
      call defvar(grid,fid,trlndi,'trlndi(ntm,dist_im,dist_jm,nhc)')
      call defvar(grid,fid,trdwnimp,
     &     'trdwnimp(ntm,dist_im,dist_jm)')
#ifdef TRACERS_OCEAN
      call defvar(grid,fid,traccpda,'traccpda(ntm)')
      call defvar(grid,fid,traccpdg,'traccpdg(ntm)')
c      call defvar(grid,fid,tricbimp,'tricbimp(ntm,two)')
#endif
#endif
      call declare_conserv_diags( grid, fid, 'wlani(dist_im,dist_jm)' )
      call declare_conserv_diags( grid, fid, 'elani(dist_im,dist_jm)' )
      call declare_conserv_diags( grid, fid, 'wiceb(dist_im,dist_jm)' )
      call declare_conserv_diags( grid, fid, 'eiceb(dist_im,dist_jm)' )
      return
      end subroutine def_rsf_landice

      subroutine new_io_landice(fid,iaction)
!@sum  new_io_landice read/write landice arrays from/to restart files
!@auth M. Kelley
!@ver  beta new_ prefix avoids name clash with the default version
      use model_com, only : ioread,iowrite
      use domain_decomp_atm, only : grid
      use pario, only : write_dist_data,read_dist_data,
     &     write_data,read_data
      use landice_com
      use landice
      use conserv_diags
      implicit none
      integer fid   !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
      external conserv_MLI, conserv_MICB, conserv_HLI, conserv_HICB
      select case (iaction)
      case (iowrite)            ! output to restart file
        call write_dist_data(grid,fid,'fhc',fhc)
        call write_dist_data(grid,fid,'elevhc',elevhc)
        call write_dist_data(grid,fid,'snowli',snowli)
        call write_dist_data(grid,fid,'tlandi',tlandi,jdim=3)
        call write_dist_data(grid,fid,'mdwnimp',mdwnimp)
        call write_dist_data(grid,fid,'edwnimp',edwnimp)
        call write_data(grid,fid,'accpda',accpda)
        call write_data(grid,fid,'eaccpda',eaccpda)
        call write_data(grid,fid,'accpdg',accpdg)
        call write_data(grid,fid,'eaccpdg',eaccpdg)
        call write_data(grid,fid,'micbimp',micbimp)
        call write_data(grid,fid,'eicbimp',eicbimp)
#ifdef TRACERS_WATER
        call write_dist_data(grid,fid,'trsnowli',trsnowli,jdim=3)
        call write_dist_data(grid,fid,'trlndi',trlndi,jdim=3)
        call write_dist_data(grid,fid,'trdwnimp',trdwnimp,jdim=3)
#ifdef TRACERS_OCEAN
        call write_data(grid,fid,'traccpda',traccpda)
        call write_data(grid,fid,'traccpdg',traccpdg)
c        call write_data(grid,fid,'tricbimp',tricbimp)
#endif
#endif
        call dump_conserv_diags( grid, fid, 'wlani', conserv_MLI )
        call dump_conserv_diags( grid, fid, 'elani', conserv_HLI )
        call dump_conserv_diags( grid, fid, 'wiceb', conserv_MICB )
        call dump_conserv_diags( grid, fid, 'eiceb', conserv_HICB )
      case (ioread)            ! input from restart file
        call read_dist_data(grid,fid,'fhc',fhc)
        call read_dist_data(grid,fid,'elevhc',elevhc)
        call read_dist_data(grid,fid,'snowli',snowli)
        call read_dist_data(grid,fid,'tlandi',tlandi,jdim=3)
c set some defaults for quantities which may not be in the
c restart file
        mdwnimp(:,:) = 0.; edwnimp(:,:) = 0.
        accpda = 0.; eaccpda = 0.
        accpdg = 0.; eaccpdg = 0.
        MICBIMP(:) = 0  ;  EICBIMP(:) = 0
        call read_dist_data(grid,fid,'mdwnimp',mdwnimp)
        call read_dist_data(grid,fid,'edwnimp',edwnimp)
        call read_data(grid,fid,'accpda',accpda,bcast_all=.true.)
        call read_data(grid,fid,'eaccpda',eaccpda,bcast_all=.true.)
        call read_data(grid,fid,'accpdg',accpdg,bcast_all=.true.)
        call read_data(grid,fid,'eaccpdg',eaccpdg,bcast_all=.true.)
        call read_data(grid,fid,'micbimp',micbimp,bcast_all=.true.)
        call read_data(grid,fid,'eicbimp',eicbimp,bcast_all=.true.)
#ifdef TRACERS_WATER
        call read_dist_data(grid,fid,'trsnowli',trsnowli,jdim=3)
        call read_dist_data(grid,fid,'trlndi',trlndi,jdim=3)
        call read_dist_data(grid,fid,'trdwnimp',trdwnimp,jdim=3)
#ifdef TRACERS_OCEAN
        call read_data(grid,fid,'traccpda',traccpda,
     &       bcast_all=.true.)
        call read_data(grid,fid,'traccpdg',traccpdg,
     &       bcast_all=.true.)
c        call read_data(grid,fid,'tricbimp',tricbimp,bcast_all=.true.)
#endif
#endif
      end select
      return
      end subroutine new_io_landice

      subroutine def_rsf_glaacc(fid,r4_on_disk)
!@sum  def_rsf_glaacc defines accumulation array structure in restart/acc files
!@auth M. Kelley
!@ver  beta
      use landice_com, only : ijhc
      use domain_decomp_atm, only : grid
      use pario, only : defvar
      implicit none
      integer fid   !@var fid file id
      logical :: r4_on_disk  !@var r4_on_disk if true, real*8 stored as real*4

      call defvar(grid,fid,ijhc,'ijhc(dist_im,dist_jm,nhc,kijhc)',
     &     r4_on_disk=r4_on_disk)

      return
      end subroutine def_rsf_glaacc

      subroutine new_io_glaacc(fid,iaction)
!@sum  new_io_glaacc read/write accumulation arrays from/to restart/acc files
!@auth M. Kelley
!@ver  beta new_ prefix avoids name clash with the default version
      use resolution, only : im,jm
      use model_com, only : ioread,iowrite,iowrite_single,idacc
      use landice_com, only : ijhc,ijhc_frac,ia_ijhc
      use domain_decomp_atm, only : grid
      use domain_decomp_1d, only : hasNorthPole, hasSouthPole
      use pario, only : write_dist_data,read_dist_data
      implicit none
      integer fid   !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
      integer :: k,l
      select case (iaction)
      case (iowrite,iowrite_single) ! output to restart or acc file
        if(iaction.eq.iowrite_single) then ! pole fills needed for acc-files
          do l=1,size(ijhc,4)
            do k=1,size(ijhc,3)
              if(hasSouthPole(grid)) then
                ijhc(2:im, 1,k,l) = ijhc(1, 1,k,l)
              endif
              if(hasNorthPole(grid)) then
                ijhc(2:im,jm,k,l) = ijhc(1,jm,k,l)
              endif
            enddo
          enddo
          do l=1,size(ijhc,4) ! mult by area fraction for masking purposes
            if(l == ijhc_frac) cycle
            ijhc(:,:,:,l) = ijhc(:,:,:,l)*
     &           (ijhc(:,:,:,ijhc_frac)/idacc(ia_ijhc(ijhc_frac)))
          enddo
        endif
        call write_dist_data(grid,fid,'ijhc',ijhc)
      case (ioread)            ! input from restart or acc file
        call read_dist_data(grid,fid,'ijhc',ijhc)
      end select
      return
      end subroutine new_io_glaacc

      subroutine def_meta_glaacc(fid)
!@sum  def_meta_glaacc defines metadata in acc files
!@auth M. Kelley
!@ver  beta
      use landice_com, only :
     &     ia_ijhc,scale_ijhc,denom_ijhc,sname_ijhc
     &     ,cdl_ijhc,cdl_ijhc_latlon
      use domain_decomp_atm, only : grid
      use pario, only : defvar,write_attr
      use cdl_mod, only : defvar_cdl
      implicit none
      integer :: fid         !@var fid file id

      call write_attr(grid,fid,'ijhc','reduction','sum')
      call write_attr(grid,fid,'ijhc','split_dim',4)
      call defvar(grid,fid,ia_ijhc,'ia_ijhc(kijhc)')
      call defvar(grid,fid,scale_ijhc,'scale_ijhc(kijhc)')
      call defvar(grid,fid,denom_ijhc,'denom_ijhc(kijhc)')
      call defvar(grid,fid,sname_ijhc,'sname_ijhc(sname_strlen,kijhc)')
      call defvar_cdl(grid,fid,cdl_ijhc,
     &     'cdl_ijhc(cdl_strlen,kcdl_ijhc)')
#ifdef CUBED_SPHERE
      call defvar_cdl(grid,fid,cdl_ijhc_latlon,
     &     'cdl_ijhc_latlon(cdl_strlen,kcdl_ijhc_latlon)')
#endif

      return
      end subroutine def_meta_glaacc

      subroutine write_meta_glaacc(fid)
!@sum  write_meta_glaacc writes metadata to acc files
!@auth M. Kelley
!@ver  beta
      use landice_com, only :
     &     ia_ijhc,scale_ijhc,denom_ijhc,sname_ijhc
     &     ,cdl_ijhc,cdl_ijhc_latlon
      use domain_decomp_atm, only : grid
      use pario, only : defvar,write_data
      use cdl_mod, only : write_cdl
      implicit none
      integer :: fid         !@var fid file id

      call write_data(grid,fid,'ia_ijhc',ia_ijhc)
      call write_data(grid,fid,'scale_ijhc',scale_ijhc)
      call write_data(grid,fid,'denom_ijhc',denom_ijhc)
      call write_data(grid,fid,'sname_ijhc',sname_ijhc)
      call write_cdl(grid,fid,'cdl_ijhc',cdl_ijhc)
#ifdef CUBED_SPHERE
      call write_cdl(grid,fid,'cdl_ijhc_latlon',cdl_ijhc_latlon)
#endif

      return
      end subroutine write_meta_glaacc

#endif /* NEW_IO */
