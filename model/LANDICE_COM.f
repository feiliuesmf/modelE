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

      IMPLICIT NONE
      SAVE
!@var SNOWLI snow amount on land ice (kg/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: SNOWLI
!@var TLANDI temperature of each land ice layer (C)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: TLANDI
!@var MDWNIMP downward implicit ice amount accumulator (kg)
!@var EDWNIMP downward implicit energy amount accumulator (J)
!@var FSHGLM = fraction of SH GMELT water; Sum[FSHGLM(:,:)] = 1
!@var FNHGLM = fraction of NH GMELT water; Sum[FNHGLM(:,:)] = 1
      Real*8,Allocatable,Dimension(:,:) ::
     *  MDWNIMP,EDWNIMP, FSHGLM,FNHGLM

#ifdef TRACERS_WATER
!@var TRSNOWLI tracer amount in land ice snow (kg/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: TRSNOWLI
!@var TRLNDI tracer amount in land ice (kg/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: TRLNDI
!@var TDWNIMP downward implicit tracer amount accumulator (kg)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: TRDWNIMP
#endif

#ifdef COUPLE_GLIMMER
!@var glint_greenland Glimmer model of the Greenland ice sheet
      type(glint_params) :: glint_greenland;
#endif

      END MODULE LANDICE_COM

      SUBROUTINE ALLOC_LANDICE_COM(grid)
!@sum  To allocate arrays whose sizes now need to be determined at
!@+    run time
!@auth NCCS (Goddard) Development Team
      USE DOMAIN_DECOMP_ATM, ONLY : DIST_GRID
      USE RESOLUTION, ONLY : IM,LM
      Use LANDICE_COM, Only: SNOWLI,TLANDI, MDWNIMP,EDWNIMP,
     *                       FSHGLM,FNHGLM
#ifdef TRACERS_WATER
      USE LANDICE_COM, ONLY : TRSNOWLI, TRLNDI, TRDWNIMP
      USE TRACER_COM, only : NTM
#ifdef TRACERS_OCEAN
      ! landice_com should really be the owner of *ACC[PDA,PDG]
      USE LANDICE, only : TRACCPDA, TRACCPDG
#endif
#endif
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(IN) :: grid

      INTEGER :: I_1H, I_0H, J_1H, J_0H
      INTEGER :: IER

      I_0H = grid%I_STRT_HALO
      I_1H = grid%I_STOP_HALO
      J_0H = grid%J_STRT_HALO
      J_1H = grid%J_STOP_HALO

      ALLOCATE( SNOWLI(I_0H:I_1H,J_0H:J_1H),
     *          TLANDI(2,I_0H:I_1H,J_0H:J_1H),
     *          MDWNIMP(I_0H:I_1H,J_0H:J_1H),
     *          EDWNIMP(I_0H:I_1H,J_0H:J_1H),
     *          FSHGLM(I_0H:I_1H,J_0H:J_1H),
     *          FNHGLM(I_0H:I_1H,J_0H:J_1H),
     *          STAT=IER)
#ifdef TRACERS_WATER
      ALLOCATE( TRSNOWLI(NTM,I_0H:I_1H,J_0H:J_1H),
     *          TRLNDI  (NTM,I_0H:I_1H,J_0H:J_1H),
     *          TRDWNIMP(NTM,I_0H:I_1H,J_0H:J_1H),
     *          STAT=IER)
#ifdef TRACERS_OCEAN
      ALLOCATE(TRACCPDA(NTM), TRACCPDG(NTM))
      TRACCPDA = 0.; TRACCPDG = 0.
#endif
#endif

      RETURN
      END SUBROUTINE ALLOC_LANDICE_COM

      SUBROUTINE io_landice(kunit,iaction,ioerr)
!@sum  io_landice reads and writes landice variables to file
!@auth Gavin Schmidt
      USE MODEL_COM, only : ioread,iowrite,lhead,irsfic,irsficno,irerun
      USE DOMAIN_DECOMP_ATM, only : grid
      USE DOMAIN_DECOMP_1D, only : GET, AM_I_ROOT
      USE DOMAIN_DECOMP_1D, only : PACK_DATA, UNPACK_DATA, PACK_COLUMN
      USE DOMAIN_DECOMP_1D, only : UNPACK_COLUMN, broadcast,
     *     BACKSPACE_PARALLEL
      USE LANDICE_COM
      USE LANDICE
      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var HEADER Character string label for individual records
      CHARACTER*80 :: HEADER, MODULE_HEADER = "GLAIC03"
!@var SNOWLI_GLOB dummy global array used for esmf i/o
      REAL*8, DIMENSION(IM,JM) :: SNOWLI_GLOB
!@var TLANDI_GLOB dummy global array used for esmf i/o
      REAL*8, DIMENSION(2,IM,JM) :: TLANDI_GLOB
      REAL*8, DIMENSION(IM,JM) :: MDWNIMP_GLOB, EDWNIMP_GLOB
      INTEGER :: J_0H,J_1H
#ifdef TRACERS_WATER
!@var TRHEADER Character string label for individual records
      CHARACTER*80 :: TRHEADER, TRMODULE_HEADER = "TRGLAC02"
!@var TRSNOWLI_GLOB  dummy global arrays used for i/o
      REAL*8, DIMENSION(NTM,IM,JM) :: TRSNOWLI_GLOB
!@var TRLNDI_GLOB  dummy global arrays used for i/o
      REAL*8, DIMENSION(NTM,IM,JM) :: TRLNDI_GLOB
      REAL*8, DIMENSION(NTM,IM,JM) :: TRDWNIMP_GLOB
      write (TRMODULE_HEADER(lhead+1:80)
     *     ,'(a7,i3,a)')'R8 dim(',NTM,',im,jm):TRSNOWLI,TRLNDI,TRDWN'
#ifdef TRACERS_OCEAN
     *     //',TRACC,TRIMP*2'
#endif
#endif

      MODULE_HEADER(lhead+1:80) = 'R8 SNOW(im,jm),T(2,im,jm),MDWN,EDWN'
     *     //',MACC,EACC,IMP*8'

      SELECT CASE (IACTION)

      CASE (:IOWRITE)            ! output to standard restart file
C**** Gather into global arrays
        CALL PACK_DATA(grid,SNOWLI,SNOWLI_GLOB)
        CALL PACK_COLUMN( grid,TLANDI,TLANDI_GLOB)
        CALL PACK_COLUMN( grid,MDWNIMP,MDWNIMP_GLOB)
        CALL PACK_COLUMN( grid,EDWNIMP,EDWNIMP_GLOB)
        IF (AM_I_ROOT())
     &       WRITE (kunit,err=10) MODULE_HEADER,SNOWLI_GLOB,TLANDI_GLOB
     *       ,MDWNIMP_GLOB,EDWNIMP_GLOB,ACCPDA,ACCPDG,EACCPDA,EACCPDG
     *       ,MICBIMP,EICBIMP
#ifdef TRACERS_WATER
C**** Gather into global arrays
          CALL PACK_COLUMN( grid,TRSNOWLI,TRSNOWLI_GLOB )
          CALL PACK_COLUMN( grid,TRLNDI,  TRLNDI_GLOB  )
          CALL PACK_COLUMN( grid,TRDWNIMP,TRDWNIMP_GLOB )
        IF (AM_I_ROOT())
     &         WRITE (kunit,err=10) TRMODULE_HEADER,TRSNOWLI_GLOB
     *         ,TRLNDI_GLOB,TRDWNIMP_GLOB
#ifdef TRACERS_OCEAN
     *         ,TRACCPDA,TRACCPDG  !  ,TRICBIMP
#endif
#endif

      CASE (IOREAD:)            ! input from restart file
        if ( AM_I_ROOT() ) then
          READ (kunit,err=10) HEADER
          CALL BACKSPACE_PARALLEL(kunit)
          IF (HEADER(1:LHEAD).EQ.MODULE_HEADER(1:LHEAD)) THEN
            READ (kunit,err=10) HEADER,SNOWLI_GLOB,TLANDI_GLOB
     *           ,MDWNIMP_GLOB,EDWNIMP_GLOB,ACCPDA,ACCPDG,EACCPDA
     *           ,EACCPDG,MICBIMP,EICBIMP
          ELSEIF (HEADER(1:LHEAD).EQ."GLAIC01" .or. HEADER(1:LHEAD).EQ
     *           ."GLAIC02") THEN
            READ (kunit,err=10) HEADER,SNOWLI_GLOB,TLANDI_GLOB
            MDWNIMP_GLOB=0 ; EDWNIMP_GLOB=0
            ACCPDA=0 ; ACCPDG=0 ; EACCPDA=0 ; EACCPDG=0
            MICBIMP(:) = 0  ;  EICBIMP(:) = 0
          ELSE
            PRINT*,"Discrepancy in module version ",HEADER,MODULE_HEADER
            GO TO 10
          END IF
        end if
C****** Get useful ESMF parameters
        CALL GET( GRID, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H )
C****** Load data into distributed arrays
        CALL UNPACK_DATA( GRID, SNOWLI_GLOB, SNOWLI)
        CALL UNPACK_COLUMN( GRID, TLANDI_GLOB, TLANDI)
        CALL UNPACK_COLUMN( GRID, MDWNIMP_GLOB, MDWNIMP)
        CALL UNPACK_COLUMN( GRID, EDWNIMP_GLOB, EDWNIMP)
        call broadcast(grid,  ACCPDA)
        call broadcast(grid,  ACCPDG)
        call broadcast(grid, EACCPDA)
        call broadcast(grid, EACCPDG)
        call broadcast(grid, MICBIMP)
        call broadcast(grid, EICBIMP)

#ifdef TRACERS_WATER
        SELECT CASE (IACTION)
        CASE (IRERUN,IOREAD,IRSFIC,IRSFICNO)    ! from reruns/restarts
          if ( AM_I_ROOT() ) then
            READ (kunit,err=10) TRHEADER
            CALL BACKSPACE_PARALLEL(kunit)
            IF (TRHEADER(1:LHEAD).EQ.TRMODULE_HEADER(1:LHEAD)) THEN
              READ (kunit,err=10) TRHEADER,TRSNOWLI_GLOB,TRLNDI_GLOB
     *             ,TRDWNIMP_GLOB
#ifdef TRACERS_OCEAN
     *             ,TRACCPDA,TRACCPDG  !  ,TRICBIMP
#endif
            ELSEIF (TRHEADER(1:LHEAD).EQ."TRGLAC01") THEN
              READ (kunit,err=10) TRHEADER,TRSNOWLI_GLOB,TRLNDI_GLOB
            ELSE
              PRINT*,"Discrepancy in module version ",TRHEADER
     *             ,TRMODULE_HEADER
              GO TO 10
            END IF
          end if                !..am_i_root
C********* Load data into distributed arrays
          CALL UNPACK_COLUMN(GRID, TRSNOWLI_GLOB, TRSNOWLI)
          CALL UNPACK_COLUMN(GRID, TRLNDI_GLOB,   TRLNDI)
          CALL UNPACK_COLUMN(GRID, TRDWNIMP_GLOB, TRDWNIMP)
#ifdef TRACERS_OCEAN
          call broadcast(grid, TRACCPDA)
          call broadcast(grid, TRACCPDG)
c          call broadcast(grid, TRICBIMP)
#endif
        END SELECT
#endif

      END SELECT

      RETURN
 10   IOERR=1
      RETURN
      END SUBROUTINE io_landice


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
      call defvar(grid,fid,snowli,'snowli(dist_im,dist_jm)')
      call defvar(grid,fid,tlandi,'tlandi(d2,dist_im,dist_jm)')
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
     &     'trsnowli(ntm,dist_im,dist_jm)')
      call defvar(grid,fid,trlndi,'trlndi(ntm,dist_im,dist_jm)')
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
#endif /* NEW_IO */
