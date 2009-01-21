#include "rundeck_opts.h"

      module veg_com
!@sum  GHY_COM contains the areas used by the Ground Hydrology routines
!@auth N. Kiang, I. Aleinov
!@ver  1.0
      use ghy_com, only : ngm,imt,nlsn
#ifdef TRACERS_WATER
      use tracer_com, only : ntm
#endif
      implicit none
      save

!---  boundary conditions (read from file)
!@var vdata(:,:,k)  fraction of gridbox of veg.type k=1-12
      real*8, ALLOCATABLE, dimension(:,:,:) :: vdata

!---  prognostic variables (saved to restart file)
!@var Cint Internal foliage CO2 concentration (mol/m3)
      real*8, ALLOCATABLE, dimension(:,:) :: Cint
!@var Qfol Foliage surface mixing ratio (kg/kg)
      real*8, ALLOCATABLE, dimension(:,:) :: Qfol
!@var cnc_ij canopy conductance
      real*8, ALLOCATABLE, dimension(:,:) :: cnc_ij

!---  work arrays (recomputed for each restart)
      real*8, ALLOCATABLE, dimension(:,:,:) :: afr
      real*8, ALLOCATABLE, dimension(:,:,:) :: ala,acs,almass !nyk almass=leaf mass
      real*8, ALLOCATABLE, dimension(:,:,:,:) :: alaf     !nyk lai components by vegtype
      real*8, ALLOCATABLE, dimension(:,:,:) :: alaif   !nyk lai by vegtype (not * vdata!)
      real*8, ALLOCATABLE, dimension(:,:) :: afb,avh,aalbveg !nyk aalbveg
      real*8, ALLOCATABLE, dimension(:,:) :: can_w_capacity
      real*8, ALLOCATABLE, dimension(:,:) :: anm,anf ! adf

      end module veg_com


      subroutine io_vegetation(kunit,iaction,ioerr)
!@sum  io_soils reads and writes soil arrays to file
!@auth I. Aleinov
!@ver  1.0
      use model_com, only : ioread,iowrite,lhead,irerun,irsfic,irsficno
      use model_com, only : im,jm
      use domain_decomp_atm, only : grid, am_i_root
      use domain_decomp_atm, only : pack_data, unpack_data
      use veg_com, only : Cint, Qfol, cnc_ij
      implicit none

      integer kunit   !@var kunit unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
!@var ioerr 1 (or -1) if there is (or is not) an error in i/o
      integer, intent(inout) :: ioerr
!@var header character string label for individual records
      character*80 :: header, module_header = "vegetation01"
!@var cint_glob work array for parallel_io
!@var qfol_glob work array for parallel_io
!@var cnc_ij_glob work array for parallel_io
      real*8, dimension(im,jm) :: cint_glob, qfol_glob, cnc_ij_glob

      write(module_header(lhead+1:80),'(a)') 'cint,qfol,cnc_ij'

      select case (iaction)
      case (:iowrite)            ! output to standard restart file
        call pack_data(grid, cint, cint_glob)
        call pack_data(grid, qfol, qfol_glob)
        call pack_data(grid, cnc_ij, cnc_ij_glob)
        if (am_i_root())
     &    write (kunit,err=10) module_header,cint_glob,qfol_glob,
     &                         cnc_ij_glob
      case (ioread:)            ! input from restart file
        if ( AM_I_ROOT() ) then
          read (kunit,err=10) header, cint_glob, qfol_glob, cnc_ij_glob
          if (header(1:lhead).ne.module_header(1:lhead)) then
            print*,"discrepancy in module version ",header,module_header
            go to 10
          end if
        end if
        call unpack_data(grid, cint_glob,   cint)
        call unpack_data(grid, qfol_glob,   qfol)
        call unpack_data(grid, cnc_ij_glob, cnc_ij)
      end select

      return
 10   ioerr=1
      return
      end subroutine io_vegetation

      SUBROUTINE ALLOC_VEG_COM(grid)
!@sum  To allocate arrays whose sizes now need to be determined at
!@+    run time
!@auth NCCS (Goddard) Development Team
!@ver  1.0
      USE VEG_COM
      USE DOMAIN_DECOMP_ATM, ONLY : DIST_GRID, GET
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(IN) :: grid

      INTEGER :: I_1H, I_0H, J_1H, J_0H
      INTEGER :: IER

C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)
      I_0H = grid%I_STRT_HALO
      I_1H = grid%I_STOP_HALO

      ALLOCATE(    vdata(I_0H:I_1H,J_0H:J_1H,12),
     *              Cint(I_0H:I_1H,J_0H:J_1H),
     *              Qfol(I_0H:I_1H,J_0H:J_1H),
     *            cnc_ij(I_0H:I_1H,J_0H:J_1H),
     *               afr(ngm,I_0H:I_1H,J_0H:J_1H),
     *         STAT=IER)

      ALLOCATE(    ala(3,  I_0H:I_1H,J_0H:J_1H),
     *             acs(3,  I_0H:I_1H,J_0H:J_1H),
     *          almass(3,  I_0H:I_1H,J_0H:J_1H),
     *            alaf(3,11,I_0H:I_1H,J_0H:J_1H),
     *           alaif(11,  I_0H:I_1H,J_0H:J_1H),
     *         STAT=IER)

      ALLOCATE(   afb(I_0H:I_1H,J_0H:J_1H),
     *            avh(I_0H:I_1H,J_0H:J_1H),
     *            aalbveg(I_0H:I_1H,J_0H:J_1H),
     *            can_w_capacity(I_0H:I_1H,J_0H:J_1H),
     *            anm(I_0H:I_1H,J_0H:J_1H),
     *            anf(I_0H:I_1H,J_0H:J_1H),
     *         STAT=IER)

      END SUBROUTINE ALLOC_VEG_COM

