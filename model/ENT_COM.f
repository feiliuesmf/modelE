#include "rundeck_opts.h"

      module ent_com
!@sum  ENT_COM contains the data needed for Dynamic Vegetation Model (ENT)
!@auth I. Aleinov
!@ver  1.0
      use model_com, only : im,jm
      use ghy_com, only : ngm,imt,nlsn
      use ent_mod
      implicit none
      save

!@var entcells structures which keep the internal state of each Ent cell
      type(entcelltype_public), allocatable :: entcells(:,:)

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

      contains

!!! I/O Not MPI - parallelized yet!!

      subroutine ent_read_state( kunit, entcells )
!@sum read ent state from the file
      type(entcelltype_public), intent(out) :: entcells(:,:)
      integer, intent(in) :: kunit
      !---
      !integer, parameter :: MAX_BUFFER=100000 ! need realistic estimate
      real*8, allocatable ::  buffer(:)
      integer ic, jc, i, j, nbuf

      ic = size(entcells,1)
      jc = size(entcells,2)

      !call openunit('ent_state',iu_entstate,.true.,.true.)
      do j=1,jc
        do i=1,ic
          read(kunit) nbuf
          allocate( buffer(nbuf) )
          read(kunit) buffer
          print *, "ent_read_state: i, j ", i, j
          print *, buffer
         ! check length of buffer : if( buffer(1) > MAX_BUFFER ) ??
          if ( buffer(1) > 0.d0 ) then ! the cell is present
            call ent_cell_construct( entcells(i,j) )
            call ent_cell_unpack(buffer, entcells(i,j))
          endif
          deallocate( buffer )

          call ent_cell_print(999,entcells(i,j))

        enddo
      enddo

      end subroutine ent_read_state


      subroutine ent_write_state( kunit, entcells )
!@sum write ent state to the file
      type(entcelltype_public), intent(in) :: entcells(:,:)
      integer, intent(in) :: kunit
      !---
      real*8, pointer :: buffer(:)
      integer ic, jc, i, j
      integer, save :: counter=0

      ic = size(entcells,1)
      jc = size(entcells,2)

      counter = mod(counter+1,2)

      !call openunit('ent_state_new',iu_entstate,.true.,.false.)
      do j=1,jc
        do i=1,ic
          call ent_cell_pack(buffer, entcells(i,j))
          print *, "ent_write_state: i, j ", i, j
          print *, buffer
          call ent_cell_print(990+counter,entcells(i,j))
          write(kunit) size(buffer)
          write(kunit) buffer
          deallocate(buffer)
        enddo
      enddo

      end subroutine ent_write_state

      end module ent_com


      subroutine io_vegetation(kunit,iaction,ioerr)
!@sum  io_soils reads and writes soil arrays to file
!@auth I. Aleinov
!@ver  1.0
      use model_com, only : ioread,iowrite,lhead,irerun,irsfic,irsficno
      use model_com, only : im,jm
      use domain_decomp, only : grid, am_i_root
      use domain_decomp, only : pack_data, unpack_data
      use ent_com, only : Cint, Qfol, cnc_ij, entcells,
     &     ent_read_state,ent_write_state
      use ent_mod
      
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
        call ent_write_state( kunit, entcells )
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
        call  ent_read_state( kunit, entcells )
      end select

      return
 10   ioerr=1
      return
      end subroutine io_vegetation

      SUBROUTINE ALLOC_ENT_COM(grid)
!@sum  To allocate arrays whose sizes now need to be determined at
!@+    run time
!@auth NCCS (Goddard) Development Team
!@ver  1.0
      USE ENT_MOD
      USE ENT_COM
      USE DOMAIN_DECOMP, ONLY : DIST_GRID, GET
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(IN) :: grid

      INTEGER :: J_1H, J_0H
      INTEGER :: IER

C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)

      ALLOCATE(    entcells(im,J_0H:J_1H),
     *         STAT=IER)
      ! initialize ent cells to something meaningful
      !call ent_cell_construct( entcells ) ! moved to init_module_ent
      call ent_cell_nullify( entcells )


      ALLOCATE(    vdata(im,J_0H:J_1H,12),
     *              Cint(im,J_0H:J_1H),
     *              Qfol(im,J_0H:J_1H),
     *            cnc_ij(im,J_0H:J_1H),
     *         STAT=IER)


      END SUBROUTINE ALLOC_ENT_COM

