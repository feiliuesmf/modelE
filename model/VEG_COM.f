#include "rundeck_opts.h"

      module veg_com
!@sum  GHYCOM contains the areas used by the Ground Hydrology routines
!@auth N. Kiang, I. Aleinov
!@ver  1.0
      use model_com, only : im,jm
      use sle001, only : ngm,imt,nlsn
#ifdef TRACERS_WATER
      use tracer_com, only : ntm
#endif
      implicit none
      save

!---  boundary conditions (read from file)
!@var vdata(:,:,k)  fraction of gridbox of veg.type k=1-11
      real*8, dimension(im,jm,11) :: vdata

!---  prognostic variables (saved to restart file)
!@var Cint Internal foliage CO2 concentration (mol/m3)
      real*8, dimension(im,jm) :: Cint
!@var Qfol Foliage surface mixing ratio (kg/kg)
      real*8, dimension(im,jm) :: Qfol
!@var cnc_ij canopy conductance
      real*8, dimension(im,jm) :: cnc_ij

!---  work arrays (recomputed for each restart)
      real*8, dimension(ngm,im,jm) :: afr
      real*8, dimension(3,im,jm) :: ala,acs,almass !nyk almass=leaf mass
      real*8, dimension(3,8,im,jm) :: alaf     !nyk lai components by vegtype
      real*8, dimension(8,im,jm) :: alaif   !nyk lai by vegtype (not * vdata!)
      real*8, dimension(im,jm) :: afb,avh,aalbveg !nyk aalbveg
      real*8, dimension(im,jm) :: can_w_capacity
      real*8, dimension(im,jm) :: anm,anf ! adf

      end module veg_com


      subroutine io_vegetation(kunit,iaction,ioerr)
!@sum  io_soils reads and writes soil arrays to file
!@auth I. Aleinov
!@ver  1.0
      use model_com, only : ioread,iowrite,lhead,irerun,irsfic,irsficno
      use veg_com, only : Cint, Qfol, cnc_ij
      implicit none

      integer kunit   !@var kunit unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
!@var ioerr 1 (or -1) if there is (or is not) an error in i/o
      integer, intent(inout) :: ioerr
!@var header character string label for individual records
      character*80 :: header, module_header = "vegetation01"

      write(module_header(lhead+1:80),'(a)') 'cint,qfol,cnc_ij'

      select case (iaction)
      case (:iowrite)            ! output to standard restart file
        write (kunit,err=10) module_header,cint,qfol,cnc_ij
      case (ioread:)            ! input from restart file
        read (kunit,err=10) header,cint,qfol,cnc_ij
        if (header(1:lhead).ne.module_header(1:lhead)) then
          print*,"discrepancy in module version ",header,module_header
          go to 10
        end if
      end select

      return
 10   ioerr=1
      return
      end subroutine io_vegetation


