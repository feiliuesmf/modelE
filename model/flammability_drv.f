#include "rundeck_opts.h"

      module flammability_com    
!@sum for routines to calculate flammability potential of surface 
!@+   vegetation. Right now this is exceedingly simple, not really
!@+   warrenting its own codes, but I set it up in the usual 
!@+   com/drv/physics way because likely one day it will have 
!@+   more complicated coupling to the MODELE vegetation, pbl,
!@+   surface, etc. codes. And then to the tracer code to alter
!@+   emissions from fires.
!@auth Greg Faluvegi based on direction from Olga Pechony
!@ver  1.0 (based on Olga's Word document Flammability.doc)

      implicit none
      save

      real*8, allocatable, dimension(:,:) :: flammability,veg_density

      end module flammability_com


      subroutine alloc_flammability(grid)
!@SUM  alllocates arrays whose sizes need to be determined
!@+    at run-time
!@auth Greg Faluvegi
!@ver  1.0
      use domain_decomp, only: dist_grid, get
      use model_com, only: im
      use flammability_com, only: flammability,veg_density
      
      implicit none
      
      type (dist_grid), intent(in) :: grid
      integer :: ier, J_1H, J_0H       

      call get( grid , J_STRT_HALO=J_0H, J_STOP_HALO=J_1H )

      allocate( flammability(IM,J_0H:J_1H) )
      allocate(    veg_density(IM,J_0H:J_1H) )

      return
      end subroutine alloc_flammability


      subroutine init_flammability
!@sum initialize flamability to zero and read the veg density file
!@auth Greg Faluvegi based on direction from Olga Pechony
!@ver  1.0 
      use model_com, only: im,jm
      use flammability_com, only: flammability, veg_density
      use domain_decomp,only: grid, get, am_i_root, write_parallel,
     & unpack_data
      use filemanager, only: openunit, closeunit

      implicit none

      real*8, allocatable,  dimension(:,:) :: temp_glob
      real*4, allocatable,  dimension(:,:) :: temp4_glob
      character*80 :: title
      character*150 :: out_line
      integer :: J_1H, J_0H, iu_data

      call get(grid,J_STRT_HALO=J_0H,J_STOP_HALO=J_1H)

      flammability(:,J_0H:J_1H) = 0.d0

      if(am_i_root( ))then
        allocate( temp_glob(im,jm) ,temp4_glob(im,jm) )
        call openunit('VEG_DENSE',iu_data,.true.,.true.)
        read (iu_data) title,temp4_glob
        call closeunit(iu_data)
        temp_glob(:,:)=dble(temp4_glob(:,:))
      endif
      write(out_line,*) trim(title),' read from VEG_DENSE'
      call write_parallel(trim(out_line))
      call unpack_data( grid, temp_glob, veg_density )
     
      if(am_i_root()) deallocate( temp_glob, temp4_glob )

      return
      end subroutine init_flammability


      subroutine flammability_drv
!@sum driver routine for flammability potential of surface
!@+   vegetation calculation.
!@auth Greg Faluvegi based on direction from Olga Pechony
!@ver  1.0 
      use model_com, only: im, jm, dtsrc, ptop, p
      use domain_decomp,only: grid, get
      use flammability_com, only: flammability, veg_density
      use pblcom, only: tsavg, qsavg
      use fluxes, only: prec
      use constant, only: sday, lhe
      use diag_com, only: ij_flam,aij=>aij_loc

      implicit none

      integer :: J_0S, J_1S, i, j
      logical :: have_south_pole, have_north_pole     
      real*8 :: qsat ! this is a function in UTILDBL.f
                              
      call get(grid, J_STRT_SKP=J_0S, J_STOP_SKP=J_1S,
     &               HAVE_SOUTH_POLE = have_south_pole,
     &               HAVE_NORTH_POLE = have_north_pole)

      if(have_north_pole)flammability(:,JM)=0.d0
      if(have_south_pole)flammability(:,1) =0.d0
      do j=J_0S,J_1S
        do i=1,im
          call calc_flammability(tsavg(i,j),sday*prec(i,j)/dtsrc,
     &    min(1.d0,qsavg(i,j)/qsat(tsavg(i,j),lhe,p(i,j)+ptop)),
     &    veg_density(i,j),flammability(i,j)) 
          aij(i,j,ij_flam)=aij(i,j,ij_flam)+flammability(i,j)
        enddo
      enddo
  
      return
      end subroutine flammability_drv

