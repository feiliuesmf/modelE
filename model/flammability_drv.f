#include "rundeck_opts.h"
      module flammability_com    
!@sum for routines to calculate flammability potential of surface 
!@+   vegetation. Optionally also altering tracer biomass sources.
!@auth Greg Faluvegi based on direction from Olga Pechony
!@ver  1.0 (based on Olga's Word document Flammability.doc)

      implicit none
      save

! main variables:
!@var flammability unitless flammability coefficient
!@var base_flam unitless flammability coefficient from a run with 
!@+   baseline climate (e.g. GFED 1997-2002 period), monthly averages
!@var vegetation density read from file for flammability calculation
      real*8, allocatable, dimension(:,:) :: flammability,veg_density
#ifdef DYNAMIC_BIOMASS_BURNING
!@param mfcc MODIS fire count calibration (goes with the input files)
!@+ units are fires/day when multiplied by the unitless flammability
      real*8, parameter :: mfcc=680.d0
!@var epfc emission per fire count, defined for ntm tracers for now
!@+ generally in kg/m2/s/fire
      real*8, allocatable, dimension(:,:,:) :: epfc 
#endif
      real*8, parameter :: missing=-1.d30

! rest is for the running average:
!@param maxHR_prec maximum number of sub-daily accumulations
!@param nday_prec number of days in running average for prec
!@var DRAfl daily running average of model prec for flammability
!@var ravg_prec period running average of model prec for flammability
!@var PRSfl period running sum of model prec for flammability
!@var HRAfl hourly running average of model prec for flammability
!@var iHfl "hourly" index for averages of model prec
!@var iDfl "daily"  index for averages of model prec
!@var i0fl ponter to current index in running sum of model prec
!@var first_prec whether in the first model averaging per. for prec
      ! next line: 1800=dtsrc, 1=#calls per step, but since those are
      ! not parameters, I don't know how to soft-code this. There is a
      ! failsafe in flammability_drv, however:
      integer, parameter :: maxHR_prec=24*1*3600/1800, nday_prec=30
      real*8, allocatable, dimension(:,:,:):: DRAfl
      real*8, allocatable, dimension(:,:)  :: ravg_prec,PRSfl,iHfl,iDfl,
     &                                        i0fl,first_prec
      real*8, allocatable, dimension(:,:,:):: HRAfl
!@var raP_acc accumulate running avg precip for SUBDD output
      real*8, allocatable, dimension(:,:):: raP_acc

      end module flammability_com


      subroutine alloc_flammability(grid)
!@SUM  alllocates arrays whose sizes need to be determined
!@+    at run-time
!@auth Greg Faluvegi
!@ver  1.0
      use domain_decomp_atm, only: dist_grid, get
      use model_com, only: im
      use flammability_com, only: flammability,veg_density,
     & first_prec,iHfl,iDfl,i0fl,DRAfl,ravg_prec,PRSfl,HRAfl,
     & nday_prec,maxHR_prec,raP_acc
#ifdef DYNAMIC_BIOMASS_BURNING
     & ,epfc
      use tracer_com, only: ntm
#endif 
 
      implicit none
      
      type (dist_grid), intent(in) :: grid
      integer :: ier, J_1H, J_0H, I_1H, I_0H

      call get( grid , J_STRT_HALO=J_0H, J_STOP_HALO=J_1H )
      I_0H=GRID%I_STRT_HALO
      I_1H=GRID%I_STOP_HALO 

      allocate( flammability(I_0H:I_1H,J_0H:J_1H) )
      allocate( veg_density (I_0H:I_1H,J_0H:J_1H) )
      allocate( first_prec  (I_0H:I_1H,J_0H:J_1H) )
      allocate( iHfl        (I_0H:I_1H,J_0H:J_1H) )
      allocate( iDfl        (I_0H:I_1H,J_0H:J_1H) )
      allocate( i0fl        (I_0H:I_1H,J_0H:J_1H) )
      allocate( DRAfl       (I_0H:I_1H,J_0H:J_1H,nday_prec) )
      allocate( ravg_prec   (I_0H:I_1H,J_0H:J_1H) )
      allocate( PRSfl       (I_0H:I_1H,J_0H:J_1H) )
      allocate( raP_acc     (I_0H:I_1H,J_0H:J_1H) )
      allocate( HRAfl       (I_0H:I_1H,J_0H:J_1H,maxHR_prec) )
#ifdef DYNAMIC_BIOMASS_BURNING
      allocate( epfc        (I_0H:I_1H,J_0H:J_1H,ntm) )
#endif

      return
      end subroutine alloc_flammability


      subroutine init_flammability
!@sum initialize flamability, veg density, etc. for fire model
!@auth Greg Faluvegi based on direction from Olga Pechony
!@ver  1.0 
      use model_com, only: im,jm,Itime,ItimeI
      use flammability_com, only: flammability, veg_density, first_prec
     & ,missing
#ifdef DYNAMIC_BIOMASS_BURNING
     & ,epfc
      use tracer_com, only: ntm,trname,do_fire
#endif
      use domain_decomp_atm,only: grid, get, am_i_root, readt_parallel
      use filemanager, only: openunit, closeunit, nameunit

      implicit none
      character*80 :: title,fname

      integer :: I_1H, I_0H, J_1H, J_0H, iu_data, n

      call get(grid,J_STRT_HALO=J_0H,J_STOP_HALO=J_1H)
      call get(grid,I_STRT_HALO=I_0H,I_STOP_HALO=I_1H)

      call openunit('VEG_DENSE',iu_data,.true.,.true.)
      call readt_parallel(grid,iu_data,nameunit(iu_data),veg_density,1)
      call closeunit(iu_data)

      if(Itime==ItimeI)then
        first_prec(:,:)=1.d0
        flammability(:,:)=missing
      endif

#ifdef DYNAMIC_BIOMASS_BURNING
      epfc(:,:,:)=0.d0 ! by default no fire emissions for tracer
      do n=1,ntm
        if(do_fire(n)) then
          fname=trim(trname(n))//'_EPFC'
          call openunit(trim(fname),iu_data,.true.,.true.)
          call readt_parallel(grid,iu_data,nameunit(iu_data),
     &         epfc(:,:,n),1)
          call closeunit(iu_data)
        endif
      enddo
#endif

      return
      end subroutine init_flammability


      subroutine io_flammability(kunit,iaction,ioerr)
!@sum  io_flammabilty reads and writes flammability variables to file
!@auth Greg Faluvegi (based on Jean Lerner io_tracer)
!@ver  1.0 
      use model_com, only: im,jm,ioread,iowrite,irsfic,irsficno,irerun
      use domain_decomp_1d, only: get,grid,am_i_root,
     &     pack_data,unpack_data
      use flammability_com, only: iHfl,iDfl,i0fl,first_prec,PRSfl,
     & DRAfl,HRAfl,maxHR_prec,nday_prec,ravg_prec,flammability

      implicit none

      integer :: kunit   !@var kunit unit number of read/write
      integer :: iaction !@var iaction flag for reading or writing to file
!@var ioerr 1 (or -1) if there is (or is not) an error in i/o
      integer, intent(inout) :: ioerr
!@var header character string label for individual records

      real*8, dimension(:,:), allocatable:: general_glob
      real*8, dimension(:,:,:), allocatable :: DRAfl_glob,HRAfl_glob
      integer :: itm
      character*80 :: header

      INTEGER :: J_0, J_1, J_1H, J_0H

      CALL GET(grid, J_STRT=J_0,     J_STOP=J_1,
     &         J_STRT_HALO=J_0H,J_STOP_HALO=J_1H)

      if(am_i_root()) allocate(
     &     general_glob(IM,JM),
     &     DRAfl_glob(IM,JM,nday_prec),
     &     HRAfl_glob(IM,JM,maxHR_prec)
     &     )

      SELECT CASE (IACTION)

      CASE (:IOWRITE) ! output to end-of-month restart file

       header='CALCULATE_FLAMMABILITY: DRAfl(i,j,days)'
        do itm=1,nday_prec
         call pack_data(grid,drafl(:,:,itm),drafl_glob(:,:,itm))
        end do
        if(am_i_root())write(kunit,err=10)header,drafl_glob
       header='CALCULATE_FLAMMABILITY: HRAfl(i,j,hours)'
        do itm=1,maxHR_prec
         call pack_data(grid,hrafl(:,:,itm),hrafl_glob(:,:,itm))
        end do
        if(am_i_root())write(kunit,err=10)header,hrafl_glob
       header='CALCULATE_FLAMMABILITY: PRSfl(i,j)'
        call pack_data(grid,prsfl(:,:),general_glob(:,:))
        if(am_i_root())write(kunit,err=10)header,general_glob
       header='CALCULATE_FLAMMABILITY: i0fl(i,j) (real)'
        call pack_data(grid,i0fl(:,:),general_glob(:,:))
        if(am_i_root())write(kunit,err=10)header,general_glob
       header='CALCULATE_FLAMMABILITY: iDfl(i,j) (real)'
        call pack_data(grid,iDfl(:,:),general_glob(:,:))
        if(am_i_root())write(kunit,err=10)header,general_glob
       header='CALCULATE_FLAMMABILITY: iHfl(i,j) (real)'
        call pack_data(grid,iHfl(:,:),general_glob(:,:))
        if(am_i_root())write(kunit,err=10)header,general_glob
       header='CALCULATE_FLAMMABILITY: first_prec(i,j) (real)'
        call pack_data(grid,first_prec(:,:),general_glob(:,:))
        if(am_i_root())write(kunit,err=10)header,general_glob
       header='CALCULATE_FLAMMABILITY: ravg_prec(i,j)'
        call pack_data(grid,ravg_prec(:,:),general_glob(:,:))
        if(am_i_root())write(kunit,err=10)header,general_glob
       header='CALCULATE_FLAMMABILITY: flammability(i,j)'
        call pack_data(grid,flammability(:,:),general_glob(:,:))
        if(am_i_root())write(kunit,err=10)header,general_glob

      CASE (IOREAD:)          ! input from restart file
        SELECT CASE (IACTION)
        CASE (ioread,irerun,irsfic,irsficno) ! restarts

          if(am_i_root())read(kunit,err=10)header,drafl_glob
          do itm=1,nday_prec
            call unpack_data(grid,drafl_glob(:,:,itm),drafl(:,:,itm))
          end do           
          if(am_i_root())read(kunit,err=10)header,hrafl_glob
          do itm=1,maxHR_prec
            call unpack_data(grid,hrafl_glob(:,:,itm),hrafl(:,:,itm))
          end do
          if(am_i_root())read(kunit,err=10)header,general_glob
          call unpack_data(grid,general_glob(:,:),prsfl(:,:))
          if(am_i_root())read(kunit,err=10)header,general_glob
          call unpack_data(grid,general_glob(:,:),i0fl(:,:))
          if(am_i_root())read(kunit,err=10)header,general_glob
          call unpack_data(grid,general_glob(:,:),iDfl(:,:))
          if(am_i_root())read(kunit,err=10)header,general_glob
          call unpack_data(grid,general_glob(:,:),iHfl(:,:))
          if(am_i_root())read(kunit,err=10)header,general_glob
          call unpack_data(grid,general_glob(:,:),first_prec(:,:))
          if(am_i_root())read(kunit,err=10)header,general_glob
          call unpack_data(grid,general_glob(:,:),ravg_prec(:,:))
          if(am_i_root())read(kunit,err=10)header,general_glob
          call unpack_data(grid,general_glob(:,:),flammability(:,:))

        END SELECT
      END SELECT

      call freemem
      return

 10   ioerr=1
      call freemem
      call stop_model('error in io_flammability',255) 
      return

      contains
      subroutine freemem
      if(am_i_root()) deallocate(general_glob,DRAfl_glob,HRAfl_glob)
      end subroutine freemem

      end subroutine io_flammability

#ifdef NEW_IO
      subroutine def_rsf_flammability(fid)
!@sum  def_rsf_flammability defines flammability array structure in 
!@+    restart files
!@auth Greg Faluvegi (directly from M. Kelley's def_rsf_lakes)
!@ver  beta
      use flammability_com
      use domain_decomp_atm, only : grid
      use pario, only : defvar
      implicit none
      integer fid   !@var fid file id

      call defvar(grid,fid,drafl,'drafl(dist_im,dist_jm,nday_prec)')
      call defvar(grid,fid,hrafl,'hrafl(dist_im,dist_jm,maxHR_prec)')
      call defvar(grid,fid,prsfl,'prsfl(dist_im,dist_jm)')
      call defvar(grid,fid,i0fl,'i0fl(dist_im,dist_jm)') ! real
      call defvar(grid,fid,iDfl,'iDfl(dist_im,dist_jm)') ! real
      call defvar(grid,fid,iHfl,'iHfl(dist_im,dist_jm)') ! real
      call defvar(grid,fid,first_prec,'first_prec(dist_im,dist_jm)')
      call defvar(grid,fid,ravg_prec,'ravg_prec(dist_im,dist_jm)')
      call defvar(grid,fid,flammability,'flammability(dist_im,dist_jm)')

      return
      end subroutine def_rsf_flammability

      subroutine new_io_flammability(fid,iaction)
!@sum  new_io_flammability read/write lake arrays from/to restart files
!@auth Greg Faluvegi (directly from M. Kelley's new_io_lakes)
!@ver  beta new_ prefix avoids name clash with the default version
      use model_com, only : ioread,iowrite
      use domain_decomp_atm, only : grid
      use pario, only : write_dist_data,read_dist_data
      use flammability_com
      implicit none
      integer fid   !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
      select case (iaction)
      case (iowrite)            ! output to restart file
        call write_dist_data(grid, fid, 'drafl', drafl )
        call write_dist_data(grid, fid, 'hrafl', hrafl )
        call write_dist_data(grid, fid, 'prsfl', prsfl )
        call write_dist_data(grid, fid, 'i0fl', i0fl )
        call write_dist_data(grid, fid, 'iDfl', iDfl )
        call write_dist_data(grid, fid, 'iHfl', iHfl )
        call write_dist_data(grid, fid, 'first_prec', first_prec )
        call write_dist_data(grid, fid, 'ravg_prec', ravg_prec )
        call write_dist_data(grid, fid, 'flammability', flammability )
      case (ioread)            ! input from restart file
        call read_dist_data(grid, fid, 'drafl', drafl )
        call read_dist_data(grid, fid, 'hrafl', hrafl )
        call read_dist_data(grid, fid, 'prsfl', prsfl )
        call read_dist_data(grid, fid, 'i0fl', i0fl )
        call read_dist_data(grid, fid, 'iDfl', iDfl )
        call read_dist_data(grid, fid, 'iHfl', iHfl )
        call read_dist_data(grid, fid, 'first_prec', first_prec )
        call read_dist_data(grid, fid, 'ravg_prec', ravg_prec )
        call read_dist_data(grid, fid, 'flammability', flammability )
      end select
      return
      end subroutine new_io_flammability
#endif /* NEW_IO */

      subroutine flammability_drv
!@sum driver routine for flammability potential of surface
!@+   vegetation calculation.
!@auth Greg Faluvegi based on direction from Olga Pechony
!@ver  1.0 
      use model_com, only: im, jm, dtsrc, ptop, p
      use domain_decomp_atm,only: grid, get
      use flammability_com, only: flammability,veg_density,ravg_prec,
     & ravg_prec,iHfl,iDfl,i0fl,first_prec,HRAfl,DRAfl,PRSfl,missing,
     & raP_acc

      use pblcom, only: tsavg, qsavg
      use fluxes, only: prec
      use constant, only: sday, lhe
      use diag_com, only: ij_flam,aij=>aij_loc

      implicit none

      integer :: J_0S, J_1S, I_0H, I_1H, i, j
      logical :: have_south_pole, have_north_pole     
      real*8 :: qsat ! this is a function in UTILDBL.f
                              
      call get(grid, J_STRT_SKP=J_0S, J_STOP_SKP=J_1S,
     &               HAVE_SOUTH_POLE = have_south_pole,
     &               HAVE_NORTH_POLE = have_north_pole)
      call get(grid,I_STRT_HALO=I_0H,I_STOP_HALO=I_1H)

      if(have_north_pole)flammability(I_0H:I_1H,JM)=missing
      if(have_south_pole)flammability(I_0H:I_1H,1) =missing

      do j=J_0S,J_1S
        do i=I_0H,I_1H

          ! update the precipitation running average, then:
          call prec_running_average(prec(i,j),ravg_prec(i,j), 
     &    iHfl(i,j),iDfl(i,j),i0fl(i,j),first_prec(i,j),HRAfl(i,j,:),
     &    DRAfl(i,j,:),PRSfl(i,j))
    
          ! for sub-daily diag purposes, accumulate the running avg:
          raP_acc(i,j)=raP_acc(i,j)+ravg_prec(i,j)

          ! if the first period has elapsed, calculate the flammability
          if(first_prec(i,j)==0)
     &    call calc_flammability(tsavg(i,j),sday*ravg_prec(i,j)/dtsrc,
     &    min(1.d0,qsavg(i,j)/qsat(tsavg(i,j),lhe,p(i,j)+ptop)),
     &    veg_density(i,j),flammability(i,j)) 

          ! update diagnostic
          aij(i,j,ij_flam)=aij(i,j,ij_flam)+flammability(i,j)

        enddo
      enddo
  
      return
      end subroutine flammability_drv

