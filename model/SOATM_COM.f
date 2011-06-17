#include "rundeck_opts.h"

      module domain_decomp_atm
      use dist_grid_mod, only : dist_grid,am_i_root,broadcast
#ifdef STANDALONE_HYCOM
      use  hycom_dim, only : grid=>ogrid
#else
      use oceanr_dim, only : grid=>ogrid
#endif
      implicit none

      end module domain_decomp_atm

      module geom
      implicit none

      real*8, dimension(:,:), allocatable ::
     &     lon2d,lat2d,lon2d_dg,lat2d_dg,sinlat2d,coslat2d,axyp

#ifdef STANDALONE_HYCOM
      integer, parameter :: im_nodup=318 ! fix hard-coding
      integer, parameter :: iseam=69
#endif

      end module geom

      subroutine alloc_geom
      use domain_decomp_atm, only : grid
      use geom
      implicit none
      integer :: i_0h, i_1h, j_1h, j_0h
      integer :: i, j, i_0, i_1, j_1, j_0
      integer :: ier

      i_0h = grid%i_strt_halo
      i_1h = grid%i_stop_halo
      j_0h = grid%j_strt_halo
      j_1h = grid%j_stop_halo

      allocate(lon2d(i_0h:i_1h,j_0h:j_1h),
     &         lat2d(i_0h:i_1h,j_0h:j_1h),
     &      sinlat2d(i_0h:i_1h,j_0h:j_1h),
     &      coslat2d(i_0h:i_1h,j_0h:j_1h),
     &      axyp(i_0h:i_1h,j_0h:j_1h),
     &     stat = ier)

      end subroutine alloc_geom

      module diag_com
      use mdiag_com, only : sname_strlen,units_strlen,lname_strlen
      use cdl_mod, only : cdl_type
#ifdef STANDALONE_HYCOM
      use geom, only : im_nodup
#endif
      use domain_decomp_1d, only : dist_grid
      implicit none
      real*8, dimension(:,:,:), allocatable, target :: aij
      integer :: kaij=0
      real*8, allocatable :: scale_aij(:)
      integer, allocatable :: denom_aij(:),ia_aij(:)
      character(len=sname_strlen), allocatable :: sname_aij(:)
      character(len=units_strlen), allocatable :: units_aij(:)
      character(len=lname_strlen), allocatable :: lname_aij(:)

      integer :: ij_foc,ij_roc,ij_rsi,ij_ocnsw,ij_sisw,
     &     ij_ocnalb,ij_sialb,ij_msi,ij_ocheat,ij_sst,ij_sss

      type(cdl_type) :: cdl_ij_template,cdl_aij

#ifdef STANDALONE_HYCOM
      real*8, dimension(:,:,:), allocatable, target :: aij_nodup
      type(dist_grid) :: grid_nodup
#else
      type(dist_grid), pointer :: grid_nodup
#endif

      real*8, dimension(:,:,:), pointer :: aij_ioptr

      end module diag_com

      subroutine alloc_adiag
      use constant, only : rhoi
      use domain_decomp_atm, only : grid
      use domain_decomp_1d, only : init_grid
      use mdiag_com
      use diag_com
      use cdl_mod, only : init_cdl_type,add_dim,add_var,add_varline
      use model_com, only : DTsrc
      use geom, only : lon2d_dg,lat2d_dg
      implicit none
      integer :: i_0h, i_1h, j_1h, j_0h
      integer :: i, j, i_0, i_1, j_1, j_0
      integer :: k,ier
      integer, parameter :: kaij_max=20
      real*8 :: scale_aij_(kaij_max)
      integer :: denom_aij_(kaij_max),ia_aij_(kaij_max)
      character(len=sname_strlen) :: sname_aij_(kaij_max)
      character(len=units_strlen) :: units_aij_(kaij_max)
      character(len=lname_strlen) :: lname_aij_(kaij_max)
      character(len=20) :: ijstr

      i_0h = grid%i_strt_halo
      i_1h = grid%i_stop_halo
      j_0h = grid%j_strt_halo
      j_1h = grid%j_stop_halo


      denom_aij_ = 0
      ia_aij_ = ia_cpl
      scale_aij_ = 1.

      k=0
c
      k=k+1
      IJ_FOC=k
      lname_aij_(k)="Ocean fraction"
      sname_aij_(k)="ocnfr"
      units_aij_(k)="%"
      scale_aij_(k)=100.
c
      k=k+1
      IJ_ROC=k
      lname_aij_(k)="Open ocean fraction"
      sname_aij_(k)="opocnfr"
      units_aij_(k)="%"
      scale_aij_(k)=100.
c
      k=k+1
      IJ_RSI=k
      lname_aij_(k)="Sea ice fraction"
      sname_aij_(k)="icefr"
      units_aij_(k)="%"
      scale_aij_(k)=100.
      denom_aij_(k)=ij_foc
c
      k=k+1
      IJ_OCNSW=k
      lname_aij_(k)="Open ocean incident solar radiation"
      sname_aij_(k)="ocn_incsw"
      units_aij_(k)="W/m^2"
      scale_aij_(k)=1.
      denom_aij_(k)=ij_roc
c
      k=k+1
      IJ_SISW=k
      lname_aij_(k)="Sea ice incident solar radiation"
      sname_aij_(k)="si_incsw"
      units_aij_(k)="W/m^2"
      scale_aij_(k)=1.
      denom_aij_(k)=ij_rsi
c
      k=k+1
      IJ_OCNALB=k
      lname_aij_(k)="Open ocean albedo"
      sname_aij_(k)="ocnalb"
      units_aij_(k)="%"
      scale_aij_(k)=100.
      denom_aij_(k)=IJ_OCNSW
c
      k=k+1
      IJ_SIALB=k
      lname_aij_(k)="Sea ice albedo"
      sname_aij_(k)="sialb"
      units_aij_(k)="%"
      scale_aij_(k)=100.
      denom_aij_(k)=IJ_SISW
c
      k=k+1
      IJ_MSI=k
      lname_aij_(k)="Sea ice thickness"
      sname_aij_(k)="zsi"
      units_aij_(k)="m"
      scale_aij_(k)=1./rhoi
      denom_aij_(k) = IJ_RSI
c
      k=k+1
      IJ_OCHEAT=k
      lname_aij_(k)="Open ocean total surface heat flux"
      sname_aij_(k)="totht"
      units_aij_(k)="W/m^2"
      scale_aij_(k)=1./DTsrc
      denom_aij_(k) = ij_roc
c
      k=k+1
      IJ_SST=k
      lname_aij_(k)="Sea surface temperature"
      sname_aij_(k)="sst"
      units_aij_(k)="C"
      scale_aij_(k)=1.
      denom_aij_(k)=ij_foc
c
      k=k+1
      IJ_SSS=k
      lname_aij_(k)="Sea surface salinity"
      sname_aij_(k)="sss"
      units_aij_(k)="psu"
      scale_aij_(k)=1.
      denom_aij_(k)=ij_foc
c
      kaij = k
      allocate(aij(i_0h:i_1h,j_0h:j_1h,kaij))
      allocate(scale_aij(kaij),denom_aij(kaij),ia_aij(kaij))
      allocate(sname_aij(kaij),units_aij(kaij),lname_aij(kaij))
      scale_aij = scale_aij_(1:kaij)
      denom_aij = denom_aij_(1:kaij)
      ia_aij    =    ia_aij_(1:kaij)
      sname_aij = sname_aij_(1:kaij)
      units_aij = units_aij_(1:kaij)
      lname_aij = lname_aij_(1:kaij)
#ifdef STANDALONE_HYCOM
      call init_grid( grid_nodup, im_nodup, grid%jm_world, 1,
     &     bc_periodic=.true. )
      allocate(aij_nodup(im_nodup,j_0h:j_1h,kaij))
#else
      grid_nodup => grid
#endif

      ijstr='(y,x) ;'
      call init_cdl_type('cdl_aij',cdl_aij)
      call add_dim(cdl_aij,'x',grid_nodup%im_world)
      call add_dim(cdl_aij,'y',grid%jm_world)
      call add_var(cdl_aij,'float lon'//trim(ijstr),
     &     units='degrees_east')
      call add_var(cdl_aij,'float lat'//trim(ijstr),
     &     units='degrees_north')

      do k=1,kaij
        if(trim(sname_aij(k)).eq.'unused') cycle
        call add_var(cdl_aij,
     &       'float '//trim(sname_aij(k))//trim(ijstr),
     &       long_name=trim(lname_aij(k)),
     &       units=trim(units_aij(k)) )
        if(denom_aij(k) .ne. 0) then
          call add_varline(cdl_aij,trim(sname_aij(k))//
     &         ':missing_value = -1.e30f ;')
        endif
        call add_varline(cdl_aij,trim(sname_aij(k))//
     &       ':coordinates = "lat lon" ;')
      enddo

      allocate(lon2d_dg(grid_nodup%im_world,j_0h:j_1h),
     &         lat2d_dg(grid_nodup%im_world,j_0h:j_1h))


      end subroutine alloc_adiag

      module fluxes
      use domain_decomp_atm, only : grid
      use exchange_types, only : atmocn_xchng_vars,atmice_xchng_vars
      implicit none

!@var atmocn,atmice derived-type strucures containing
!@+   variables (or pointers thereto) needed for atmospheric
!@+   interactions with water and floating ice
      type(atmocn_xchng_vars) :: atmocn
      type(atmice_xchng_vars) :: atmice

      end module fluxes

      subroutine alloc_fluxes
      use constant, only : radian
      use geom
      use exchange_types, only : alloc_xchng_vars
      use domain_decomp_atm, only : grid
      use fluxes
#ifdef STANDALONE_HYCOM
      use hycom_arrays_glob, only : scatter_hycom_arrays
      use hycom_arrays, only : depths,lonij,latij,scp2
#else
      use ocean, only : focean,dxypo,olon_dg,olat_dg,dlatm,sinic,cosic
#endif
      implicit none
      integer :: i_0h, i_1h, j_1h, j_0h
      integer :: i, j, i_0, i_1, j_1, j_0
      integer :: ier

      i_0h = grid%i_strt_halo
      i_1h = grid%i_stop_halo
      j_0h = grid%j_strt_halo
      j_1h = grid%j_stop_halo

      i_0 = grid%i_strt
      i_1 = grid%i_stop
      j_0 = grid%j_strt
      j_1 = grid%j_stop

      call alloc_xchng_vars(grid,atmocn)
      atmocn%grid => grid
      call alloc_xchng_vars(grid,atmice)
      atmice%grid => grid

#ifdef STANDALONE_HYCOM
      call scatter_hycom_arrays
      atmocn%focean(:,:) = 0.
      do j=grid%j_strt,grid%j_stop
      do i=grid%i_strt,grid%i_stop
        if(depths(i,j).ne.0.) atmocn%focean(i,j) = 1.
        axyp(i,j) = scp2(i,j)
        lon2d(i,j) = lonij(i,j,3)-180.
        if(lon2d(i,j).lt.0.) lon2d(i,j)=lon2d(i,j)+360.
        lon2d(i,j) = lon2d(i,j)*radian
        lat2d(i,j) = latij(i,j,3)*radian
      enddo
      enddo
      call panam_nodup(lonij(:,:,3),lon2d_dg)
      where(lon2d_dg.gt.180.) lon2d_dg=lon2d_dg-360.
      call panam_nodup(latij(:,:,3),lat2d_dg)
#else
      atmocn % dlatm = dlatm
      atmocn % sini(:) = sinic(:)
      atmocn % cosi(:) = cosic(:)
      do j=grid%j_strt,grid%j_stop
        atmocn%focean(:,j) = focean(:,j)
        axyp(:,j) = dxypo(j)
        lon2d_dg(:,j) = olon_dg(:,1)
        lat2d_dg(:,j) = olat_dg(j,1)
        lon2d(:,j) = (180.+olon_dg(:,1))*radian
        lat2d(:,j) = olat_dg(j,1)*radian
      enddo
#endif

      do j=grid%j_strt,grid%j_stop
      do i=grid%i_strt,grid%i_stop
        atmocn%lat(i,j) = lat2d(i,j)
        sinlat2d(i,j) = sin(lat2d(i,j))
        coslat2d(i,j) = cos(lat2d(i,j))
      enddo
      enddo

      atmice%focean = atmocn%focean
      atmice%lat = atmocn%lat
      atmice%dlatm = atmocn%dlatm     ! needed?
      atmice%sini(:) = atmocn%sini(:) ! needed?
      atmice%cosi(:) = atmocn%cosi(:) ! needed?
      deallocate(atmice%prec);  atmice%prec => atmocn%prec
      deallocate(atmice%eprec); atmice%eprec => atmocn%eprec
      deallocate(atmice%srfp);  atmice%srfp => atmocn%srfp

      atmocn%modd5s = -999 ! ocean should not compute conserv diags

      atmocn%gmelt = 0.; atmocn%egmelt = 0. ! no icebergs

      return
      end subroutine alloc_fluxes

#ifdef STANDALONE_HYCOM
      subroutine panam_nodup(xin,xout)
      use domain_decomp_atm, only : grid
      use domain_decomp_1d, only : am_i_root,pack_data,broadcast
      use geom, only : im_nodup,iseam
      implicit none
      real*8 ::
     &      xin(grid%im_world,grid%j_strt_halo:grid%j_stop_halo),
     &     xout(     im_nodup,grid%j_strt_halo:grid%j_stop_halo)
      real*8, allocatable :: xtmp_loc(:,:),xtmp(:,:)
      integer :: i,j,jj
      integer :: jm
      jm = grid%jm_world
      allocate(xtmp_loc(iseam,grid%j_strt_halo:grid%j_stop_halo))
      xtmp_loc(1:iseam,:) = xin(1:iseam,:)
      xout(:,:) = xin(iseam+1:,:)
      allocate(xtmp(iseam,jm))
      call pack_data(grid,xtmp_loc,xtmp)
      call broadcast(grid,xtmp)
      ! fix hard-coding
      do j=max(grid%j_strt,71),min(grid%j_stop,250)
        jj = 141-j; if(jj.lt.1) jj=jj+jm
        xout(1:iseam,j) = xtmp(iseam:1:-1,jj)
      enddo
      deallocate(xtmp_loc,xtmp)
      return
      end subroutine panam_nodup
#endif
