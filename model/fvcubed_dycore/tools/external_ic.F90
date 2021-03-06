 module external_ic_mod

   use constants_mod,     only: pi, omega, grav, kappa, rdgas, rvgas, cp_air
   use fv_arrays_mod,     only: fv_atmos_type
   use fv_mapz_mod,       only: mappm
   use mpp_domains_mod,   only: mpp_get_tile_id, domain2d
   use init_hydro_mod,    only: p_var
   use mpp_mod,           only: mpp_error, FATAL, NOTE
   use fms_mod,           only: file_exist, read_data
   use mpp_domains_mod,   only: mpp_update_domains
   use fms_io_mod,        only: get_tile_string, field_size
   use fv_diagnostics_mod,only: prt_maxmin
   use fv_mp_mod,         only: gid, domain, tile, ng,         &
                                is,js,ie,je, isd,jsd,ied,jed, fill_corners, YDir
   use fv_grid_utils_mod, only: ptop_min, fc, f0, ew, es, g_sum, vlon, vlat,  &
                                edge_vect_s,edge_vect_n,edge_vect_w,edge_vect_e
   use fv_grid_tools_mod, only: grid, agrid, cubed_sphere,  &
                                dx,dy, dxa,dya, dxc,dyc, area, rarea
   use tracer_manager_mod, only: get_tracer_names, get_number_tracers, get_tracer_index
   use field_manager_mod,  only: MODEL_ATMOS
   use fv_surf_map_mod,    only: surfdrv
   use fv_io_mod,          only: fv_io_read_tracers 

   implicit none
   private
   public get_external_ic, cubed_a2d
   real, parameter:: zvir = rvgas/rdgas - 1.

contains

   subroutine get_external_ic( Atm, fv_domain )

      type(fv_atmos_type), intent(inout) :: Atm(:)
      type(domain2d),      intent(inout) :: fv_domain
      real:: alpha = 0.
      integer i,j, nq

! * Initialize coriolis param:
 
      do j=jsd,jed+1
         do i=isd,ied+1
            fc(i,j) = 2.*omega*( -1.*cos(grid(i,j,1))*cos(grid(i,j,2))*sin(alpha) + &
                                     sin(grid(i,j,2))*cos(alpha) )
         enddo
      enddo

      do j=jsd,jed
         do i=isd,ied
            f0(i,j) = 2.*omega*( -1.*cos(agrid(i,j,1))*cos(agrid(i,j,2))*sin(alpha) + &
                                     sin(agrid(i,j,2))*cos(alpha) )
         enddo
      enddo

      call mpp_update_domains( f0, domain )
      if ( cubed_sphere ) call fill_corners(f0, Atm(1)%npx, Atm(1)%npy, YDir)
 
! Read in cubed_sphere terrain
 
      call get_cubed_sphere_terrain(Atm, fv_domain)
 
! Read in the specified external dataset and do all the needed transformation
      if ( Atm(1)%ncep_ic ) then
           nq = 1
           call get_ncep_ic( Atm, fv_domain, nq )
#ifndef NO_FV_TRACERS
           call fv_io_read_tracers(fv_domain, Atm)
           if(gid==0) write(6,*) 'All tracers except sphum replaced by FV IC'
#endif
      else
           nq = 4
           call get_fv_ic( Atm, fv_domain, nq )
      endif


      call prt_maxmin('T', Atm(1)%pt, is, ie, js, je, ng, Atm(1)%npz, 1., gid==0)

      call p_var(Atm(1)%npz,  is, ie, js, je, Atm(1)%ak(1),  ptop_min,     &
                 Atm(1)%delp, Atm(1)%delz, Atm(1)%pt, Atm(1)%ps,   &
                 Atm(1)%pe,   Atm(1)%peln, Atm(1)%pk, Atm(1)%pkz,  &
                 kappa, Atm(1)%q, ng, Atm(1)%ncnst, Atm(1)%dry_mass,  &
                 Atm(1)%adjust_dry_mass, .true., .true.,         &
                 Atm(1)%hydrostatic, Atm(1)%k_top, Atm(1)%nwat)

  end subroutine get_external_ic



  subroutine get_cubed_sphere_terrain( Atm, fv_domain )
    type(fv_atmos_type), intent(inout) :: Atm(:)
    type(domain2d),      intent(inout) :: fv_domain
    integer              :: ntileMe
    integer, allocatable :: tile_id(:)
    character(len=64)    :: fname
    integer              ::  n
    real ftop

    ntileMe = size(Atm(:))  ! This will have to be modified for mult tiles per PE
                            ! always one at this point

    allocate( tile_id(ntileMe) )
    tile_id = mpp_get_tile_id( fv_domain )
 
    do n=1,ntileMe

       call get_tile_string(fname, 'INPUT/fv_core.res.tile', tile_id(n), '.nc' )

       if( file_exist(fname) ) then
          call read_data(fname, 'phis', Atm(n)%phis(is:ie,js:je),      &
                         domain=fv_domain, tile_count=n)
       else
          call surfdrv(  Atm(n)%npx, Atm(n)%npy, grid, agrid,   &
                         area, dx, dy, dxc, dyc, Atm(n)%phis, gid==0 ) 
          call mpp_error(NOTE,'terrain datasets generated using USGS data')
       endif

    end do
 
    call mpp_update_domains( Atm(1)%phis, domain )
    ftop = g_sum(Atm(1)%phis(is:ie,js:je), is, ie, js, je, ng, area, 1)
 
    call prt_maxmin('ZS', Atm(1)%phis,  is, ie, js, je, ng, 1, 1./grav, gid==0)
    if(gid==0) write(6,*) 'mean terrain height (m)=', ftop/grav
 
    deallocate( tile_id )

  end subroutine get_cubed_sphere_terrain



  subroutine get_ncep_ic( Atm, fv_domain, nq )
      type(fv_atmos_type), intent(inout) :: Atm(:)
      type(domain2d),      intent(inout) :: fv_domain
      integer, intent(in):: nq

      character(len=128) :: fname
      real, allocatable:: ps0(:,:), gz0(:,:), t0(:,:,:), dp0(:,:,:), q0(:,:,:)
      real, allocatable:: ua(:,:,:), va(:,:,:)
      real, allocatable:: lat(:), lon(:), ak0(:), bk0(:)
      integer :: i, j, k, im, jm, km, npz
      integer tsize(4)
      integer  sphum, liq_wat, ice_wat, cld_amt       ! GFDL AM2 physics
      logical found
      real dak, dbk

      npz = Atm(1)%npz

! Zero out all initial tracer fields:
      Atm(1)%q = 0.

! Read in lat-lon FV core restart file
!     fname = "/archive/sjl/CUBED/NCEP20050825/20050825_00Z_T382.nc"
      fname = Atm(1)%res_latlon_dynamics

      if( file_exist(fname) ) then
          call field_size(fname, 'T', tsize, field_found=found)
          if(gid==0) write(*,*) 'Using NCEP restart:', fname 

          if ( found ) then
               im = tsize(1); jm = tsize(2); km = tsize(3)
               if(gid==0)  write(*,*) 'External IC dimensions:', tsize
          else
               call mpp_error(FATAL,'==> Error from get_external_ic:        &
                              field not found')
          endif

! Define the lat-lon coordinate:
          allocate (  lon(im) )
          allocate (  lat(jm) )
 
          allocate ( ak0(km+1) )
          allocate ( bk0(km+1) )
          allocate ( gz0(im,jm) )
          allocate ( ps0(im,jm) )
          allocate (  ua(im,jm,km) )
          allocate (  va(im,jm,km) )
          allocate (  t0(im,jm,km) )
          allocate ( dp0(im,jm,km) )
          allocate (  q0(im,jm,km) )

          call read_data (fname, 'LAT', lat)
          call read_data (fname, 'LON', lon)

          do i=1,im
             lon(i) = lon(i) * (pi/180.)  ! lon(1) = 0.
          enddo
          do j=1,jm
             lat(j) = lat(j) * (pi/180.)
          enddo

          call read_data (fname, 'hyai', ak0)
          call read_data (fname, 'hybi', bk0)

! Limiter to prevent NAN at top during remapping 
          ak0(1) = min(1.e-5, 0.1*bk0(2)*1.E5)

          call read_data (fname, 'PS', ps0)
          if(gid==0) call pmaxmin( 'PS_ncep', ps0, im,    jm, 0.01)

          call read_data (fname, 'PHIS', gz0)
          if(gid==0) call pmaxmin( 'ZS_ncep', gz0, im,    jm, 1./grav)

          call read_data (fname, 'U',     ua)
          call read_data (fname, 'V',     va)
          if(gid==0) call pmaxmin( 'U_ncep',   ua, im*jm, km, 1.)
          if(gid==0) call pmaxmin( 'V_ncep',   va, im*jm, km, 1.)

          call read_data (fname, 'T',     t0)
          if(gid==0) call pmaxmin( 'T_ncep',   t0, im*jm, km, 1.)

! Compute dp0
          do k=1,km
             dak = ak0(k+1)-ak0(k)
             dbk = bk0(k+1)-bk0(k)
             do j=1,jm
                do i=1,im
                   dp0(i,j,k) = dak + dbk*ps0(i,j)
                enddo
             enddo
          enddo

! Read in tracers: only sphum at this point
!         q0 = 0.
          sphum   = get_tracer_index(MODEL_ATMOS, 'sphum')
          call read_data (fname, 'Q',  q0)
      else
          call mpp_error(FATAL,'==> Error from get_external_ic:        &
                         Expected file '//trim(fname)//' does not exist')
      endif

! Horizontal interpolation to the cubed sphere grid center
! remap vertically with terrain adjustment

! nq is assumed to be 1 here:
      call remap_xyz( im, jm, km, npz, nq, nq, lon, lat, ak0, bk0, ps0,   &
                      gz0, ua, va, t0, q0, dp0, Atm )

      deallocate ( ak0 ) 
      deallocate ( bk0 ) 
      deallocate ( ps0 ) 
      deallocate ( gz0 ) 
      deallocate ( t0 ) 
      deallocate ( q0 ) 
      deallocate ( dp0 ) 
      deallocate ( ua ) 
      deallocate ( va ) 
      deallocate ( lat ) 
      deallocate ( lon ) 

  end subroutine get_ncep_ic



  subroutine get_fv_ic( Atm, fv_domain, nq )
      type(fv_atmos_type), intent(inout) :: Atm(:)
      type(domain2d),      intent(inout) :: fv_domain
      integer, intent(in):: nq

      character(len=128) :: fname
      real, allocatable:: ps0(:,:), gz0(:,:), u0(:,:,:), v0(:,:,:), t0(:,:,:), dp0(:,:,:), q0(:,:,:,:)
      real, allocatable:: ua(:,:,:), va(:,:,:)
      real, allocatable:: lat(:), lon(:), ak0(:), bk0(:)
      integer :: i, j, k, im, jm, km, npz
      integer tsize(4)
      integer sphum, liq_wat, ice_wat, cld_amt       ! GFDL AM2 physics
      logical found

      npz = Atm(1)%npz

! Zero out all initial tracer fields:
      Atm(1)%q = 0.

! Read in lat-lon FV core restart file
      fname = Atm(1)%res_latlon_dynamics

      if( file_exist(fname) ) then
          call field_size(fname, 'T', tsize, field_found=found)
          if(gid==0) write(*,*) 'Using lat-lon FV restart:', fname 

          if ( found ) then
               im = tsize(1); jm = tsize(2); km = tsize(3)
               if(gid==0)  write(*,*) 'External IC dimensions:', tsize
          else
               call mpp_error(FATAL,'==> Error from get_external_ic:        &
                              field not found')
          endif

! Define the lat-lon coordinate:
          allocate (  lon(im) )
          allocate (  lat(jm) )

          do i=1,im
             lon(i) = (0.5 + real(i-1)) * 2.*pi/real(im)
          enddo

          do j=1,jm
             lat(j) = -0.5*pi + real(j-1)*pi/real(jm-1)   ! SP to NP 
          enddo
 
          allocate ( ak0(1:km+1) )
          allocate ( bk0(1:km+1) )
          allocate ( ps0(1:im,1:jm) )
          allocate ( gz0(1:im,1:jm) )
          allocate (  u0(1:im,1:jm,1:km) )
          allocate (  v0(1:im,1:jm,1:km) )
          allocate (  t0(1:im,1:jm,1:km) )
          allocate ( dp0(1:im,1:jm,1:km) )

          call read_data (fname, 'ak', ak0)
          call read_data (fname, 'bk', bk0)
          call read_data (fname, 'Surface_geopotential', gz0)
          call read_data (fname, 'U',     u0)
          call read_data (fname, 'V',     v0)
          call read_data (fname, 'T',     t0)
          call read_data (fname, 'DELP', dp0)

! Share the load
          if(gid==0) call pmaxmin( 'ZS_data', gz0, im,    jm, 1./grav)
          if(gid==1) call pmaxmin( 'U_data',   u0, im*jm, km, 1.)
          if(gid==1) call pmaxmin( 'V_data',   v0, im*jm, km, 1.)
          if(gid==2) call pmaxmin( 'T_data',   t0, im*jm, km, 1.)
          if(gid==3) call pmaxmin( 'DEL-P',   dp0, im*jm, km, 0.01)


      else
          call mpp_error(FATAL,'==> Error from get_external_ic:        &
                         Expected file '//trim(fname)//' does not exist')
      endif

! Read in tracers: only AM2 "physics tracers" at this point
      fname = Atm(1)%res_latlon_tracers

      if( file_exist(fname) ) then
          if(gid==0) write(*,*) 'Using lat-lon tracer restart:', fname 

          allocate ( q0(im,jm,km,Atm(1)%ncnst) )
          q0 = 0.

          sphum   = get_tracer_index(MODEL_ATMOS, 'sphum')
          liq_wat = get_tracer_index(MODEL_ATMOS, 'liq_wat')
          ice_wat = get_tracer_index(MODEL_ATMOS, 'ice_wat')
          cld_amt = get_tracer_index(MODEL_ATMOS, 'cld_amt')

          if(gid==0) write(*,*) 'Tracer indices', sphum, liq_wat, ice_wat, cld_amt 

          call read_data (fname, 'SPHUM',   q0(1:im,1:jm,1:km,sphum))

          if ( nq>=4 ) then
          call read_data (fname, 'LIQ_WAT', q0(1:im,1:jm,1:km,liq_wat))
          call read_data (fname, 'ICE_WAT', q0(1:im,1:jm,1:km,ice_wat))
          call read_data (fname, 'CLD_AMT', q0(1:im,1:jm,1:km,cld_amt))
          endif
      else
          call mpp_error(FATAL,'==> Error from get_external_ic:        &
                         Expected file '//trim(fname)//' does not exist')
      endif

! D to A transform on lat-lon grid:
      allocate (  ua(im,jm,km) )
      allocate (  va(im,jm,km) )

      call d2a3d(u0, v0,  ua,  va, im, jm, km, lon)

      deallocate ( u0 ) 
      deallocate ( v0 ) 

      if(gid==4) call pmaxmin( 'UA', ua, im*jm, km, 1.)
      if(gid==4) call pmaxmin( 'VA', va, im*jm, km, 1.)

      do j=1,jm
         do i=1,im
            ps0(i,j) = ak0(1)
         enddo
      enddo

      do k=1,km
         do j=1,jm
            do i=1,im
               ps0(i,j) = ps0(i,j) + dp0(i,j,k)
            enddo
         enddo
      enddo

  if (gid==0) call pmaxmin( 'PS_data (mb)', ps0, im, jm, 0.01)

! Horizontal interpolation to the cubed sphere grid center
! remap vertically with terrain adjustment

      call remap_xyz( im, jm, km, npz, nq, Atm(1)%ncnst, lon, lat, ak0, bk0,   &
                      ps0,  gz0, ua, va, t0, q0, dp0, Atm )

      deallocate ( ak0 ) 
      deallocate ( bk0 ) 
      deallocate ( ps0 ) 
      deallocate ( gz0 ) 
      deallocate ( t0 ) 
      deallocate ( q0 ) 
      deallocate ( dp0 ) 
      deallocate ( ua ) 
      deallocate ( va ) 
      deallocate ( lat ) 
      deallocate ( lon ) 

  end subroutine get_fv_ic


  subroutine remap_xyz( im, jm, km, npz, nq, ncnst, lon, lat, ak0, bk0, ps0, gz0,   &
                        ua, va, ta, qa, dp, Atm )

  type(fv_atmos_type), intent(inout) :: Atm(:)
  integer, intent(in):: im, jm, km, npz, nq, ncnst
  real,    intent(in):: lon(im), lat(jm), ak0(km+1), bk0(km+1)
  real,    intent(in):: gz0(im,jm), ps0(im,jm)
  real,    intent(in), dimension(im,jm,km):: ua, va, ta,  dp
  real,    intent(in), dimension(im,jm,km,ncnst):: qa
! local:
  real, dimension(isd:ied,jsd:jed,npz):: ut, vt   ! winds 
  real, dimension(is:ie,km):: up, vp, tp
  real, dimension(is:ie,km+1):: pe0, pn0
  real pt0(km), gz(km+1), pk0(km+1)
  real qp(is:ie,km,ncnst)
  real, dimension(is:ie,npz):: qn1
  real, dimension(is:ie,npz+1):: pe1, pn1
  real :: rdlon(im)
  real :: rdlat(jm)
  real:: a1, b1, c1, c2, c3, c4
  real:: gzc, psc, pst
  integer i,j,k, i1, i2, jc, i0, j0, iq
  integer  sphum, liq_wat, ice_wat, cld_amt

  sphum   = get_tracer_index(MODEL_ATMOS, 'sphum')
  liq_wat = get_tracer_index(MODEL_ATMOS, 'liq_wat')
  ice_wat = get_tracer_index(MODEL_ATMOS, 'ice_wat')
  cld_amt = get_tracer_index(MODEL_ATMOS, 'cld_amt')

   if ( sphum/=1 ) then
        call mpp_error(FATAL,'SPHUM must be 1st tracer')
   endif

! call prt_maxmin('C_LON', agrid(isd,jsd,1), is, ie, js, je, ng, 1, 180./pi, gid==0)
! call prt_maxmin('C_LAT', agrid(isd,jsd,2), is, ie, js, je, ng, 1, 180./pi, gid==0)

  pk0(1) = ak0(1)**kappa 

  do i=1,im-1
     rdlon(i) = 1. / (lon(i+1) - lon(i))
  enddo
     rdlon(im) = 1. / (lon(1) + 2.*pi - lon(im))

  do j=1,jm-1
     rdlat(j) = 1. / (lat(j+1) - lat(j))
  enddo

! * Interpolate to cubed sphere cell center
  do 5000 j=js,je

     do i=is,ie
        pe0(i,1) = ak0(1)
        pn0(i,1) = log(ak0(1))
     enddo

     do i=is,ie

       if ( agrid(i,j,1)>lon(im) ) then
            i1 = im;     i2 = 1
            a1 = (agrid(i,j,1)-lon(im)) * rdlon(im)
       elseif ( agrid(i,j,1)<lon(1) ) then
            i1 = im;     i2 = 1
            a1 = (agrid(i,j,1)+2.*pi-lon(im)) * rdlon(im)
       else
            do i0=1,im-1
            if ( agrid(i,j,1)>=lon(i0) .and. agrid(i,j,1)<=lon(i0+1) ) then
               i1 = i0;  i2 = i0+1
               a1 = (agrid(i,j,1)-lon(i1)) * rdlon(i0)
               go to 111
            endif
            enddo
       endif

111    continue

       if ( agrid(i,j,2)<lat(1) ) then
            jc = 1
            b1 = 0.
       elseif ( agrid(i,j,2)>lat(jm) ) then
            jc = jm-1
            b1 = 1.
       else
          do j0=1,jm-1
          if ( agrid(i,j,2)>=lat(j0) .and. agrid(i,j,2)<=lat(j0+1) ) then
               jc = j0
               b1 = (agrid(i,j,2)-lat(jc)) * rdlat(jc)
               go to 222
          endif
          enddo
       endif
222    continue

#ifndef DEBUG_REMAP
       if ( a1<0.0 .or. a1>1.0 .or.  b1<0.0 .or. b1>1.0 ) then
            write(*,*) i,j,a1, b1
       endif
#endif
       c1 = (1.-a1) * (1.-b1)
       c2 =     a1  * (1.-b1)
       c3 =     a1  *     b1
       c4 = (1.-a1) *     b1

! Interpolated surface pressure
       psc = c1*ps0(i1,jc  ) + c2*ps0(i2,jc  ) +    &
             c3*ps0(i2,jc+1) + c4*ps0(i1,jc+1)

! Interpolated surface geopotential
       gzc = c1*gz0(i1,jc  ) + c2*gz0(i2,jc  ) +    &
             c3*gz0(i2,jc+1) + c4*gz0(i1,jc+1)

! 3D fields:
       do iq=1,ncnst
          if ( iq==sphum .or. iq==liq_wat .or. iq==ice_wat .or. iq==cld_amt ) then
          do k=1,km
             qp(i,k,iq) = c1*qa(i1,jc,  k,iq) + c2*qa(i2,jc,  k,iq) +  &
                          c3*qa(i2,jc+1,k,iq) + c4*qa(i1,jc+1,k,iq)
          enddo
          endif
       enddo

       do k=1,km
          up(i,k) = c1*ua(i1,jc,  k) + c2*ua(i2,jc,  k) +  &
                    c3*ua(i2,jc+1,k) + c4*ua(i1,jc+1,k)
          vp(i,k) = c1*va(i1,jc,  k) + c2*va(i2,jc,  k) +  &
                    c3*va(i2,jc+1,k) + c4*va(i1,jc+1,k)
          tp(i,k) = c1*ta(i1,jc,  k) + c2*ta(i2,jc,  k) +  &
                    c3*ta(i2,jc+1,k) + c4*ta(i1,jc+1,k)
! Virtual effect:
          tp(i,k) = tp(i,k)*(1.+zvir*qp(i,k,sphum))
       enddo
! Tracers:

       do k=2,km+1
          pe0(i,k) = ak0(k) + bk0(k)*psc
          pn0(i,k) = log(pe0(i,k))
          pk0(k) = pe0(i,k)**kappa
       enddo

#ifdef USE_DATA_ZS
       Atm(1)%  ps(i,j) = psc
       Atm(1)%phis(i,j) = gzc
#else

! * Adjust interpolated ps to model terrain
       gz(km+1) = gzc 
       do k=km,1,-1
           gz(k) = gz(k+1) + rdgas*tp(i,k)*(pn0(i,k+1)-pn0(i,k)) 
       enddo
! Only lowest layer potential temp is needed
          pt0(km) = tp(i,km)/(pk0(km+1)-pk0(km))*(kappa*(pn0(i,km+1)-pn0(i,km)))
       if( Atm(1)%phis(i,j)>gzc ) then
           do k=km,1,-1
              if( Atm(1)%phis(i,j) <  gz(k)  .and.    &
                  Atm(1)%phis(i,j) >= gz(k+1) ) then
                  pst = pk0(k) + (pk0(k+1)-pk0(k))*(gz(k)-Atm(1)%phis(i,j))/(gz(k)-gz(k+1))
                  go to 123
              endif
           enddo
       else
! Extrapolation into the ground
           pst = pk0(km+1) + (gzc-Atm(1)%phis(i,j))/(cp_air*pt0(km))
       endif

123    Atm(1)%ps(i,j) = pst**(1./kappa)
#endif
     enddo   !i-loop
 

! * Compute delp from ps
     do i=is,ie
        pe1(i,1) = Atm(1)%ak(1)
        pn1(i,1) = log(pe1(i,1))
     enddo
     do k=2,npz+1
       do i=is,ie
          pe1(i,k) = Atm(1)%ak(k) + Atm(1)%bk(k)*Atm(1)%ps(i,j)
          pn1(i,k) = log(pe1(i,k))
       enddo
     enddo

     do k=1,npz
        do i=is,ie
           Atm(1)%delp(i,j,k) = pe1(i,k+1) - pe1(i,k)
        enddo
     enddo
 
! Use kord=4 for winds; kord=8 for others
!------
! map u
!------
      call mappm(km, pe0, up, npz, pe1, qn1, is,ie, -1, 4)
      do k=1,npz
         do i=is,ie
            ut(i,j,k) = qn1(i,k)
         enddo
      enddo
!------
! map v
!------
      call mappm(km, pe0, vp, npz, pe1, qn1, is,ie, -1, 4)
      do k=1,npz
         do i=is,ie
            vt(i,j,k) = qn1(i,k)
         enddo
      enddo

!---------------
! map tracers
!----------------
      do iq=1,ncnst
! Note: AM2 physics tracers only
         if ( iq==sphum .or. iq==liq_wat .or. iq==ice_wat .or. iq==cld_amt ) then
         call mappm(km, pe0, qp(is,1,iq), npz, pe1,  qn1, is,ie, 0, 8)
         do k=1,npz
            do i=is,ie
               Atm(1)%q(i,j,k,iq) = qn1(i,k)
            enddo
         enddo
         endif
      enddo

!-------------------------------------------------------------
! map virtual temperature using geopotential conserving scheme.
!-------------------------------------------------------------
      call mappm(km, pn0, tp, npz, pn1, qn1, is,ie, 1, 8)
      do k=1,npz
         do i=is,ie
            Atm(1)%pt(i,j,k) = qn1(i,k)/(1.+zvir*Atm(1)%q(i,j,k,sphum))
         enddo
      enddo

5000 continue

  call prt_maxmin('PS_model', Atm(1)%ps, is, ie, js, je, ng, 1, 0.01, gid==0)
  call prt_maxmin('UT', ut, is, ie, js, je, ng, npz, 1., gid==0)
  call prt_maxmin('VT', vt, is, ie, js, je, ng, npz, 1., gid==0)

!----------------------------------------------
! winds: lat-lon ON A to Cubed-D transformation:
!----------------------------------------------
  call cubed_a2d(Atm(1)%npx, Atm(1)%npy, npz, ut, vt, Atm(1)%u, Atm(1)%v )

  if (gid==0) write(*,*) 'done remap_xyz'

 end subroutine remap_xyz


 subroutine cubed_a2d( npx, npy, npz, ua, va, u, v )

! Purpose; Transform wind on A grid to D grid

  use mpp_domains_mod,    only: mpp_update_domains

  integer, intent(in):: npx, npy, npz
  real, intent(inout), dimension(isd:ied,jsd:jed,npz):: ua, va
  real, intent(out):: u(isd:ied,  jsd:jed+1,npz)
  real, intent(out):: v(isd:ied+1,jsd:jed  ,npz)
! local:
  real v3(is-1:ie+1,js-1:je+1,3)
  real ue(is-1:ie+1,js:je+1,3)    ! 3D winds at edges
  real ve(is:ie+1,js-1:je+1,  3)    ! 3D winds at edges
  real, dimension(is:ie):: ut1, ut2, ut3
  real, dimension(js:je):: vt1, vt2, vt3
  integer i, j, k, im2, jm2

  call mpp_update_domains(ua, domain, whalo=1, ehalo=1, shalo=1, nhalo=1, complete=.false.)
  call mpp_update_domains(va, domain, whalo=1, ehalo=1, shalo=1, nhalo=1, complete=.true.)

    im2 = (npx-1)/2
    jm2 = (npy-1)/2

    do k=1, npz
! Compute 3D wind on A grid
       do j=js-1,je+1
          do i=is-1,ie+1
             v3(i,j,1) = ua(i,j,k)*vlon(i,j,1) + va(i,j,k)*vlat(i,j,1)
             v3(i,j,2) = ua(i,j,k)*vlon(i,j,2) + va(i,j,k)*vlat(i,j,2)
             v3(i,j,3) = ua(i,j,k)*vlon(i,j,3) + va(i,j,k)*vlat(i,j,3)
          enddo
       enddo

! A --> D
! Interpolate to cell edges
       do j=js,je+1
          do i=is-1,ie+1
             ue(i,j,1) = 0.5*(v3(i,j-1,1) + v3(i,j,1))
             ue(i,j,2) = 0.5*(v3(i,j-1,2) + v3(i,j,2))
             ue(i,j,3) = 0.5*(v3(i,j-1,3) + v3(i,j,3))
          enddo
       enddo

       do j=js-1,je+1
          do i=is,ie+1
             ve(i,j,1) = 0.5*(v3(i-1,j,1) + v3(i,j,1))
             ve(i,j,2) = 0.5*(v3(i-1,j,2) + v3(i,j,2))
             ve(i,j,3) = 0.5*(v3(i-1,j,3) + v3(i,j,3))
          enddo
       enddo

! --- E_W edges (for v-wind):
     if ( is==1 ) then
       i = 1
       do j=js,je
        if ( j>jm2 ) then
             vt1(j) = edge_vect_w(j)*ve(i,j-1,1)+(1.-edge_vect_w(j))*ve(i,j,1)
             vt2(j) = edge_vect_w(j)*ve(i,j-1,2)+(1.-edge_vect_w(j))*ve(i,j,2)
             vt3(j) = edge_vect_w(j)*ve(i,j-1,3)+(1.-edge_vect_w(j))*ve(i,j,3)
        else
             vt1(j) = edge_vect_w(j)*ve(i,j+1,1)+(1.-edge_vect_w(j))*ve(i,j,1)
             vt2(j) = edge_vect_w(j)*ve(i,j+1,2)+(1.-edge_vect_w(j))*ve(i,j,2)
             vt3(j) = edge_vect_w(j)*ve(i,j+1,3)+(1.-edge_vect_w(j))*ve(i,j,3)
        endif
       enddo
       do j=js,je
          ve(i,j,1) = vt1(j)
          ve(i,j,2) = vt2(j)
          ve(i,j,3) = vt3(j)
       enddo
     endif

     if ( (ie+1)==npx ) then
       i = npx
       do j=js,je
        if ( j>jm2 ) then
             vt1(j) = edge_vect_e(j)*ve(i,j-1,1)+(1.-edge_vect_e(j))*ve(i,j,1)
             vt2(j) = edge_vect_e(j)*ve(i,j-1,2)+(1.-edge_vect_e(j))*ve(i,j,2)
             vt3(j) = edge_vect_e(j)*ve(i,j-1,3)+(1.-edge_vect_e(j))*ve(i,j,3)
        else
             vt1(j) = edge_vect_e(j)*ve(i,j+1,1)+(1.-edge_vect_e(j))*ve(i,j,1)
             vt2(j) = edge_vect_e(j)*ve(i,j+1,2)+(1.-edge_vect_e(j))*ve(i,j,2)
             vt3(j) = edge_vect_e(j)*ve(i,j+1,3)+(1.-edge_vect_e(j))*ve(i,j,3)
        endif
       enddo
       do j=js,je
          ve(i,j,1) = vt1(j)
          ve(i,j,2) = vt2(j)
          ve(i,j,3) = vt3(j)
       enddo
     endif

! N-S edges (for u-wind):
     if ( js==1 ) then
       j = 1
       do i=is,ie
        if ( i>im2 ) then
             ut1(i) = edge_vect_s(i)*ue(i-1,j,1)+(1.-edge_vect_s(i))*ue(i,j,1)
             ut2(i) = edge_vect_s(i)*ue(i-1,j,2)+(1.-edge_vect_s(i))*ue(i,j,2)
             ut3(i) = edge_vect_s(i)*ue(i-1,j,3)+(1.-edge_vect_s(i))*ue(i,j,3)
        else
             ut1(i) = edge_vect_s(i)*ue(i+1,j,1)+(1.-edge_vect_s(i))*ue(i,j,1)
             ut2(i) = edge_vect_s(i)*ue(i+1,j,2)+(1.-edge_vect_s(i))*ue(i,j,2)
             ut3(i) = edge_vect_s(i)*ue(i+1,j,3)+(1.-edge_vect_s(i))*ue(i,j,3)
        endif
       enddo
       do i=is,ie
          ue(i,j,1) = ut1(i)
          ue(i,j,2) = ut2(i)
          ue(i,j,3) = ut3(i)
       enddo
     endif

     if ( (je+1)==npy ) then
       j = npy
       do i=is,ie
        if ( i>im2 ) then
             ut1(i) = edge_vect_n(i)*ue(i-1,j,1)+(1.-edge_vect_n(i))*ue(i,j,1)
             ut2(i) = edge_vect_n(i)*ue(i-1,j,2)+(1.-edge_vect_n(i))*ue(i,j,2)
             ut3(i) = edge_vect_n(i)*ue(i-1,j,3)+(1.-edge_vect_n(i))*ue(i,j,3)
        else
             ut1(i) = edge_vect_n(i)*ue(i+1,j,1)+(1.-edge_vect_n(i))*ue(i,j,1)
             ut2(i) = edge_vect_n(i)*ue(i+1,j,2)+(1.-edge_vect_n(i))*ue(i,j,2)
             ut3(i) = edge_vect_n(i)*ue(i+1,j,3)+(1.-edge_vect_n(i))*ue(i,j,3)
        endif
       enddo
       do i=is,ie
          ue(i,j,1) = ut1(i)
          ue(i,j,2) = ut2(i)
          ue(i,j,3) = ut3(i)
       enddo
     endif

     do j=js,je+1
        do i=is,ie
           u(i,j,k) =  ue(i,j,1)*es(1,i,j,1) +  &
                       ue(i,j,2)*es(2,i,j,1) +  &
                       ue(i,j,3)*es(3,i,j,1)
        enddo
     enddo
     do j=js,je
        do i=is,ie+1
           v(i,j,k) = ve(i,j,1)*ew(1,i,j,2) +  &
                      ve(i,j,2)*ew(2,i,j,2) +  &
                      ve(i,j,3)*ew(3,i,j,2)
        enddo
     enddo
 
   enddo         ! k-loop

 end subroutine cubed_a2d



 subroutine d2a3d(u, v,  ua,   va,  im,  jm, km, lon)
      integer, intent(in):: im, jm, km           ! Dimensions
      real, intent(in ) :: lon(im)
      real, intent(in ), dimension(im,jm,km):: u, v
      real, intent(out), dimension(im,jm,km):: ua, va
! local
      real :: coslon(im),sinlon(im)    ! Sine and cosine in longitude
      integer i, j, k
      integer imh
      real un, vn, us, vs

      integer :: ks, ke

      imh = im/2

      do i=1,im
         sinlon(i) = sin(lon(i))
         coslon(i) = cos(lon(i))
      enddo

      do k=1,km
         do j=2,jm-1
            do i=1,im
               ua(i,j,k) = 0.5*(u(i,j,k) + u(i,j+1,k))
            enddo
         enddo

         do j=2,jm-1
            do i=1,im-1
               va(i,j,k) = 0.5*(v(i,j,k) + v(i+1,j,k))
            enddo
            va(im,j,k) = 0.5*(v(im,j,k) + v(1,j,k))
         enddo

! Projection at SP
             us = 0.
             vs = 0.
             do i=1,imh
                us = us + (ua(i+imh,2,k)-ua(i,2,k))*sinlon(i)      &
                     + (va(i,2,k)-va(i+imh,2,k))*coslon(i)
                vs = vs + (ua(i+imh,2,k)-ua(i,2,k))*coslon(i)      &
                     + (va(i+imh,2,k)-va(i,2,k))*sinlon(i)
             enddo
             us = us/im
             vs = vs/im
             do i=1,imh
                ua(i,1,k)   = -us*sinlon(i) - vs*coslon(i)
                va(i,1,k)   =  us*coslon(i) - vs*sinlon(i)
                ua(i+imh,1,k)   = -ua(i,1,k)
                va(i+imh,1,k)   = -va(i,1,k)
             enddo

! Projection at NP
             un = 0.
             vn = 0.
             do i=1,imh
                un = un + (ua(i+imh,jm-1,k)-ua(i,jm-1,k))*sinlon(i)    &
                     + (va(i+imh,jm-1,k)-va(i,jm-1,k))*coslon(i)
                vn = vn + (ua(i,jm-1,k)-ua(i+imh,jm-1,k))*coslon(i)    &
                     + (va(i+imh,jm-1,k)-va(i,jm-1,k))*sinlon(i)
             enddo

             un = un/im
             vn = vn/im
             do i=1,imh
                ua(i,jm,k) = -un*sinlon(i) + vn*coslon(i)
                va(i,jm,k) = -un*coslon(i) - vn*sinlon(i)
                ua(i+imh,jm,k) = -ua(i,jm,k)
                va(i+imh,jm,k) = -va(i,jm,k)
             enddo
      enddo

  end subroutine d2a3d



  subroutine pmaxmin( qname, a, im, jm, fac )

      character*(*)  qname
      integer im, jm
      integer i, j
      real a(im,jm)

      real qmin(jm), qmax(jm)
      real pmax, pmin
      real fac                     ! multiplication factor

      do j=1,jm
         pmax = a(1,j)
         pmin = a(1,j)
         do i=2,im
            pmax = max(pmax, a(i,j))
            pmin = min(pmin, a(i,j))
         enddo
         qmax(j) = pmax
         qmin(j) = pmin
      enddo
!
! Now find max/min of amax/amin
!
            pmax = qmax(1)
            pmin = qmin(1)
         do j=2,jm
            pmax = max(pmax, qmax(j))
            pmin = min(pmin, qmin(j))
         enddo

      write(*,*) qname, ' max = ', pmax*fac, ' min = ', pmin*fac

 end subroutine pmaxmin



 end module external_ic_mod

