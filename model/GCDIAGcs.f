#include "rundeck_opts.h"
!@sum  This file contains routines which calculate atmospheric
!@+    transport and energetics diagnostics using modelE fields
!@+    on a domain-decomposed cubed sphere grid.
!@+    The routines are built upon functionality provided by the
!@+    cs2ll_utils module that facilitates non-local computations
!@+    such as eddy statistics and spectral analyses.
!@auth M. Kelley

      module gcdiag
      use diag_zonal, only : imlon,imlonh,jmlat
      use cs2ll_utils
      implicit none
      save

!@var cs2ll an object containing info needed for cube-to-latlon
!@+   interpolation and coordination of processors
      type(cs2ll_type) :: cs2ll

!@var gclon,gclat,dxyp lons,lats,areas of the latlon grid used to
!@+   compute AGC diagnostics
      real*8, dimension(imlon) :: gclon
      real*8, dimension(jmlat) :: gclat,dxyp,dxp

!@var gc_xxx indices for AGC accumulations
      integer ::
     &     gc_nptsavg,gc_dp,
     &     gc_u,gc_v,gc_temp,gc_hght,gc_q,gc_theta,gc_vvel,
     &     gc_totntmom,gc_eddntmom,
     &     gc_totntsh,gc_eddntsh,
     &     gc_totntlh,gc_eddntlh,
     &     gc_totntgeo,gc_eddntgeo,
     &     gc_totke,gc_eddke,
     &     gc_totntke,
     &     gc_totvtlh,gc_eddvtlh,
     &     gc_totvtdse,gc_eddvtdse,
     &     gc_totvtmom,gc_eddvtmom,
     &     gc_eddvtgeo,
     &     gc_psi

!@var jlat_eq,jlat_45n,jlat_50n special latitudes for spectral analyses
      integer :: jlat_eq,jlat_50n,jlat_45n

      end module gcdiag

      subroutine gc_defs
!@sum  gc_defs defines metadata for AGC diagnostics and
!@+    initializes the cube-to-latlon infrastructure object cs2ll
!@auth M. Kelley
      use CONSTANT, only : rgas,lhe,bygrav,pi,radius
      use MODEL_COM, only : qcheck
      use GCDIAG
      use DIAG_COM
      use domain_decomp_atm, only : grid,am_i_root
      implicit none
      integer :: ilon,jlat,k,kk
      real*8 :: byim,dlon,dlat,sins,sinn
      character(len=10) :: zstr

      byim = 1./imlon

      dlon = 2.*pi/imlon
      dlat = pi/jmlat
      do ilon=1,imlon
        gclon(ilon) = -pi + dlon*(real(ilon,kind=8)-.5)
      enddo
      sins = -1.
      do jlat=1,jmlat
        gclat(jlat) = -pi/2. + dlat*(real(jlat,kind=8)-.5)
        lat_gc(jlat) = gclat(jlat)*180./pi
        sinn = sin(-pi/2. + dlat*real(jlat,kind=8))
        dxyp(jlat) = radius*radius*dlon*(sinn-sins)
        dxp(jlat) = radius*dlon*cos(gclat(jlat))
        sins = sinn
      enddo
      jlat_eq = jmlat/2 ! slightly off the equator
      jlat_45n = nint(real(jmlat,kind=8)*.75d0)
      jlat_50n = nint(real(jmlat,kind=8)*140d0/180d0)

      call init_cs2ll_type(grid,imlon,jmlat,cs2ll,gclon,gclat)

c
      do k=1,kagcx
         write(sname_gc(k),'(a3,i3.3)') 'AGC',k
         lname_gc(k) = 'unused'
         units_gc(k) = 'unused'
         pow_gc(k) = 0
         scale_gc(k) = 1.
         ia_gc(k) = ia_dga
         denom_gc(k) = 0
         jgrid_gc(k) = 1
         lgrid_gc(k) = ctr_cp
      enddo
c
      k=0
c
      k=k+1
      gc_dp = k
      sname_gc(k) = 'dp'
      lname_gc(k) =  'PRESSURE DIFFERENCES (CP)'
      units_gc(k) = 'mb'
      scale_gc(k) = byim
c
      k=k+1
      gc_nptsavg = k
      sname_gc(k) = 'npts_avg'
      lname_gc(k) = 'NUMBER OF GRIDPOINTS INCLUDED IN AVERAGE'
      units_gc(k) = '1'
c
      k=k+1
      gc_temp = k
      sname_gc(k) = 'temp'
      lname_gc(k) = 'TEMPERATURE'
      units_gc(k) = 'K'
      denom_gc(k) = gc_dp
c
      k=k+1
      gc_hght = k
      sname_gc(k) = 'height'
      lname_gc(k) = 'HEIGHT'
      units_gc(k) = 'm'
      scale_gc(k) = BYGRAV
      denom_gc(k) = gc_dp
c
      k=k+1
      gc_q = k
      sname_gc(k) = 'q'
      lname_gc(k) = 'SPECIFIC HUMIDITY'
      units_gc(k) = 'kg/kg'
      denom_gc(k) = gc_dp
c
      k=k+1
      gc_theta = k
      sname_gc(k) = 'pot_temp'
      lname_gc(k) = 'POTENTIAL TEMPERATURE'
      units_gc(k) = 'K'
      scale_gc(k) = p1000k
      denom_gc(k) = gc_dp
c
      k=k+1
      gc_u = k
      sname_gc(k) = 'u'
      lname_gc(k) = 'ZONAL WIND (U COMPONENT)'
      units_gc(k) = 'm/s'
      denom_gc(k) = gc_dp
c
      k=k+1
      gc_v = k
      sname_gc(k) = 'v'
      lname_gc(k) = 'MERIDIONAL WIND (V COMPONENT)'
      units_gc(k) = 'm/s'
      denom_gc(k) = gc_dp
c
      k=k+1
      gc_vvel = k
      sname_gc(k) = 'vvel'
      lname_gc(k) = 'VERTICAL VELOCITY (POSITIVE UPWARD)'
      units_gc(k) = 'mb/s'
      scale_gc(k) = -byim
      lgrid_gc(k) = edg_cp
c
      k=k+1
      gc_eddke = k
      sname_gc(k) = 'edd_ke'
      lname_gc(k) = 'EDDY KINETIC ENERGY'
      units_gc(k) = 'm^2/s^2'
      scale_gc(k) = .5
      denom_gc(k) = gc_dp
c
      k=k+1
      gc_totke = k
      sname_gc(k) = 'tot_ke'
      lname_gc(k) = 'TOTAL KINETIC ENERGY'
      units_gc(k) = 'm^2/s^2'
      scale_gc(k) = .5
      denom_gc(k) = gc_dp
c
      k=k+1
      gc_eddntsh = k
      sname_gc(k) = 'edd_nt_sh'
      lname_gc(k) = 'NORTH. TRANS. SENS. HT. BY EDDIES'
      units_gc(k) = 'K*m/s'
      denom_gc(k) = gc_dp
c
      k=k+1
      gc_totntsh = k
      sname_gc(k) = 'tot_nt_sh'
      lname_gc(k) = 'TOT NORTH. TRANS. SENS. HT'
      units_gc(k) = 'K*m/s'
      denom_gc(k) = gc_dp
c
      k=k+1
      gc_eddntgeo = k
      sname_gc(k) = 'edd_nt_geo'
      lname_gc(k) = 'NORTH. TRANS. GEOPOT. BY EDDIES'
      units_gc(k) = 'm^3/s^3'
      denom_gc(k) = gc_dp
c
      k=k+1
      gc_totntgeo = k
      sname_gc(k) = 'tot_nt_geo'
      lname_gc(k) = 'TOT NORTH. TRANS. GEOPOT.'
      units_gc(k) = 'm^3/s^3'
      denom_gc(k) = gc_dp
c
      k=k+1
      gc_eddntmom = k
      sname_gc(k) = 'edd_nt_mom'
      lname_gc(k) = 'NORTH TRANS ZON. MOM. BY EDDIES'
      units_gc(k) = 'm^2/s^2'
      denom_gc(k) = gc_dp
c
      k=k+1
      gc_totntmom = k
      sname_gc(k) = 'tot_nt_mom'
      lname_gc(k) = 'TOT NORTH. TRANS. ZON. MOM.'
      units_gc(k) = 'm^2/s^2'
      denom_gc(k) = gc_dp
c
      k=k+1
      gc_eddntlh = k
      sname_gc(k) = 'edd_nt_lh'
      lname_gc(k) = 'NORTH TRANS HUMIDITY BY EDDIES'
      units_gc(k) = 'm/s'
      denom_gc(k) = gc_dp
c
      k=k+1
      gc_totntlh = k
      sname_gc(k) = 'tot_nt_lh'
      lname_gc(k) = 'TOTAL NORTHWARD TRANSPORT OF HUMIDITY'
      units_gc(k) = 'm/s'
      denom_gc(k) = gc_dp
c
      k=k+1
      gc_totntke = k
      sname_gc(k) = 'tot_nt_ke'
      lname_gc(k) = 'TOTAL NORTHWARD TRANSPORT OF KINETIC ENERGY'
      units_gc(k) = 'm^3/s^3'
      scale_gc(k) = .5
      denom_gc(k) = gc_dp
c
      k=k+1
      gc_eddvtdse = k
      sname_gc(k) = 'edd_vt_dse'
      lname_gc(k) = 'VERT TRANS DRY STAT. ENER. BY EDDIES'
      units_gc(k) = 'W/m^2'
      scale_gc(k) = -100.*BYGRAV*BYIM
      lgrid_gc(k) = edg_cp
c
      k=k+1
      gc_totvtdse = k
      sname_gc(k) = 'tot_vt_dse'
      lname_gc(k) = 'TOTAL LGE SCALE VERT. TRANS. OF DRY STAT ENRG'
      units_gc(k) = 'W/m^2'
      scale_gc(k) = -100.*BYGRAV*BYIM
      lgrid_gc(k) = edg_cp
c
      k=k+1
      gc_eddvtlh = k
      sname_gc(k) = 'edd_vt_lh'
      lname_gc(k) = 'VERT TRANS LATENT HEAT BY EDDIES'
      units_gc(k) = 'W/m^2'
      scale_gc(k) = -100.*BYGRAV*BYIM*LHE
      lgrid_gc(k) = edg_cp
c
      k=k+1
      gc_totvtlh = k
      sname_gc(k) = 'tot_vt_lh'
      lname_gc(k) = 'TOTAL LGE SCALE VERT. TRANS. OF LATENT HEAT'
      units_gc(k) = 'W/m^2'
      scale_gc(k) = -100.*BYGRAV*BYIM*LHE
      lgrid_gc(k) = edg_cp
c
      k=k+1
      gc_eddvtgeo = k
      sname_gc(k) = 'vt_geopot_eddy'
      lname_gc(k) = 'VERT. TRANS. OF GEOPOT. ENERGY BY EDDIES'
      units_gc(k) = 'W/m^2'
      scale_gc(k) = -100.*BYGRAV*BYIM
      lgrid_gc(k) = edg_cp
c
      k=k+1
      gc_eddvtmom = k
      sname_gc(k) = 'vt_u_eddy'
      lname_gc(k) = 'EDDY VERTICAL ZONAL MOM. FLUX'
      units_gc(k) = 'm^2/s^2'
      scale_gc(k) = -BYGRAV*BYIM
      lgrid_gc(k) = edg_cp
c
      k=k+1
      gc_totvtmom = k
      sname_gc(k) = 'tot_vt_u'
      lname_gc(k) = 'TOTAL VERTICAL ZONAL MOM. FLUX'
      units_gc(k) = 'm^2/s^2'
      scale_gc(k) = -BYGRAV*BYIM
      lgrid_gc(k) = edg_cp
c
      k=k+1
      gc_psi = k  ! derived
      sname_gc(k) = 'psi'
      lname_gc(k) = 'STREAMFUNCTION'
      units_gc(k) = 'kg/s'
      scale_gc(k) = 100.*BYGRAV
      lgrid_gc(k) = edg_cp
c
      if(k.gt.kagc) then
        if(am_i_root()) write (6,*)
     &       'gc_defs: Increase kagc=',kagc,' to at least ',k
        call stop_model( 'kagc too small', 255 )
      endif
      if (AM_I_ROOT()) then
         write (6,*) 'Number of AGC diagnostics defined: kagcmax=',k
         if(qcheck) then
           do kk=1,k
             write (6,'(i4,'':'',a)') kk,trim(lname_gc(kk))
           end do
         endif
      end if

c
c Declare the dimensions and metadata of AGC output fields using
c netcdf CDL notation.  The C convention for dimension ordering
c must be used (reversed wrt Fortran).  Information needed for
c printing ASCII tables of the output is stored here as well.
c
      cdl_gc = ''
      cdl_gc(1:2)(:) = (/
     &     'netcdf xxx { ', 'dimensions:  ' /)
      write(cdl_gc(3),'(a,i3,a)') '   lat_gc = ',jmlat,' ;'
      write(cdl_gc(4),'(a,i3,a)') '   lat_gc_plus3 = ',jmlat+3,' ;'
      write(cdl_gc(5),'(a,i3,a)') '   plm = ',lm,' ;'
      write(cdl_gc(6),'(a,i3,a)') '   ple = ',lm,' ;'
      write(cdl_gc(7),'(a,i3,a)') '   shnhgm = 3 ;'
      cdl_gc(8:14)(:) = (/
     &     'variables:                           ',
     &     'float lat_gc(lat_gc) ;               ',
     &     '   lat_gc:units = "degrees_north" ;  ',
     &     'float plm(plm) ;                     ',
     &     '   plm:units = "mb" ;                ',
     &     'float ple(ple) ;                     ',
     &     '   ple:units = "mb" ;                '
     &     /)
      kk = count(len_trim(cdl_gc).gt.0)
      do k=1,kagc
        if(trim(units_gc(k)).eq.'unused') cycle
        if(lgrid_gc(k).eq.ctr_ml .or. lgrid_gc(k).eq.ctr_cp) then
          zstr='plm'
        else
          zstr='ple'
        endif
        kk = kk + 1
        cdl_gc(kk) = 'float '//trim(sname_gc(k))//'('//
     &       trim(zstr)//',lat_gc) ;'
        kk = kk + 1
        cdl_gc(kk) = '   '//trim(sname_gc(k))//':long_name = "'//
     &       trim(lname_gc(k))//'" ;'
        kk = kk + 1
        cdl_gc(kk) = '   '//trim(sname_gc(k))//':units = "'//
     &       trim(units_gc(k))//'" ;'
        if(pow_gc(k).ne.0) then
          kk = kk + 1
          write(cdl_gc(kk),'(a,i3,a)')
     &         '   '//trim(sname_gc(k))//':prtpow = ',pow_gc(k),' ;'
        endif
        kk = kk + 1
        cdl_gc(kk) = 'float '//trim(sname_gc(k))//'_hemis('//
     &       trim(zstr)//',shnhgm) ;'
        if(denom_gc(k).gt.0) then
          kk = kk + 1
          cdl_gc(kk) = 'float '//trim(sname_gc(k))//
     &         '_vmean(lat_gc_plus3) ;'
        endif
      enddo
      kk = kk + 1
      cdl_gc(kk) = '}'

      return
      end subroutine gc_defs

      subroutine diagb
!@sum DIAGB calculate constant-pressure latitude-circle diagnostics.
!@+   This version is for a cubed sphere grid.
!@auth M. Kelley
      use constant, only : sha,tf,teeny,rgas
      use model_com, only : dtsrc,mdiag,ls1,lm,sige,ptop,t,p,q
      use geom, only : byaxyp
      use gcdiag
      use diag_com, only : agc=>agc_loc,speca,nspher,klayer,ple
      use dynamics, only : phi,sda,pit,ualij,valij
      use diag_loc, only : tx
      use domain_decomp_atm, only : get, grid, am_i_root, sumxpe,
     &     halo_update
      use fftw_com
      implicit none
      real*8, dimension(imlonh+1,nspher) :: ke,ke_part
      real*8, dimension(imlonh+1) :: xu,xv,xke

      integer :: i,j,l,ldn,lup,kq,nq
      integer :: jj,jlat,jcalc,khemi,kspher,kspher_eq,kspher_45n
      integer, dimension(lm) :: klayer_eq,klayer_45n
      integer, parameter :: nqmax=31
      real*8, dimension(lm,nqmax), target :: lonsums
      real*8, dimension(:), pointer ::
     &     dpsum,dpusum,dpvsum,dptxsum,dpqsum,dpthsum,dpphisum,
     &     nsum,nlesum,wlesum,txlesum,ulesum,
     &     philesum,qlesum
      real*8, dimension(cs2ll%isd:cs2ll%ied,lm) ::
     &     dp,dpu,dpv,dptx,dpq,dpth,dpphi,
     &     wle,txle,ule,vle,phile,qle
      real*8, dimension(cs2ll%isd:cs2ll%ied,lm) ::
     &     txll,thll,qll,phill,omgll
      real*8, dimension(lm,cs2ll%isd:cs2ll%ied) ::
     &     uli,vli
      real*8, dimension(cs2ll%isd:cs2ll%ied) :: psll
      real*8, dimension(lm) :: bydpsum
      real*8, dimension(imlon,lm) :: usqrtdp,vsqrtdp

      real*8 :: dpl,wtdn,wtup,bydp,vi,kedp,bysqrtdp,now
      real*8 :: dpuvi,dptvi,dpqvi,dpphivi
     &     ,dpkei,dpkevi,wui,wti,wzi,wqi
      real*8 :: uzm,vzm,tzm,qzm,phizm,wzm,byrhozm
      real*8 :: eddwui,eddwqi,eddwti,eddwzi
      real*8 :: dpeddkei,dpedduvi,dpeddtvi,dpeddqvi,dpeddphivi
      real*8 :: uprime,vprime,wprime,tprime,qprime,phiprime
      real*8, dimension(lm+1) :: pecp,pedge
      real*8, dimension(lm) :: pmid,lwdn,lwup
      integer, parameter :: lmxmax=2*lm
      integer :: lmx
      real*8, dimension(lmxmax) :: dpx
      integer, dimension(lmxmax) :: lmod,lcp

      real*8, dimension(grid%i_strt_halo:grid%i_stop_halo,
     &                  grid%j_strt_halo:grid%j_stop_halo,lm) :: omega

      integer :: i_0,i_1,j_0,j_1

      j_0=grid%j_strt
      j_1=grid%j_stop
      i_0=grid%i_strt
      i_1=grid%i_stop
c
c set up pointers that facilitate longitudinal sums across PEs
c
      kq=0
c layer midpoint quantities
      kq=kq+1; dpsum => lonsums(:,kq)
      kq=kq+1; nsum => lonsums(:,kq)
      kq=kq+1; dpusum => lonsums(:,kq)
      kq=kq+1; dpvsum => lonsums(:,kq)
      kq=kq+1; dptxsum => lonsums(:,kq)
      kq=kq+1; dpqsum => lonsums(:,kq)
      kq=kq+1; dpthsum => lonsums(:,kq)
      kq=kq+1; dpphisum => lonsums(:,kq)
c edge quantities
      kq=kq+1; nlesum => lonsums(:,kq)
      kq=kq+1; wlesum => lonsums(:,kq)
      kq=kq+1; txlesum => lonsums(:,kq)
      kq=kq+1; ulesum => lonsums(:,kq)
      kq=kq+1; philesum => lonsums(:,kq)
      kq=kq+1; qlesum => lonsums(:,kq)
      nq = kq

c
c some arrays to make the kspher logic more robust
c
      klayer_eq(:) = klayer(:) + 2
      klayer_45n(:) = klayer(:) + 3

      ke_part(:,:)=0.

c
c obtain vertical velocity in mb/s units
c
      do l=1,ls1-2
      do j=j_0,j_1
      do i=i_0,i_1
c        omega(i,j,l) = (sda(i,j,l)/dtsrc+pit(i,j)*sige(l+1))*byaxyp(i,j)
        omega(i,j,l) = (sda(i,j,l)/dtsrc)*byaxyp(i,j)
      enddo
      enddo
      enddo
      do l=ls1-1,lm-1
      do j=j_0,j_1
      do i=i_0,i_1
        omega(i,j,l) = sda(i,j,l)*byaxyp(i,j)/dtsrc
      enddo
      enddo
      enddo
      omega(:,:,lm) = 0d0

c
c Halo updates and fills of cube corners in preparation for interpolation.
c Some quantities may already have been haloed - check later.
c
      call halo_update(grid,p)
      call corner_fill_3D(grid,p,1)
      call halo_update(grid,t)
      call corner_fill_3D(grid,t,lm)
      call halo_update(grid,tx)
      call corner_fill_3D(grid,tx,lm)
      call halo_update(grid,q)
      call corner_fill_3D(grid,q,lm)
      call halo_update(grid,phi)
      call corner_fill_3D(grid,phi,lm)
      call halo_update(grid,omega)
      call corner_fill_3D(grid,omega,lm)
      call halo_update(grid,ualij,jdim=3)
      call corner_fill_4D(grid,ualij,lm,1)
      call halo_update(grid,valij,jdim=3)
      call corner_fill_4D(grid,valij,lm,1)

c
c Loop over latitudes
c
      do jj=1,jmlat

      jlat = cs2ll%jlat_sched(jj) ! the schedule facilitates parallelism
      if(cs2ll%ni(jlat).eq.0) cycle ! no work for this PE at this jlat

c
c horizontal interpolation to this latitude circle
c
      call interp_to_jlat_4D(grid,cs2ll,ualij,uli,lm,1,jlat)
      call interp_to_jlat_4D(grid,cs2ll,valij,vli,lm,1,jlat)

      call interp_to_jlat_3D(grid,cs2ll,p,psll,1,jlat)

      call interp_to_jlat_3D(grid,cs2ll,t,thll,lm,jlat)
      call interp_to_jlat_3D(grid,cs2ll,tx,txll,lm,jlat)
      call interp_to_jlat_3D(grid,cs2ll,q,qll,lm,jlat)
      call interp_to_jlat_3D(grid,cs2ll,phi,phill,lm,jlat)
      call interp_to_jlat_3D(grid,cs2ll,omega,omgll,lm,jlat)

      pecp(2:lm+1) = ple(1:lm)
      pecp(1) = 1d30 ! ensure that all column mass is included
      pedge(ls1:lm+1) = pecp(ls1:lm+1)
      do l=1,lm
        dpsum(l) = 0d0
        nsum(l) = 0d0
        dpusum(l) = 0d0
        dpvsum(l) = 0d0
        dptxsum(l) = 0d0
        dpqsum(l) = 0d0
        dpthsum(l) = 0d0
        dpphisum(l) = 0d0
        nlesum(l) = 0d0
        wlesum(l) = 0d0
        txlesum(l) = 0d0
        ulesum(l) = 0d0
        philesum(l) = 0d0
        qlesum(l) = 0d0
      enddo
      do i=cs2ll%is(jlat),cs2ll%ie(jlat)
        if(psll(i).eq.0.) cycle ! valid lons sometimes noncontiguous

c
c conservative vertical regridding to the midpoint of cp layers
c
        pedge(1:ls1-1) = sige(1:ls1-1)*psll(i)+ptop
        call get_dx_intervals(pedge,lm,pecp,lm,dpx,lmod,lcp,lmx,lmxmax)

        do l=1,lm
          dp(i,l) = 0d0
          dpu(i,l) = 0d0
          dpv(i,l) = 0d0
          dptx(i,l) = 0d0
          dpq(i,l) = 0d0
          dpth(i,l) = 0d0
          dpphi(i,l) = 0d0
        enddo
        do l=1,lmx
          dp(i,lcp(l))    = dp(i,lcp(l))    + dpx(l)
          dpu(i,lcp(l))   = dpu(i,lcp(l))   + dpx(l)*uli(lmod(l),i)
          dpv(i,lcp(l))   = dpv(i,lcp(l))   + dpx(l)*vli(lmod(l),i)
          dptx(i,lcp(l))  = dptx(i,lcp(l))  + dpx(l)*txll(i,lmod(l))
          dpq(i,lcp(l))   = dpq(i,lcp(l))   + dpx(l)*qll(i,lmod(l))
          dpth(i,lcp(l))  = dpth(i,lcp(l))  + dpx(l)*thll(i,lmod(l))
          dpphi(i,lcp(l)) = dpphi(i,lcp(l)) + dpx(l)*phill(i,lmod(l))
        enddo

c
c partial longitudinal sums of cp layer midpoint quantities
c
        do l=1,lm
          dpsum(l) = dpsum(l) + dp(i,l)
          if(dp(i,l).gt.0d0) nsum(l) = nsum(l) + 1d0
          dpusum(l) = dpusum(l) + dpu(i,l)
          dpvsum(l) = dpvsum(l) + dpv(i,l)
          dptxsum(l) = dptxsum(l) + dptx(i,l)
          dpqsum(l) = dpqsum(l) + dpq(i,l)
          dpthsum(l) = dpthsum(l) + dpth(i,l)
          dpphisum(l) = dpphisum(l) + dpphi(i,l)
        enddo

c
c vertical interpolation to cp layer edges
c if a layer edge is underground, wtdn==wtup==0
c
        pmid(:) = .5*(pedge(1:lm)+pedge(2:lm+1))
        call xintrp(pmid,lm,pecp(2:lm),lm-1,lwdn,lwup)
        do l=1,lm-1
          ldn = lwdn(l)
          lup = lwup(l)
          wtdn = lwdn(l)-ldn
          wtup = lwup(l)-lup
          ule(i,l) = wtdn*uli(ldn,i)+wtup*uli(lup,i)
          vle(i,l) = wtdn*vli(ldn,i)+wtup*vli(lup,i)
          txle(i,l) = wtdn*txll(i,ldn)+wtup*txll(i,lup)
          qle(i,l) = wtdn*qll(i,ldn)+wtup*qll(i,lup)
          phile(i,l) = wtdn*phill(i,ldn)+wtup*phill(i,lup)
c
c partial longitudinal sums of cp layer edge quantities
c
          if(wtdn+wtup.gt.0.) nlesum(l) = nlesum(l) + 1d0
          txlesum(l) = txlesum(l) + txle(i,l)
          ulesum(l) = ulesum(l) + ule(i,l)
          philesum(l) = philesum(l) + phile(i,l)
          qlesum(l) = qlesum(l) + qle(i,l)
        enddo
        call xintrp(pedge(2:lm),lm-1,pecp(2:lm),lm-1,lwdn,lwup)
        do l=1,lm-1
          ldn = lwdn(l)
          lup = lwup(l)
          wtdn = lwdn(l)-ldn
          wtup = lwup(l)-lup
          wle(i,l) = wtdn*omgll(i,ldn)+wtup*omgll(i,lup)
          wlesum(l) = wlesum(l) + wle(i,l)
        enddo

      enddo ! end loop over the longitudes on this PE

c
c Add the contributions from different PEs to the sums around this
c latitude circle.  After this call, pointers to partial longitudinal
c sums are pointers to full longitudinal sums.
c
      call sumxpe_zonal(cs2ll,jlat,lm*nq,lonsums)

      bydpsum(:) = 1d0/(dpsum(:)+teeny)

c
c store zonal sums of order-1 quantities on the root PE for this jlat
c
      if(cs2ll%am_i_rootj(jlat)) then
        do l=1,lm
          agc(jlat,l,gc_dp)  = agc(jlat,l,gc_dp) + dpsum(l)
          agc(jlat,l,gc_nptsavg)=agc(jlat,l,gc_nptsavg)+nsum(l)
          agc(jlat,l,gc_u)   = agc(jlat,l,gc_u) + dpusum(l)
          agc(jlat,l,gc_v)   = agc(jlat,l,gc_v) + dpvsum(l)
          agc(jlat,l,gc_temp) = agc(jlat,l,gc_temp) + dptxsum(l)
          agc(jlat,l,gc_hght) = agc(jlat,l,gc_hght) + dpphisum(l)
          agc(jlat,l,gc_q) = agc(jlat,l,gc_q) + dpqsum(l)
          agc(jlat,l,gc_theta) = agc(jlat,l,gc_theta) + dpthsum(l)
          agc(jlat,l,gc_vvel) = agc(jlat,l,gc_vvel) + wlesum(l)
        enddo
      endif

c
c partial longitudinal sums of layer-midpoint higher-order quantities
c
      do l=1,lm
        dpuvi = 0d0
        dptvi = 0d0
        dpqvi = 0d0
        dpphivi = 0d0
        dpkevi = 0d0
        dpkei = 0d0
        uzm = dpusum(l)*bydpsum(l)
        vzm = dpvsum(l)*bydpsum(l)
        tzm = dptxsum(l)*bydpsum(l)
        qzm = dpqsum(l)*bydpsum(l)
        phizm = dpphisum(l)*bydpsum(l)
        dpeddkei = 0d0
        dpedduvi = 0d0
        dpeddtvi = 0d0
        dpeddqvi = 0d0
        dpeddphivi = 0d0
c        dpeddkevi = 0d0
        do i=cs2ll%is(jlat),cs2ll%ie(jlat)
          if(psll(i).eq.0.) cycle ! valid lons sometimes noncontiguous
          bydp = 1d0/(dp(i,l)+teeny)
          vi = dpv(i,l)*bydp
          vprime = vi-vzm
          uprime = dpu(i,l)*bydp-uzm
          tprime = dptx(i,l)*bydp-tzm
          qprime = dpq(i,l)*bydp-qzm
          phiprime = dpphi(i,l)*bydp-phizm
          dpuvi = dpuvi + dpu(i,l)*vi
          dptvi = dptvi + dptx(i,l)*vi
          dpqvi = dpqvi + dpq(i,l)*vi
          dpphivi = dpphivi + dpphi(i,l)*vi
          kedp = (dpu(i,l)**2+dpv(i,l)**2)*bydp
          dpkevi = dpkevi + kedp*vi
          dpkei = dpkei + kedp
c          dpsqi = dpsqi + dp(i,l)*dp(i,l)
          dpeddkei = dpeddkei + dp(i,l)*(uprime**2+vprime**2)
          dpedduvi = dpedduvi + dp(i,l)*uprime*vprime
          dpeddtvi = dpeddtvi + dp(i,l)*tprime*vprime
          dpeddqvi = dpeddqvi + dp(i,l)*qprime*vprime
          dpeddphivi = dpeddphivi + dp(i,l)*phiprime*vprime
        enddo
        agc(jlat,l,gc_totke)    = agc(jlat,l,gc_totke)+dpkei
c        agc(jlat,l,gc_dpsqr)    = agc(jlat,l,gc_dpsqr)+dpsqi
        agc(jlat,l,gc_totntmom) = agc(jlat,l,gc_totntmom) + dpuvi
        agc(jlat,l,gc_totntsh)  = agc(jlat,l,gc_totntsh) + dptvi
        agc(jlat,l,gc_totntke)  = agc(jlat,l,gc_totntke)+dpkevi
        agc(jlat,l,gc_totntgeo) = agc(jlat,l,gc_totntgeo)+dpphivi
        agc(jlat,l,gc_totntlh)  = agc(jlat,l,gc_totntlh)+dpqvi

        agc(jlat,l,gc_eddke)    = agc(jlat,l,gc_eddke) + dpeddkei
        agc(jlat,l,gc_eddntmom) = agc(jlat,l,gc_eddntmom) + dpedduvi
        agc(jlat,l,gc_eddntsh)  = agc(jlat,l,gc_eddntsh) + dpeddtvi
        agc(jlat,l,gc_eddntlh)  = agc(jlat,l,gc_eddntlh) + dpeddqvi
        agc(jlat,l,gc_eddntgeo) = agc(jlat,l,gc_eddntgeo) + dpeddphivi
      enddo

c
c partial longitudinal sums of vertical fluxes: total and eddy
c
      do l=1,lm-1
        wui = 0d0
        wti = 0d0
        wzi = 0d0
        wqi = 0d0
        uzm = ulesum(l)/(nlesum(l)+teeny)
        qzm = qlesum(l)/(nlesum(l)+teeny)
        wzm = wlesum(l)/(nlesum(l)+teeny)
        tzm = txlesum(l)/(nlesum(l)+teeny)
        phizm = philesum(l)/(nlesum(l)+teeny)
        byrhozm = rgas*tzm/ple(l)
        eddwui = 0d0
        eddwqi = 0d0
        eddwti = 0d0
        eddwzi = 0d0
        do i=cs2ll%is(jlat),cs2ll%ie(jlat)
          if(psll(i).eq.0.) cycle ! valid lons sometimes noncontiguous
          wui = wui + wle(i,l)*ule(i,l)
          wti = wti + wle(i,l)*txle(i,l)
          wqi = wqi + wle(i,l)*qle(i,l)
          wzi = wzi + wle(i,l)*phile(i,l)
          wprime = wle(i,l)-wzm
          uprime = ule(i,l)-uzm
          qprime = qle(i,l)-qzm
          tprime = txle(i,l)-tzm
          phiprime = phile(i,l)-phizm
          eddwui = eddwui + uprime*wprime
          eddwqi = eddwqi + qprime*wprime
          eddwti = eddwti + tprime*wprime
          eddwzi = eddwzi + phiprime*wprime
        enddo

        agc(jlat,l,gc_totvtdse) = agc(jlat,l,gc_totvtdse)+sha*wti+wzi
        agc(jlat,l,gc_totvtlh) = agc(jlat,l,gc_totvtlh)+wqi
        agc(jlat,l,gc_totvtmom) = agc(jlat,l,gc_totvtmom)+wui*byrhozm

        agc(jlat,l,gc_eddvtdse) = agc(jlat,l,gc_eddvtdse)+
     &       sha*eddwti+eddwzi
        agc(jlat,l,gc_eddvtmom) = agc(jlat,l,gc_eddvtmom)+eddwui*byrhozm
        agc(jlat,l,gc_eddvtlh)  = agc(jlat,l,gc_eddvtlh)+eddwqi
        agc(jlat,l,gc_eddvtgeo) = agc(jlat,l,gc_eddvtgeo)+eddwzi

      enddo

c
c spectral analysis of KE
c
      do l=1,lm
        do i=cs2ll%is(jlat),cs2ll%ie(jlat)
          if(psll(i).eq.0.) then ! valid lons sometimes noncontiguous
            dpu(i,l) = 0.; dpv(i,l) = 0. ! for pack_zonal
          else
            bysqrtdp = 1d0/sqrt(dp(i,l)+teeny)
            dpu(i,l) = dpu(i,l)*bysqrtdp
            dpv(i,l) = dpv(i,l)*bysqrtdp
          endif
        enddo
      enddo
      call pack_zonal(cs2ll,jlat,lm,dpu,usqrtdp)
      call pack_zonal(cs2ll,jlat,lm,dpv,vsqrtdp)
      jcalc = cs2ll%jlat_calc(jlat)
      if(jcalc.gt.0) then
        if(jcalc.gt.jmlat/2) then
          khemi = 2
        else
          khemi = 1
        endif
        do l=1,lm
          kspher = klayer(l) + khemi-1
          kspher_eq = klayer_eq(l)
          kspher_45n = klayer_45n(l)
          call ffte(usqrtdp(:,l),xu)
          call ffte(vsqrtdp(:,l),xv)
          xke(:) = (xu(:)+xv(:))*dxyp(jcalc)
          ke_part(:,kspher) = ke_part(:,kspher) + xke(:)
          if(jcalc.eq.jlat_eq) then
            ke_part(:,kspher_eq) = ke_part(:,kspher_eq) + xke(:)
          elseif(jcalc.eq.jlat_45n) then
            ke_part(:,kspher_45n) = ke_part(:,kspher_45n) + xke(:)
          endif
        enddo
      endif

      enddo ! end latitude loop

c
c sum up spectral analysis contributions
c
      call sumxpe(ke_part,ke)
      if(am_i_root()) speca(:,18,:)=speca(:,18,:)+ke(:,:)

      call timer (now,mdiag)

      return

      contains

      subroutine xintrp(xsrc,nsrc,xdst,ndst,lw1,lw2)
      implicit none
      integer :: nsrc,ndst
      real*8, dimension(nsrc) :: xsrc
      real*8, dimension(ndst) :: xdst
      real*8, dimension(ndst) :: lw1,lw2
      real*8 :: wt1,wt2,dx,s
      integer :: ldst,lsrc,lout
      s = sign(1d0,xsrc(2)-xsrc(1))
      lout = 1
      do while(s*xdst(lout).lt.s*xsrc(1))
        lw1(lout) = 1
        lw2(lout) = 1
        lout = lout+1
      enddo
      lsrc = 1
      do ldst=lout,ndst
        if(s*xdst(ldst).gt.s*xsrc(nsrc)) exit
        do while(s*xsrc(lsrc+1).lt.s*xdst(ldst))
          lsrc = lsrc+1
        enddo
        if(xsrc(lsrc+1).eq.xdst(ldst)) then
          lw1(ldst) = real(lsrc+1,kind=8)+.5
          lw2(ldst) = real(lsrc+1,kind=8)+.5
        else
          dx = xsrc(lsrc+1)-xsrc(lsrc)
          wt2 = (xdst(ldst)-xsrc(lsrc))/dx
          wt1 = 1.-wt2
          lw1(ldst) = real(lsrc  ,kind=8)+wt1
          lw2(ldst) = real(lsrc+1,kind=8)+wt2
        endif
      enddo
      do lout=ldst,ndst
        lw1(lout) = nsrc
        lw2(lout) = nsrc
      enddo
      return
      end subroutine xintrp

      end subroutine diagb

      SUBROUTINE DIAG5A (M5,NDT)
!@sum DIAG5A calculate KE/PE spectral diagnostics around latitude circles.
!@+   This version is for a cubed sphere grid.
!@auth cube-enabled by M. Kelley
C****
C**** THIS DIAGNOSTICS ROUTINE PRODUCES A SPECTRAL ANALYSIS OF KINETIC
C**** AND AVAILABLE POTENTIAL ENERGIES AND THEIR TRANSFER RATES BY
C**** VARIOUS ATMOSPHERIC PROCESSES.
C****
C**** THE PARAMETER M INDICATES WHAT IS STORED IN SPECA(N,M,KSPHER),
C**** IT ALSO INDICATES WHEN DIAG5A IS BEING CALLED.
C**** M=1  MEAN STANDING KINETIC ENERGY            BEFORE SOURCES
C****   2  MEAN KINETIC ENERGY                     BEFORE DYNAMICS
C****   3  MEAN POTENTIAL ENERGY
C****   4  CONVERSION OF K.E. BY ADVECTION         AFTER ADVECTION
C****   5  CONVERSION OF K.E. BY CORIOLIS FORCE    AFTER CORIOLIS TERM
C****   6  CONVERSION FROM P.E. INTO K.E.          AFTER PRESS GRAD FORC
C****   7  CHANGE OF K.E. BY DYNAMICS              AFTER DYNAMICS
C****   8  CHANGE OF P.E. BY DYNAMICS
C****   9  CHANGE OF K.E. BY CONDENSATION          AFTER CONDENSATION
C****  10  CHANGE OF P.E. BY CONDENSATION
C****  11  CHANGE OF P.E. BY RADIATION             AFTER RADIATION
C****  12  CHANGE OF K.E. BY SURFACE               AFTER SURFACE
C****  13  CHANGE OF P.E. BY SURFACE
C****  14  CHANGE OF K.E. BY FILTER                AFTER FILTER
C****  15  CHANGE OF P.E. BY FILTER
C****  16  CHANGE OF K.E. BY DAILY                 AFTER DAILY
C****  17  CHANGE OF P.E. BY DAILY
C****  18  UNUSED
C****  19  LAST KINETIC ENERGY
C****  20  LAST POTENTIAL ENERGY
C****
      use constant, only : sha
      use model_com, only : lm,dsig,idacc,mdiag,
     &     p,ptop,sig,t,zatmo
      use geom, only : areag,axyp,lat2d
      use diag_com, only : speca,atpe,nspher,kspeca,klayer
c      use diag_com, only : ajl=>ajl_loc,jl_ape
      use diag_loc, only : lupa,ldna
      use dynamics, only : ualij,valij,pk,pdsig,sqrtp
      use domain_decomp_atm, only : grid,am_i_root,sumxpe,esmf_bcast,
     &     halo_update
      use gcdiag
      use fftw_com
      implicit none
      integer :: m5,ndt

      real*8, dimension(imlonh+1) :: xu,xv,xke,xape
      real*8, dimension(imlonh+1,nspher) :: ke,ape,ke_part,ape_part
      real*8, dimension(2) :: tpe,tpe_part
      real*8, dimension(lm) :: thgm,gmean,thgm_part,gmean_part
      integer, dimension(kspeca), parameter ::
     &     mtpeof=(/0,0,1,0,0,0,0,2,0,3,  4,0,5,0,6,0,7,0,0,8/)
      integer :: i,j,k,ks,kspher,l,ldn,lup,mape,mke,mtpe,n
      integer :: jj,jlat,jcalc,khemi,kspher_eq,kspher_45n
      real*8 :: sumt, now
      integer, dimension(lm) :: klayer_eq,klayer_45n
      logical :: do_ke
      real*8, dimension(cs2ll%isd:cs2ll%ied,lm) ::
     &     tll,u_tmp,v_tmp,t_tmp
      real*8, dimension(lm,cs2ll%isd:cs2ll%ied) :: uli,vli
      real*8, dimension(cs2ll%isd:cs2ll%ied) :: psll,sqrtpa
      real*8, dimension(imlon,lm) :: u_il,v_il,t_il

      integer :: i_0,i_1,j_0,j_1

      j_0=grid%j_strt
      j_1=grid%j_stop
      i_0=grid%i_strt
      i_1=grid%i_stop

C****
C**** Note: KSPHER has been re-arranged from previous models to better
C****       deal with optionally resolved stratosphere. The higher
C****       values are only used if the model top is high enough.
C****
C**** KSPHER=1 SOUTHERN TROPOSPHERE         2 NORTHERN TROPOSPHERE
C****        3 EQUATORIAL TROPOSPHERE       4 45 DEG NORTH TROPOSPHERE
C****
C****        5 SOUTHERN LOW STRATOSPHERE    6 NORTHERN LOW STRATOSPHERE
C****        7 EQUATORIAL LOW STRATOSPHERE  8 45 DEG NORTH LOW STRATOSPH
C****
C****        9 SOUTHERN MID STRATOSPHERE   10 NORTHERN MID STRATOSPHERE
C****       11 EQUATORIAL MID STRATOSPHERE 12 45 DEG NORTH MID STRATOSPH
C****
C****       13 SOUTHERN UPP STRATOSPHERE   14 NORTHERN UPP STRATOSPHERE
C****       15 EQUATORIAL UPP STRATOSPHERE 16 45 DEG NORTH UPP STRATOSPH
C****
      select case (m5)
      case (1,2,7,9,12,14,16) ! both winds and temperature have changed
        do_ke = .true.
        MKE=M5 ! mke,mape are not used for cases 1,2
        MAPE=M5+1
      case (11) ! radiation does not change winds
        do_ke = .false.
        MAPE=M5
      case default
        if(am_i_root()) then
          write(6,*) 'INCORRECT VALUE OF M5 IN DIAG5A.  M5 = ',m5
        endif
        call stop_model('INCORRECT VALUE OF M5 IN DIAG5A.',255)
      end select

C**** CURRENT TOTAL POTENTIAL ENERGY BY HEMISPHERE
c version w/o zatmo is computed elsewhere, so why the duplication?
      tpe_part(:) = 0d0
      do j=j_0,j_1
      do i=i_0,i_1
        sumt=0d0
        do l=1,lm
          sumt = sumt + t(i,j,l)*pk(l,i,j)*pdsig(l,i,j)
        end do
        sumt=(zatmo(i,j)*(p(i,j)+ptop)+sumt*sha)*axyp(i,j)
        if(lat2d(i,j).lt.0.) then ! southern hemisphere
          tpe_part(1) = tpe_part(1) + sumt
        else                      ! northern hemisphere
          tpe_part(2) = tpe_part(2) + sumt
        endif
      enddo
      enddo
      call sumxpe(tpe_part, tpe)

c
c For APE, first calculate global means for each layer of
c pot. temp (thgm) and static stability (gmean).  For convenience,
c use native-grid thermodynamic variables.  If interpolation
c to the latlon grid is not conservative, there will be
c a small inconsistency that is less important than interpretational
c ambiguities.
c
      do l=1,lm
        ldn=ldna(l)
        lup=lupa(l)
        thgm_part(l)=0d0
        gmean_part(l)=0d0
        do j=j_0,j_1
        do i=i_0,i_1
          thgm_part(l) = thgm_part(l)+t(i,j,l)*sqrtp(i,j)*axyp(i,j)
          gmean_part(l) = gmean_part(l) +axyp(i,j)*
     *         (p(i,j)*sig(l)+ptop)*(t(i,j,lup)-t(i,j,ldn))
     *         /(p(i,j)*pk(l,i,j))
        enddo
        enddo
      enddo
      call sumxpe(thgm_part,thgm)
      call sumxpe(gmean_part,gmean)
      call esmf_bcast(grid,thgm)
      call esmf_bcast(grid,gmean)
      do l=1,lm
        thgm(l)=thgm(l)/areag
        ldn=ldna(l)
        lup=lupa(l)
        gmean(l)=areag*(sig(ldn)-sig(lup))/gmean(l)
      enddo

c
c some arrays to make the kspher logic more robust
c
      klayer_eq(:) = klayer(:) + 2
      klayer_45n(:) = klayer(:) + 3

c
c Halo updates and fills of cube corners in preparation for interpolation.
c Some quantities may already have been haloed - check later.
c
      call halo_update(grid,p)
      call corner_fill_3D(grid,p,1)
      call halo_update(grid,t)
      call corner_fill_3D(grid,t,lm)
      call halo_update(grid,ualij,jdim=3)
      call corner_fill_4D(grid,ualij,lm,1)
      call halo_update(grid,valij,jdim=3)
      call corner_fill_4D(grid,valij,lm,1)

c
c Loop over latitudes
c
      ke_part(:,:)=0.
      ape_part(:,:)=0.

      do jj=1,jmlat

      jlat = cs2ll%jlat_sched(jj) ! the schedule facilitates parallelism
      if(cs2ll%ni(jlat).eq.0) cycle ! no work for this PE at this jlat

c
c interpolate p,t to this latitude circle
c
      call interp_to_jlat_3D(grid,cs2ll,t,tll ,lm,jlat)
      call interp_to_jlat_3D(grid,cs2ll,p,psll,1 ,jlat)
      do i=cs2ll%is(jlat),cs2ll%ie(jlat)
        if(psll(i).eq.0.) cycle ! valid lons sometimes noncontiguous
        sqrtpa(i) = sqrt(psll(i)*dxyp(jlat))
      enddo

c
c send temperature data to the root PE for this jlat
c
      do l=1,lm
      do i=cs2ll%is(jlat),cs2ll%ie(jlat)
        if(psll(i).eq.0.) then ! valid lons sometimes noncontiguous
          t_tmp(i,l) = 0. ! for pack_zonal
        else
          t_tmp(i,l)=tll(i,l)*sqrtpa(i)-thgm(l)
        endif
      enddo
      enddo
      call pack_zonal(cs2ll,jlat,lm,t_tmp,t_il)

      if(do_ke) then
c
c interpolate u,v to this latitude circle and
c send wind data to the root PE for this jlat
c
        call interp_to_jlat_4D(grid,cs2ll,ualij,uli,lm,1,jlat)
        call interp_to_jlat_4D(grid,cs2ll,valij,vli,lm,1,jlat)
        do l=1,lm
        do i=cs2ll%is(jlat),cs2ll%ie(jlat)
          if(psll(i).eq.0.) then ! valid lons sometimes noncontiguous
            u_tmp(i,l) = 0.; v_tmp(i,l) = 0. ! for pack_zonal
          else
            u_tmp(i,l) = uli(l,i)*sqrtpa(i)
            v_tmp(i,l) = vli(l,i)*sqrtpa(i)
          endif
        enddo
        enddo
        call pack_zonal(cs2ll,jlat,lm,u_tmp,u_il)
        call pack_zonal(cs2ll,jlat,lm,v_tmp,v_il)
      endif

c
c if at a jlat for calculation, do the spectral analysis
c
      jcalc = cs2ll%jlat_calc(jlat)
      if(jcalc.gt.0) then
        if(jcalc.gt.jmlat/2) then
          khemi = 2
        else
          khemi = 1
        endif
        do l=1,lm
          kspher = klayer(l) + khemi-1
          kspher_eq = klayer_eq(l)
          kspher_45n = klayer_45n(l)
c kinetic energy
          if(do_ke) then
            call ffte(u_il(:,l),xu)
            call ffte(v_il(:,l),xv)
            xke(:) = (xu(:)+xv(:))*dsig(l)
            ke_part(:,kspher) = ke_part(:,kspher) + xke(:)
            if(jcalc.eq.jlat_eq) then
              ke_part(:,kspher_eq) = ke_part(:,kspher_eq) + xke(:)
            elseif(jcalc.eq.jlat_45n) then
              ke_part(:,kspher_45n) = ke_part(:,kspher_45n) + xke(:)
            endif
          endif
c potential energy
          call ffte(t_il(:,l),xape)
          xape(:) = xape(:)*dsig(l)*gmean(l)
          ape_part(:,kspher) = ape_part(:,kspher) + xape(:)
          if(jcalc.eq.jlat_eq) then
            ape_part(:,kspher_eq) = ape_part(:,kspher_eq) + xape(:)
          elseif(jcalc.eq.jlat_45n) then
            ape_part(:,kspher_45n) = ape_part(:,kspher_45n) + xape(:)
          endif
        enddo
      endif

      enddo ! end latitude loop

c
c combine partial sums onto global root
c
      call sumxpe(ke_part,ke)
      call sumxpe(ape_part,ape)

c
c on global root, place sums and differences into SPECA
c
      if(am_i_root()) then

        if (ndt /= 0) then
c**** energy transfer rates as differences
          if(do_ke) then
            speca(:,mke,:)=speca(:,mke,:)+(ke(:,:)-speca(:,19,:))/ndt
          endif
          speca(:,mape,:)=speca(:,mape,:)+(ape(:,:)-speca(:,20,:))/ndt
          mtpe=mtpeof(mape)
          atpe(mtpe,:)=atpe(mtpe,:)+(tpe(:)-atpe(8,:))/ndt
        end if
c**** store latest values
        if(do_ke) speca(:,19,:)=ke(:,:)
        speca(:,20,:)=ape(:,:)
        atpe(8,:)=tpe(:)

        if (m5.eq.2) then
c**** accumulate mean kinetic energy and mean potential energy
          speca(:,2,:)=speca(:,2,:)+ke(:,:)
          speca(:,3,:)=speca(:,3,:)+ape(:,:)
          atpe(1,:)=atpe(1,:)+tpe(:)
        end if

      endif ! am_i_root

      call timer (now,mdiag)

      return
      end subroutine diag5a

      subroutine diag7a
!@sum DIAG7A calculate wave power around selected latitude circles.
!@+   This version is for a cubed sphere grid.
!@auth cube-enabled by M. Kelley
c****
c**** this routine accumulates a time sequence for selected
c**** quantities and from that prints a table of wave frequencies.
c****
      use constant, only : bygrav
      use model_com, only : lm,idacc,mdiag,p
      use dynamics, only : ualij,valij,phi
      use diag_com, only : nwav_dag,wave,max12hr_sequ
     &     ,kwp,re_and_im,ia_12hr
      use diag_loc, only : ldex
      use gcdiag
      use domain_decomp_atm, only : grid,sumxpe,am_i_root,
     &     halo_update
      use fftw_com
      implicit none

      real*8, dimension(0:imlonh) :: an,bn
      integer, parameter :: km=6,kqmax=12
      integer :: nmax=nwav_dag
      real*8, dimension(imlon,km) :: htrd,htrd_loc
      real*8, dimension(imlon,lm) :: ueq,ueq_loc,veq,veq_loc
      real*8, dimension(lm,cs2ll%isd:cs2ll%ied) :: uli,vli
      real*8, dimension(cs2ll%isd:cs2ll%ied,lm) :: uil,vil
      real*8, dimension(cs2ll%isd:cs2ll%ied,km) :: htrd_tmp
      real*8, dimension(cs2ll%isd:cs2ll%ied,lm) :: phill
      real*8, dimension(cs2ll%isd:cs2ll%ied) :: psll
      real*8, dimension(km), parameter ::
     &     pmb=(/922.,700.,500.,300.,100.,10./),
     &     ght=(/500.,2600.,5100.,8500.,15400.,30000./)
      real*8, dimension(lm) :: p00,aml,pdsigl,pmidl
      real*8, dimension(lm+1) :: pednl
      real*8 :: slope, now
      integer i,jlat,idacc9,k,kq,l,n

      idacc9=idacc(ia_12hr)+1
      idacc(ia_12hr)=idacc9
      if (idacc9.gt.max12hr_sequ) return

c
c Halo updates and fills of cube corners in preparation for interpolation.
c Some quantities may already have been haloed - check later.
c
      call halo_update(grid,p)
      call corner_fill_3D(grid,p,1)
      call halo_update(grid,phi)
      call corner_fill_3D(grid,phi,lm)
      call halo_update(grid,ualij,jdim=3)
      call corner_fill_4D(grid,ualij,lm,1)
      call halo_update(grid,valij,jdim=3)
      call corner_fill_4D(grid,valij,lm,1)

c
c interpolate winds to the equator and send to global root
c
      jlat = jlat_eq
      ueq_loc=0d0
      veq_loc=0d0
      if(cs2ll%ni(jlat).gt.0) then ! this processor has valid lons at jlat
        call interp_to_jlat_4D(grid,cs2ll,ualij,uli,lm,1,jlat)
        call interp_to_jlat_4D(grid,cs2ll,valij,vli,lm,1,jlat)
        do l=1,lm
          uil(:,l) = uli(l,:)
          vil(:,l) = vli(l,:)
        enddo
        call pack_zonal(cs2ll,jlat,lm,uil,ueq_loc)
        call pack_zonal(cs2ll,jlat,lm,vil,veq_loc)
      endif
      call sumxpe(ueq_loc, ueq)
      call sumxpe(veq_loc, veq)

c
c interpolate geopotential to 50 N and send to global root
c
      jlat = jlat_50n
      htrd_loc(:,:)=0d0
      if(cs2ll%ni(jlat).gt.0) then ! this processor has valid lons at jlat
        call interp_to_jlat_3D(grid,cs2ll,phi,phill,lm,jlat)
        call interp_to_jlat_3D(grid,cs2ll,p  ,psll ,1 ,jlat)
        do i=cs2ll%is(jlat),cs2ll%ie(jlat)
          if(psll(i).eq.0.) then ! valid lons sometimes noncontiguous
            htrd_tmp(i,:) = 0. ! for pack_zonal
            cycle
          endif
          call calc_vert_amp(psll(i),lm,p00,aml,pdsigl,pednl,pmidl)
          l=2
          do k=1,km
            do while(pmb(k).lt.pmidl(l) .and. l.lt.lm)
              l=l+1
            enddo
c**** assume that phi is linear in log p
            slope=(phill(i,l-1)-phill(i,l))/log(pmidl(l-1)/pmidl(l))
            htrd_tmp(i,k)=
     &           (phill(i,l)+slope*log(pmb(k)/pmidl(l)))*bygrav-ght(k)
          enddo
        enddo
        call pack_zonal(cs2ll,jlat,km,htrd_tmp,htrd_loc)
      endif
      call sumxpe(htrd_loc, htrd)

c
c do the fourier analyses on global root
c
      if(am_i_root()) then
        do kq=1,3
          call fft(ueq(1,ldex(kq)),an,bn)
          do n=1,nmax
            wave(1,idacc9,n,2*kq-1)=an(n)
            wave(2,idacc9,n,2*kq-1)=bn(n)
          enddo
          call fft(veq(1,ldex(kq)),an,bn)
          do n=1,nmax
            wave(1,idacc9,n,2*kq)=an(n)
            wave(2,idacc9,n,2*kq)=bn(n)
          enddo
        enddo
        do kq=7,kqmax
          call fft(htrd(1,kq-6),an,bn)
          do n=1,nmax
            wave(1,idacc9,n,kq)=an(n)
            wave(2,idacc9,n,kq)=bn(n)
          end do
        end do
      endif

      call timer (now,mdiag)
      return
      end subroutine diag7a

      subroutine diaggc_prep
      use model_com, only : lm
      use domain_decomp_atm, only : am_i_root
      use diag_com, only : kagc,agc,hemis_gc,vmean_gc
      use gcdiag
      implicit none
      integer :: j,j1,j2,k,l
      real*8 :: hemfac

      if(.not.am_i_root()) return

c stream function
      do j=1,jmlat
        agc(j,lm,gc_psi) = 0d0
        do l=lm-1,1,-1
          agc(j,l,gc_psi)=agc(j,l+1,gc_psi)-agc(j,l+1,gc_v)*dxp(j)
        enddo
      enddo

c
c compute hemispheric/global means and vertical sums
c
      hemfac = 2./sum(dxyp)
      do k=1,kagc
        do l=1,lm
          j1 = 1; j2 = jmlat/2
          hemis_gc(1,l,k) = hemfac*sum(agc(j1:j2,l,k)*dxyp(j1:j2))
          j1 = jmlat/2+1; j2 = jmlat
          hemis_gc(2,l,k) = hemfac*sum(agc(j1:j2,l,k)*dxyp(j1:j2))
          hemis_gc(3,l,k) = .5*(hemis_gc(1,l,k)+hemis_gc(2,l,k))
        enddo
        vmean_gc(1:jmlat,1,k) = sum(agc(:,:,k),dim=2)
        vmean_gc(jmlat+1:jmlat+3,1,k) = sum(hemis_gc(:,:,k),dim=2)
      enddo

      return
      end subroutine diaggc_prep
