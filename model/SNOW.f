!@sum contains code related to 3 layers snow model
!@auth I.Aleinov
!@ver  1.0
      MODULE SNOW_MODEL
!@sum SNOW_MODEL does column physics for 3 layers snow model
!@auth I.Aleinov
!@ver  1.0
      USE FILEMANAGER, only: openunit
      USE CONSTANT, only :
     $     rho_water => rhow,   ! (kg m-3)
     $     rho_ice => rhoi,     ! (kg m-3)
     $     lhm_kg => lhm,       ! (J kg-1)
     $     lhe_kg => lhe,       ! (J kg-1)
     $     shi_kg => shi,       ! (J kg-1 C-1)
     $     tfrz => tf,          ! (K)
     $     sigma => stbo        ! (W m-2 K-4)

      IMPLICIT NONE
      SAVE
      PRIVATE

      PUBLIC snow_adv,i_earth, j_earth

ccc physical parameters
!@var rho_fresh_snow density of fresh snow (kg m-3)
      real*8, parameter :: rho_fresh_snow =  150.d0 ! (kg m-3)
!@var lat_fusion volumetric latent heat of fusion for water  (J m-3)
      real*8, parameter :: lat_fusion = lhm_kg * rho_water ! (J m-3)
!@var lat_evap volumetric latent heat of evaporation for water  (J m-3)
      real*8, parameter :: lat_evap = lhe_kg * rho_water ! (J m-3)
!@var max_fract_water max fraction of free water in a wet snow (1)
      real*8, parameter :: max_fract_water = .055d0

ccc model parameters
!@var EPS small number (used as cutoff parameter) (1)
      real*8, parameter :: EPS = 1.d-8  ! was 1.d-12
!@var MAX_NL maximal number of snow layers (actually only 3 are used) (1)
      integer, parameter :: MAX_NL = 16
!@var MIN_SNOW_THICKNESS minimal thickness of snow (m)
ccc trying to increase MIN_SNOW_THICKNESS for stability
ccc      real*8, parameter :: MIN_SNOW_THICKNESS =  0.01d0  ! was 0.09d0
      real*8, parameter :: MIN_SNOW_THICKNESS =  0.1d0
!@var MIN_FRACT_COVER miminal areal fraction of snow cover (cutoff param) (1)
      real*8, parameter :: MIN_FRACT_COVER = 0.0001d0

!@var DEB_CH channel for debug output
      integer, save :: DEB_CH = 0
!@var i_earth, j_earth coordinate of current point (for debugging)
      integer i_earth, j_earth

      CONTAINS

      subroutine pass_water( wsn, hsn, dz, nl,
     &                       water_down, heat_down,
     &                       lat_fusion, max_fract_water,
     &                       rho_water, rho_fresh_snow)
!@sum  pass extra water to the lower layers
!@auth I.Aleinov
!@ver  1.0
      implicit none
      integer nl
      real*8  wsn(nl), hsn(nl), dz(nl)
      real*8 water_down, heat_down, lat_fusion, max_fract_water
      real*8 rho_water, rho_fresh_snow
      integer n
      real*8 ice, free_water, ice_old
      integer pass_water_flag

      pass_water_flag = 0
      do n=1,nl
        ice_old = min(wsn(n),-hsn(n)/lat_fusion)
        ice = ice_old
        wsn(n) = wsn(n) + water_down
        hsn(n) = hsn(n) + heat_down
        water_down = 0.d0
        heat_down = 0.d0
        if( hsn(n) .ge. 0.d0 .or. wsn(n) .le. 0.d0 ) then ! all melted
          water_down = wsn(n)
          heat_down = hsn(n)
          wsn(n) = 0.d0
          hsn(n) = 0.d0
          dz(n) = 0.d0
          pass_water_flag = 1
        else if (hsn(n).gt.-wsn(n)*lat_fusion) then !/* some melted */
          ice = -hsn(n)/lat_fusion
          free_water = wsn(n) - ice
c!!!                may be ice*max_fract_water ??
          water_down = max(0d0, free_water - wsn(n)*max_fract_water)
          wsn(n) = wsn(n) - water_down
          !/* should I reduce dz here ??? */
          !/* dz(n) = dz(n) * min( ice/ice_old, 1.d0 )  */
          dz(n) = min( dz(n), ice*rho_water/rho_fresh_snow )
          pass_water_flag = 2
        else !/* all frozen */
          if (wsn(n)+EPS .lt. ice_old) dz(n) = dz(n)*wsn(n)/ice_old
          dz(n) = min( dz(n), wsn(n)*rho_water/rho_fresh_snow )
          endif
        enddo
      return
      end subroutine pass_water


      subroutine snow_redistr(dz, wsn, hsn, nl,
     &                         dzo, wsno, hsno, nlo, fract_cover)
!@sum  redistributes snow between the layers
!@auth I.Aleinov
!@ver  1.0
      implicit none
      integer nl,nlo
      real*8 dz(nl+1), wsn(nl), hsn(nl)
      real*8 dzo(nl+1), wsno(nl), hsno(nl)
      real*8 fract_cover
      integer n, no
      real*8 total_dz, ddz, delta, fract
      real*8 fract_cover_new, fract_cover_ratio

      total_dz = 0.d0
      do n=1,nlo
        total_dz = total_dz + dzo(n)
        enddo

      fract_cover_new =
     &          min( 1.d0, total_dz*fract_cover/MIN_SNOW_THICKNESS )

      if( fract_cover_new .lt. MIN_FRACT_COVER ) then
        fract_cover = 0.d0
        nl = 1
        dz(1) = 0.d0
        wsn(1) = 0.d0
        hsn(1) = 0.d0
        return
        endif

      fract_cover_ratio = fract_cover/fract_cover_new
      total_dz = total_dz*fract_cover_ratio

      if ( total_dz .gt. MIN_SNOW_THICKNESS*6.d0 ) then
        nl = 3
      else
        nl = 1
        endif

c!!!  trying one layer model ...
c!!      nl = 1

      if ( total_dz .gt. 0.16d0 .and. nl.gt.1 ) then
        dz(1) = 0.05d0   ! trying thinner upper layer
c!!        dz(1) = 0.10d0
        ddz = (total_dz-dz(1))/(nl-1)
        do n=2,nl
          dz(n) = ddz
          enddo
      else
        ddz = total_dz/nl
        do n=1,nl
          dz(n) = ddz
          enddo
        endif

      do n=1,nl
        wsn(n) = 0.d0
        hsn(n) = 0.d0
        enddo

      no = 0
      delta = 0.d0
      do n=1,nl
        do while ( delta.lt.dz(n) .and. no.lt.nlo )
          no = no + 1
          delta = delta + dzo(no)*fract_cover_ratio
          wsn(n) = wsn(n) + wsno(no)*fract_cover_ratio
          hsn(n) = hsn(n) + hsno(no)*fract_cover_ratio
          enddo
        ddz = delta - dz(n)
        fract = ddz/(dzo(no)*fract_cover_ratio)
ccc the following is just for check
          if ( fract.lt.-EPS .or. fract.gt.1.d0+EPS ) then
            print *, 'internal error 3 in snow_redistr'
            print *, 'fract= ', fract
            call abort
            stop 251
            endif
        wsn(n) = wsn(n) - fract*wsno(no)*fract_cover_ratio
        hsn(n) = hsn(n) - fract*hsno(no)*fract_cover_ratio
        if ( n.lt.nl ) then
          wsn(n+1) = wsn(n+1) + fract*wsno(no)*fract_cover_ratio
          hsn(n+1) = hsn(n+1) + fract*hsno(no)*fract_cover_ratio
          delta = ddz
        else
          if ( abs(fract).gt.EPS ) then
            print *,'internal error 1 in snow_redistr'
            print *, 'fract= ', fract
            call abort
            stop 251
            endif
          endif
        enddo

      fract_cover = fract_cover_new
      return
      end subroutine snow_redistr


      subroutine snow_adv(dz, wsn, hsn, nl,
     &    srht, trht, snht, htpr, evaporation, pr, dt,
     &    t_ground, dz_ground, fract_cover,
     &    tsn_surf, water_to_ground, heat_to_ground,
     &    radiation_out, snsh_dt, evap_dt )
      implicit none
!@sum  a wrapper that calles real snow_adv (introduced for debugging)
!@auth I.Aleinov
!@ver  1.0
ccc input:
      integer nl
      real*8 srht, trht, snht, htpr, evaporation
      real*8 pr, dt, t_ground, dz_ground, fract_cover
      real*8 snsh_dt, evap_dt
ccc output:
      real*8 water_to_ground, heat_to_ground, tsn_surf
      real*8 radiation_out
ccc data arrays
      real*8 dz(nl+1), wsn(nl), hsn(nl)


      COMMON /SOILIN/ PR_i,PRS_i,HTPR_i,HTPRS_i,
     &  PRES_i,TS_i,QS_i,  CH_i,VSM_i,RHO_i,
     *  SRHT_i,TRHT_i,COST_i,SINT_i,EDDY_i,Q1_i,QM1_i,Z1_i,ZS_i,
     &  IGCM_i
      real*8 PR_i,PRS_i,HTPR_i,HTPRS_i
      real*8 PRES_i,TS_i,QS_i,  CH_i,VSM_i,RHO_i
      real*8 SRHT_i,TRHT_i,COST_i,SINT_i,EDDY_i,Q1_i,QM1_i,Z1_i,ZS_i
      integer IGCM_i

ccc for debug
      real*8 total_energy  !, lat_evap
      integer i, retcode
      common /snowmodel_debug/ hsn_old, dz_old, fract_cover_old,
     &     wsn_old, retcode
      real*8 hsn_old, dz_old, fract_cover_old, wsn_old

ccc for debug
      ! lat_evap =        2.5d9    !/* J m-3  */ defined globally
      hsn_old = hsn(1)
      wsn_old = wsn(1)
      dz_old = dz(1)
      fract_cover_old = fract_cover
      total_energy = 0.d0
      do i=1,nl
         total_energy = total_energy - hsn(i)*fract_cover
         enddo

      call snow_adv_1(dz, wsn, hsn, nl,
     &    srht, trht, snht, htpr, evaporation, pr, dt,
     &    t_ground, dz_ground, fract_cover,
     &    tsn_surf, water_to_ground, heat_to_ground,
     &    radiation_out, snsh_dt, evap_dt, retcode )

ccc for debug
      do i=1,nl
         total_energy = total_energy + hsn(i)*fract_cover
         enddo
      total_energy = total_energy -
     &    (srht+trht-snht-lat_evap*evaporation+htpr)*dt
      total_energy = total_energy + heat_to_ground + radiation_out*dt
      if ( abs(total_energy) .gt. 1.d0 ) then
        print*, "total energy error",i_earth, j_earth,total_energy,
     *       heat_to_ground,radiation_out*dt
         stop 'total energy error' ! call abort
       end if
ccc   just for debug:
ccc      water_to_ground =  pr
ccc      heat_to_ground = (srht+trht)*dt
      return
      end subroutine snow_adv


      subroutine snow_adv_1(dz, wsn, hsn, nl,
     &    srht, trht, snht, htpr, evaporation, pr, dt,
     &    t_ground, dz_ground, fract_cover,
     &    tsn_surf, water_to_ground, heat_to_ground,
     &    radiation_out, snsh_dt, evap_dt, retcode )
      implicit none
!@sum main program that does column snow physics
!@auth I.Aleinov
!@ver  1.0
ccc input:
      integer nl
      real*8 srht, trht, snht, htpr, evaporation
      real*8 pr, dt, t_ground, dz_ground, fract_cover
      real*8 snsh_dt, evap_dt

ccc constants: (now defined as global params)
      real*8 k_ground, c_ground
ccc output:
      real*8 water_to_ground, heat_to_ground, tsn_surf
      real*8 radiation_out
      integer retcode

ccc main parameters: layer thickness, water equivalent, heat content
      real*8 dz(nl+1), wsn(nl), hsn(nl)
      real*8 dzo(MAX_NL+1), wsno(MAX_NL), hsno(MAX_NL)
      real*8 tsn(MAX_NL+1), csn(MAX_NL), ksn(MAX_NL+1), isn(MAX_NL)

      real*8 water_down, heat_down, fresh_snow, rho_snow, flux_in
      real*8 mass_layer, mass_above, scale_rho
      real*8 fract_cover_old, flux_in_deriv, flux_corr
      real*8 delta_tsn_impl ! shft of temperature due to implicit method
      integer have_snow, n, nlo

ccc common for debug

ccc!!!  check if lat_evap shoud be replaced by heat of sublimation

ccc !!! ground properties should be passed as formal parameters !!!
      k_ground =        3.4d0    !/* W K-1 m */    /* --??? */
      c_ground =        1.d5     !/* J m-3 K-1 */  /* -- ??? */

ccc fluxes in/out
      heat_to_ground = 0.d0
      water_to_ground = 0.d0

c!!!  this is for debugging
      do n=1,nl
        if(fract_cover .gt. EPS) then
          if(wsn(n).gt.EPS .and.
     &       wsn(n)/dz(n)*rho_water+EPS.lt.rho_fresh_snow) then
            print*,"wsn error 1",i_earth, j_earth,n,wsn(n),
     *       wsn(n)/dz(n)*rho_water
     *           ,rho_fresh_snow
            stop 'wsn error 1' ! call abort
          end if
          endif
        enddo

ccc the following lines fix the problem with initial call
ccc with fract_cover==1 and dz(1)<MIN_SNOW_THICKNESS
      if( fract_cover .ge. 1.d0-EPS .and. nl .eq. 1
     &    .and. dz(1) .lt. MIN_SNOW_THICKNESS ) then
ccc condition nl .eq. 1  was put to the statement above
ccc the following if should be removed if it works ok with thicker snow
        if( nl .ne. 1 ) then
          print *, 'OOPS: nl= ',nl,' fract_cover= ',fract_cover
     &              ,'dz= ', dz
          stop 251
          endif
        fract_cover = dz(1)/MIN_SNOW_THICKNESS
ccc        if( fract_cover .lt. EPS ) then
c!!!  use fract_cover == 1 for debug only !!!!
        if( fract_cover .lt. MIN_FRACT_COVER ) then
          fract_cover = 0.d0
          heat_to_ground = hsn(1)
          water_to_ground = wsn(1)
          dz(1) = 0.d0
          hsn(1) = 0.d0
          wsn(1) = 0.d0
          else
          dz(1) = MIN_SNOW_THICKNESS
          hsn(1) = hsn(1)/fract_cover
          wsn(1) = wsn(1)/fract_cover
          endif
        endif

ccc compute amount of fresh snow
ccc !!! insert evaporation into computation of the amount of fresh snow
      fresh_snow = rho_water/rho_fresh_snow
     &                * min(pr*dt-evaporation*dt, -htpr*dt/lat_fusion)
      if(fresh_snow.lt.0.d0) fresh_snow = 0.d0

      if( fract_cover .lt. MIN_FRACT_COVER .and.
     & fresh_snow.lt.MIN_SNOW_THICKNESS*MIN_FRACT_COVER) then ! no snow
        fract_cover = 0.d0
        tsn_surf = t_ground
        heat_to_ground = heat_to_ground
     &      +(srht+trht-sigma*(tsn_surf+tfrz)**4
     &      -lat_evap*evaporation-snht+htpr)*dt
        radiation_out = sigma*(tsn_surf+tfrz)**4
        water_to_ground = water_to_ground + (pr-evaporation)*dt
        retcode = 1
        return
      endif

      if ( fract_cover.lt.1.d0 .and.
     &   fresh_snow .lt. MIN_SNOW_THICKNESS*(1.d0-fract_cover) ) then
ccc partial cover
        hsn(1) = hsn(1)*fract_cover
        wsn(1) = wsn(1)*fract_cover
        fract_cover =
     &           (dz(1)*fract_cover + fresh_snow)/MIN_SNOW_THICKNESS

        if (fract_cover.gt.1.d0 .or. fract_cover.lt.EPS) then
          if ( DEB_CH == 0 )
     $         call openunit("snow_debug", DEB_CH, .false., .false.)
          write(DEB_CH,*) 'EROOR: fract_cover= ', fract_cover
          write(DEB_CH,*) dz(1), fresh_snow
          if ( fract_cover.gt.1.d0 ) fract_cover = 1.d0
          if ( fract_cover.lt.EPS ) fract_cover = EPS
        endif

        dz(1) = MIN_SNOW_THICKNESS
        hsn(1) = hsn(1)/fract_cover
        wsn(1) = wsn(1)/fract_cover
ccc all precipitation falls on the snow
        water_down = pr*dt/fract_cover
        heat_down = htpr*dt/fract_cover
      else
ccc full cover
        hsn(1) = hsn(1)*fract_cover
        wsn(1) = wsn(1)*fract_cover
        dz(1) = dz(1)*fract_cover
        fract_cover = 1.d0
        dz(1) = dz(1) + fresh_snow
        water_down = pr*dt
        heat_down = htpr*dt
      endif

ccc !!! subtract evaporation somewhere here
      water_down = water_down  - evaporation*dt
c!!      heat_down = heat_down - lat_evap*evaporation*dt

      water_to_ground = water_to_ground -
     &                  evaporation*dt*(1.d0-fract_cover)

c!!!  this is for debugging
      do n=1,nl
        if((wsn(n)+water_down+evaporation*dt)/dz(n)*rho_water+EPS
     &      .lt.rho_fresh_snow) then
          print*,"wsn error 2",i_earth, j_earth,n,wsn(n),(wsn(n)+
     *       water_down+evaporation*dt)/dz(n)*rho_water,rho_fresh_snow
          stop 'wsn error 2' ! call abort
        end if
      enddo

      call pass_water( wsn, hsn, dz, nl, water_down, heat_down,
     &                       lat_fusion, max_fract_water,
     &                       rho_water, rho_fresh_snow)
      heat_to_ground = heat_to_ground + heat_down*fract_cover
      water_to_ground = water_to_ground + water_down*fract_cover

c!!!  this is for debugging
      do n=1,nl
        if(wsn(n).gt.EPS .and.
     &     wsn(n)/dz(n)*rho_water+EPS.lt.rho_fresh_snow) then
          print*,"wsn error 3",i_earth, j_earth,n,wsn(n),
     *    wsn(n)/dz(n)*rho_water
     *         ,rho_fresh_snow
          stop 'wsn error 3' ! call abort
        end if
      enddo

ccc redistribute snow over the layers
      fract_cover_old = fract_cover
      nlo = nl
      do n=1,nl
        dzo(n) = dz(n)
        hsno(n) = hsn(n)
        wsno(n) = wsn(n)
      enddo

      call snow_redistr(dz, wsn, hsn, nl,
     &                         dzo, wsno, hsno, nlo, fract_cover)
      dz(nl+1) = dz_ground

ccc check if there is any snow at all ?
ccc      if ( fract_cover.lt.EPS ) then   !/* no snow left ... */
c!!! disable fractional cover for debugging
      if ( fract_cover.lt.MIN_FRACT_COVER ) then   !/* debugging !!!*/
        tsn_surf = t_ground
        heat_to_ground = heat_to_ground +
     &       ( srht+trht-sigma*(tsn_surf+tfrz)**4
     &      -lat_evap*evaporation-snht )*dt
        radiation_out = sigma*(tsn_surf+tfrz)**4
        do n=1,nlo
           heat_to_ground = heat_to_ground + hsno(n)*fract_cover_old
           water_to_ground = water_to_ground + wsno(n)*fract_cover_old
        enddo
        retcode = 2
        return
      endif

ccc compute spec. heat and thermal conductivity
      do n=1,nl
        rho_snow = wsn(n)*rho_water/dz(n)
        !csn(n) = (1.9d6) * rho_snow / rho_ice
        csn(n) = shi_kg * rho_snow
        ksn(n) = (3.22d-6) * rho_snow**2
      enddo
      ksn(nl+1) = k_ground

ccc compute temperature of the layers (and amount of ice)
      do n=1,nl
        if ( hsn(n) .gt. 0.d0 ) then
          print *, 'OOPS, No snow in layer ', n
          stop 251
        else if ( hsn(n) .gt. -wsn(n)*lat_fusion ) then
          tsn(n) = 0.d0
          isn(n) = -hsn(n)/lat_fusion
        else
          tsn(n) = (hsn(n)+wsn(n)*lat_fusion)/(csn(n)*dz(n))
          isn(n) = wsn(n)
        endif
      enddo
      tsn(nl+1) = t_ground

c!!! this is for debugging
c      if(tsn(1).lt.-120.d0) call abort

ccc compute incomming heat flux (from atm.)
ccc include all fluxes except htpr (which is already included)
      flux_in =
     &      srht
     &      +trht - sigma*(tsn(1)+tfrz)**4
     &      -lat_evap*evaporation
     &      -snht

      flux_in_deriv =
     &     -4.d0*sigma*(tsn(1)+tfrz)**3 - lat_evap*evap_dt - snsh_dt

ccc      heat_to_ground = heat_to_ground + flux_in*dt*(1.d0-fract_cover)
      heat_to_ground = heat_to_ground + (1.d0-fract_cover)*dt
     & * ( srht+trht-sigma*(t_ground+tfrz)**4
     &   -lat_evap*(evaporation+evap_dt*(t_ground-tsn(1)))
     &   -snht-snsh_dt*(t_ground-tsn(1)) )

      radiation_out = fract_cover*sigma*(tsn(1)+tfrz)**4
     &   +  (1.d0-fract_cover)
     & * ( sigma*(t_ground+tfrz)**4
     &   +lat_evap*evap_dt*(t_ground-tsn(1))
     &   +snsh_dt*(t_ground-tsn(1)) )

ccc solve heat equation

      call heat_eq(
     &     dz, tsn, hsn, csn, ksn, nl,
     &     flux_in, flux_in_deriv, flux_corr, dt )

      heat_to_ground = heat_to_ground + flux_in*dt*fract_cover

ccc now redistribute extra flux proportionallly to input derivatives
ccc snht_out_cor, delta_tsn_impl
      delta_tsn_impl = flux_corr / flux_in_deriv
ccc      radiation_out = radiation_out - flux_corr*fract_cover
      radiation_out = radiation_out -
     &     (-4.d0*sigma*(tsn(1)+tfrz)**3)*delta_tsn_impl*fract_cover
      snht = snht + snsh_dt * delta_tsn_impl * fract_cover
      evaporation = evaporation + evap_dt * delta_tsn_impl * fract_cover

ccc and now remove (add) water due to extra evaporation.
c!!! this may make wsn(1) negative, the only way I see now to prevent it
c!!! is to keep minimal thickness of snow big enough

      water_down = - evap_dt * delta_tsn_impl

c!!!  this is for debugging
      do n=1,nl
        if(wsn(n).gt.EPS .and.
     &     wsn(n)/dz(n)*rho_water+EPS.lt.rho_fresh_snow) then
          print*,"wsn error 4",i_earth, j_earth,n,wsn(n),
     *    wsn(n)/dz(n)*rho_water
     *         ,rho_fresh_snow
          stop 'wsn error 4' ! call abort
        end if
      enddo

ccc pass extra water down
ccc      water_down = 0.d0
      heat_down = 0.d0
      call pass_water( wsn, hsn, dz, nl, water_down, heat_down,
     &                       lat_fusion, max_fract_water,
     &                       rho_water, rho_fresh_snow)
      heat_to_ground = heat_to_ground + heat_down*fract_cover
      water_to_ground = water_to_ground + water_down*fract_cover

ccc update dz
      do n=1,n
        if( hsn(n).gt.0.d0 .or. isn(n).lt.EPS ) then
          dz(n) = 0.d0
        else if( hsn(n) .gt. -wsn(n)*lat_fusion ) then
          dz(n) = dz(n) * min( (-hsn(n)/lat_fusion)/isn(n), 1.d0)
        endif
      enddo

ccc redistribute snow over the layers
      fract_cover_old = fract_cover
      nlo = nl
      do n=1,nl
        dzo(n) = dz(n)
        hsno(n) = hsn(n)
        wsno(n) = wsn(n)
      enddo

      call snow_redistr(dz, wsn, hsn, nl,
     &                         dzo, wsno, hsno, nlo, fract_cover)
      dz(nl+1) = dz_ground

ccc check if there is any snow at all ?
      if ( fract_cover.lt.MIN_FRACT_COVER ) then   !/* debudgging !!!*/
        tsn_surf = t_ground
        do n=1,nl
           heat_to_ground = heat_to_ground + hsno(n)*fract_cover_old
           water_to_ground = water_to_ground + wsno(n)*fract_cover_old
        enddo
        retcode = 3
        return
      endif

ccc compute temperature of the layers
      do n=1,nl
        if ( hsn(n) .gt. 0.d0 ) then
          print *, 'OOPS, No snow in layer ', n
          stop 251
        else if ( hsn(n) .gt. -wsn(n)*lat_fusion ) then
          tsn(n) = 0.d0
        else
          tsn(n) = (hsn(n)+wsn(n)*lat_fusion)/(csn(n)*dz(n))
        endif
      enddo
      tsn(nl+1) = t_ground

ccc repack the layers
      mass_above = 0.d0
      do n=1,nl
        if( dz(n) .gt. EPS ) then
          mass_layer = wsn(n)*rho_water
          mass_above = mass_above + .5d0 * mass_layer
          scale_rho = .5d-7 * 9.8d0 * mass_above
     &     * exp(14.643d0 - 4000.d0/(tsn(n)+tfrz)
     &                  - .02d0 * mass_layer/dz(n))
     &     * dt
          scale_rho = 1.d0 + scale_rho
          dz(n) = dz(n) / scale_rho
          if( dz(n).lt.mass_layer/rho_ice )
     &                   dz(n) = mass_layer/rho_ice
          mass_above = mass_above + .5d0 * mass_layer
        endif
      enddo

      tsn_surf = tsn(1)

c!!!  this is for debugging
      do n=1,nl
        if(wsn(n).gt.EPS .and.
     &     wsn(n)/dz(n)*rho_water+EPS.lt.rho_fresh_snow) then
          print*,"wsn error 5",i_earth, j_earth,n,wsn(n),
     *    wsn(n)/dz(n)*rho_water
     *         ,rho_fresh_snow
          stop 'wsn error 5' ! call abort
        end if
      enddo

c!!! this is for debugging
ccc      if(tsn(1).lt.-120.d0) call abort

      if ( i_earth.eq.6300 .and. j_earth.eq.3200) then
        if ( DEB_CH == 0 )
     $       call openunit("snow_debug", DEB_CH, .false., .false.)
         write(DEB_CH,*) tsn(1), hsn(1), wsn(1), dz(1), fract_cover, nl
      endif

c!!! this is for debugging
      if(tsn(1).lt.-120.d0) then
        print*,"tsn error",i_earth, j_earth,1,tsn(1:nl)
        stop 'tsn error' ! call abort
      end if

      retcode = 0
      return
      end subroutine snow_adv_1


      subroutine heat_eq(
     &     dz, tsn, hsn, csn, ksn, nl,
     &     flux_in, flux_in_deriv, flux_corr, dt)
      implicit none
!@sum solves heat transport equation for snow layers
!@auth I.Aleinov
!@ver  1.0
!@alg This solver is currently using a half-implicit scheme.
!@+   Implicitness is controlled by the parameters:
!@+       gamma - for upper boundary
!@+       eta(MAX_NL+1) - for the rest of the boundaries
!@+   In general we are trying to use half-implicit method (gamma = .5d0)
!@+   at the surface. But this may introduce a systematic error when
!@+   temperature of the snow is close to 0 C. 
!@+   When option DO_EXPLIC_0 is enabled the program checks for possible
!@+   error and if necessary recomputes the solution with reduced gamma,
!@+   i.e. in more explicit way.
!@+   DO_EXPLIC_0 is enabled by default. One may try to turn it off if
!@+   there are serious stability problems.
!@calls :TRIDIAG
      integer nl
      real*8 dz(nl+1), tsn(nl+1), hsn(nl)
      real*8 csn(nl), ksn(nl+1)
ccc n+1 corresponds to the first layer of the ground
      real*8 flux_in, flux_in_deriv, flux_corr, dt, syst_flux_err

      real*8 eta(MAX_NL+1), tnew(MAX_NL), gamma
ccc arrays for coefficients:
      real*8 a(MAX_NL), b(MAX_NL), c(MAX_NL), f(MAX_NL)
      real*8 left, right, dt_to_cdz
      integer n, iter
      integer :: itermax=2

ccc setting implicit/explicit choice for each layer
ccc eta=1 - implicit, eta=0 - explicit
ccc for now using explicit method for layers with T=0.d0

      do n=1,nl
         eta(n) = .5d0
ccc        eta(n) = 1.d0
c!!        if( tsn(n) .ge. 0.d0 ) eta(n) = 0.d0
        enddo
      eta(nl+1) = 0.d0
ccc      gamma = .5d0
c!! I am trying to make it nearly expicit for a test !!!
      gamma = .5d0

ccc In general we are trying to use half-implicit method (gamma = .5d0)
ccc on the interface. But this introduces a systematic error when
ccc tnew(1) > 0. So in the case DO_EXPLIC_0 we move to explicit method
ccc (i.e. reduce gamma)

ccc DO_EXPLIC_0 enables correction of the error which is introduced
ccc by implicit scheme when tsn(1) == 0 C
ccc (of course one has to use preprocessing to do it ...)
ccc it is enabled by default

#define DO_EXPLIC_0
#ifdef DO_EXPLIC_0
      do iter=1,itermax
#endif
ccc equation for upper snow layer
        n = 1
        dt_to_cdz = dt/(csn(n)*dz(n))
        right = 2.d0*dt_to_cdz/( dz(n)/ksn(n) + dz(n+1)/ksn(n+1) )
        a(n) = 1.d0 + right*eta(n) - dt_to_cdz*flux_in_deriv*gamma
        b(n) = 0.d0
        c(n) = -right*eta(n+1)
        f(n) = tsn(n)*( 1.d0 - right*(1.d0-eta(n))
     &                  - dt_to_cdz*flux_in_deriv*gamma )
     &         + tsn(n+1)*right*(1.d0-eta(n+1))
     &         + dt_to_cdz * flux_in

ccc all other snow layers including the lower one
        do n=2,nl

          dt_to_cdz = dt/(csn(n)*dz(n))
          right = 2.d0*dt_to_cdz/( dz(n)/ksn(n) + dz(n+1)/ksn(n+1) )
          left  = 2.d0*dt_to_cdz/( dz(n)/ksn(n) + dz(n-1)/ksn(n-1) )

          a(n) = 1.d0 + ( left + right )*eta(n)
          b(n) = -left*eta(n-1)
          c(n) = -right*eta(n+1)

          f(n) = tsn(n)*( 1.d0 - ( left + right )*(1.d0-eta(n)) )
     &         + tsn(n-1)*left*(1.d0-eta(n-1))
     &         + tsn(n+1)*right*(1.d0-eta(n+1))
          enddo

c        call sweep3diag(b, c, a, f, tnew, nl)
          call TRIDIAG(b,a,c,f,tnew,nl)

ccc flux_corr is the energy wich should be returned to the atmosphere
        flux_corr = flux_in_deriv*( tnew(1) - tsn(1) )*gamma

#ifdef DO_EXPLIC_0
        syst_flux_err = flux_in_deriv*( tnew(1) - 0.d0 )*gamma
        ! if( iter/=2 .and. tnew(1)>0.d0 .and. flux_in_deriv<0.d0 ) then
        ! back to 77 :-L
        if ( iter.ne.itermax .and.
     &         tnew(1).gt.0.d0 .and. flux_in_deriv.lt.0.d0 ) then
          syst_flux_err = flux_in_deriv*( tnew(1) - 0.d0 )*gamma
          gamma = (1.d0 - syst_flux_err/flux_corr) *gamma
        else
          ! exit
          ! back to 77 :-L
          goto 77
        endif
      enddo
 77   continue
#endif

      do n=1,nl
        hsn(n) = hsn(n) + (tnew(n)-tsn(n))*csn(n)*dz(n)
        enddo

ccc flux to the ground :
      n = nl
      flux_in = - ( tsn(n+1)-tnew(n)*eta(n)-tsn(n)*(1.d0-eta(n)) ) *
     &           2.d0/( dz(n)/ksn(n) + dz(n+1)/ksn(n+1) )

      return
      end subroutine heat_eq

      END MODULE SNOW_MODEL

