 
      module snow_drvm

      contains
      subroutine snow_drv(
ccc input
     &     fm, evap, snsh, srht, trht, canht,
     &     drips, dripw, htdrips, htdripw,
     &     devap_dt, dsnsh_dt, dts,
     &     tp_soil, dz_soil, nlsn,
ccc updated
     &       dzsn, wsn, hsn, nsn,
     &       fr_snow,
ccc output
     &     flmlt, fhsng, flmlt_scale, fhsng_scale, thrmsn,
     &     flux_snow
     &     )

      use snow_model, only: snow_adv, snow_fraction, snow_redistr
      implicit none

ccc input
!@var canht heat flux from canopy, i.e. stbo*t_can**4 (w/m**2)
      real*8 fm, evap, snsh, srht, trht, canht
      real*8 drips, dripw, htdrips, htdripw
      real*8 devap_dt, dsnsh_dt, dts
      real*8 tp_soil, dz_soil
      integer nlsn
ccc updated
      real*8 dzsn(nlsn+1), wsn(nlsn), hsn(nlsn), fr_snow
      integer nsn
ccc output
      real*8 flmlt, fhsng, flmlt_scale, fhsng_scale, thrmsn
      real*8 flux_snow(0:nlsn)

ccc  local vars:

      real*8 epotsn, srhtsn, trhtsn, snshsn
      real*8 prsn, htprsn
ccc      real*8 xkthsn(2),cthsn(2)

      real*8 devap_sn_dt, dsnsh_sn_dt
      real*8 fr_snow_old,epot_sn_old,tr_flux(0:nlsn)


      epotsn = fm*evap
      snshsn = fm*snsh
      srhtsn = fm*srht
      trhtsn = fm*trht + (1.d0-fm)*canht! i.e. thrm_can

ccc!!! xkthsn cthsn are actually not used at the moment - should include
c!!! should pass ground properties to snow_adv
c!      xkthsn(1)=xkh(1,1)
c!      xkthsn(2)=xkh(1,2)
c!      cthsn(1)=shc(1,1)
c!      cthsn(2)=shc(1,2)

      devap_sn_dt = fm*devap_dt
      dsnsh_sn_dt = fm*dsnsh_dt

    !  do ibv=i_bare,i_vege
      fr_snow_old = fr_snow
      epot_sn_old = epotsn

      call snow_fraction(dzsn, nsn, drips, dts,
     &     fr_snow_old, fr_snow )

      if ( fr_snow <= 0.d0 ) then ! no snow
        flmlt_scale = drips
        fhsng_scale = htdrips
        flmlt = 0.d0
        fhsng = 0.d0
        thrmsn = 0.d0           ! just in case, is * 0 anyway
        flux_snow(:) = 0.d0
        if ( fr_snow_old > 0.d0 ) then ! had snow on prev. step
          flmlt_scale  = flmlt_scale
     &         + sum( wsn(1:nsn) )*fr_snow_old/dts
          fhsng_scale = fhsng_scale
     &         + sum( hsn(1:nsn) )*fr_snow_old/dts
          wsn (1:nsn) = 0.d0
          hsn (1:nsn) = 0.d0
          dzsn(1:nsn) = 0.d0
          fr_snow = 0.d0
          nsn = 1
        endif
        return  ! nothing else to do ...
      endif
         
ccc have snow - process it

      tr_flux = 0
      call snow_redistr( dzsn, wsn, hsn,
     &     nsn, fr_snow_old/fr_snow, tr_flux, dts )

ccc for tracers
      flux_snow(:) = tr_flux(:)

      ! all snow falls on fr_snow, but all liquid H2O uniformly
      prsn = drips/fr_snow + dripw
      htprsn = htdrips/fr_snow + htdripw

!!! debug : check if  snow is redistributed correctly
      if ( dzsn(1) > 0.d0 .and.
     &     dzsn(1) + prsn*dts*100.d0/15.d0 < .099d0) then
        call stop_model("set_snow: error in dz",255)
      endif

      call snow_adv(dzsn, wsn, hsn, nsn,
     &     srhtsn, trhtsn, snshsn, htprsn,
     $     epotsn,
     $     prsn, dts,
     &     tp_soil, dz_soil,
     &     flmlt, fhsng,
     &     thrmsn, dsnsh_sn_dt, devap_sn_dt, ! fbfv(ibv),
     &     tr_flux )

      flux_snow(:) = flux_snow(:) + tr_flux(:)
      flmlt_scale = 0.d0
      fhsng_scale = 0.d0

ccc make sure we don''t get negative 
      if ( flmlt < -1.d-16 ) then
        print *,'snow_drv: melt water flux < 0'
        print *,'flmlt  = ',flmlt
        print *,'pr     = ',drips, dripw
        print *,'epotsn = ',epotsn,epot_sn_old
        print *,'frac   = ',fr_snow,fr_snow_old
        print *
        !call stop_model('flmlt(ibv) < 0.d0',255)
      endif
      flmlt = max( flmlt, 0.d0 )
 
ccc update fluxes that could have been changed by snow model
      if ( fm > 0.d0 ) then
        evap = epotsn/fm
        snsh = snshsn/fm
      endif

      return
      end subroutine snow_drv

      end module snow_drvm


