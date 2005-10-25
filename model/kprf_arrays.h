      integer, dimension (1:idm,1:jdm) ::
     & klstsv,        ! copy of k-index for momtum
     & jerlov         ! jerlov water type 1-5

      common /kpp_integ_arrays/
     .  klstsv,jerlov

      real, dimension (1:idm,1:jdm,2) :: 
     & t1sav,         ! upper sublayer temperature
     & s1sav,         ! upper sublayer salinity
     & tmlb,          ! temp in lyr. containing mlb.
     & smlb           ! saln in lyr. containing mlb

      real, dimension (1:idm,1:jdm) ::
     & hekman,        ! ekman layer thickness
     & hmonob,        ! monin-obukhov length
     & dpbl,          ! turbulent boundary layer depth
     & dpbbl,         ! bottom turbulent boundary layer depth
     & dpmold,        ! mixed layer depth at old time step
     & tmix,          ! mixed layer temperature
     & smix,          ! mixed layer salinity
     & thmix,         ! mixed layer potential density 
     & umix,  vmix    ! mixed layer velocity

      real, dimension (5) ::
     & betard,        ! red  extinction coefficient
     & betabl,        ! blue extinction coefficient
     & redfac         ! fract. of penetr. red light

      real, dimension (1:idm,1:jdm,4) ::
     & akpar          ! photosynthetically available radiation coefficent

c --- kpp variables
      real, dimension (1:idm,1:jdm,kdm+1) ::
     & zgrid          !  grid levels in meters
     &,vcty           !  vert. viscosity coefficient
     &,difs           !  vert. scalar diffusivity
     &,dift           !  vert. temperature diffusivity
     &,ghats          !  nonlocal transport

      real, dimension (1:idm,1:jdm) ::
     & buoflx,        ! surface buoyancy flux
     & bhtflx,        ! surface buoyancy flux from heat
     & mixflx,        ! mixed layer thermal energy flux
     & sswflx         ! short-wave rad'n flux at surface

      common /kpp_real_arrays/
     .  t1sav,s1sav,tmlb,smlb,hekman,hmonob,dpbl,dpbbl,
     .  dpmold,tmix,smix,thmix,umix,vmix,betard,betabl,redfac,akpar,
     .  zgrid,vcty,difs,dift,ghats,buoflx,bhtflx,mixflx,sswflx
c
      integer,parameter :: lookup=762
      integer
     & irimax(-lookup:lookup)
     &,nb

      common/giss_int2/
     &        irimax,nb
      save  /giss_int2/
c
      real
     & ribtbl(-lookup:lookup)
     &,ridb(  -lookup:lookup)
     &,slq2b( -lookup:lookup,-lookup:lookup)
     &,dri
     &,smb(   -lookup:lookup,-lookup:lookup)
     &,shb(   -lookup:lookup,-lookup:lookup)
     &,ssb(   -lookup:lookup,-lookup:lookup)
     &,back_ra_r(-39:117)
     &,sisamax(  -39:117)
     &,ra_rmax(  -39:117)
     &,c_y_r0(   -39:117)
     &,sm_r1(    -39:117)
     &,sh_r1(    -39:117)
     &,ss_r1(    -39:117)
     &,slq2_r1(  -39:117)
     &,b1,visc_cbu_limit,diff_cbt_limit
     &,theta_rcrp,theta_rcrn

      common/giss_real2/
     &        ribtbl,ridb,slq2b,dri,smb,shb,ssb,back_ra_r,
     &        sisamax,ra_rmax,c_y_r0,sm_r1,sh_r1,ss_r1,slq2_r1,
     &        b1,visc_cbu_limit,diff_cbt_limit,
     &        theta_rcrp,theta_rcrn
      save  /giss_real2/
