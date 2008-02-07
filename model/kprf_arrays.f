      module kprf_arrays
      use hycom_dim

      implicit none

      private

      public alloc_kprf_arrays

      public
     & klstsv,        ! copy of k-index for momtum
     & jerlov,         ! jerlov water type 1-5
     & t1sav,         ! upper sublayer temperature
     & s1sav,         ! upper sublayer salinity
     & tmlb,          ! temp in lyr. containing mlb.
     & smlb,           ! saln in lyr. containing mlb
     & hekman,        ! ekman layer thickness
     & hmonob,        ! monin-obukhov length
     & dpbl,          ! turbulent boundary layer depth
     & dpbbl,         ! bottom turbulent boundary layer depth
     & dpmold,        ! mixed layer depth at old time step
     & tmix,          ! mixed layer temperature
     & smix,          ! mixed layer salinity
     & thmix,         ! mixed layer potential density 
     & umix,  vmix,    ! mixed layer velocity
     & betard,        ! red  extinction coefficient
     & betabl,        ! blue extinction coefficient
     & redfac,         ! fract. of penetr. red light
     & akpar,          ! photosynthetically available radiation coefficent
     & zgrid,          !  grid levels in meters
     & vcty,           !  vert. viscosity coefficient
     & difs,           !  vert. scalar diffusivity
     & dift,           !  vert. temperature diffusivity
     & ghats,          !  nonlocal transport
     & buoflx,        ! surface buoyancy flux
     & bhtflx,        ! surface buoyancy flux from heat
     & mixflx,        ! mixed layer thermal energy flux
     & sswflx,         ! short-wave rad'n flux at surface
     & lookup,
     & irimax,
     & nb,
     & ribtbl,
     & ridb,
     & slq2b
     &,dri
     &,smb
     &,shb
     &,ssb
     &,back_ra_r
     &,sisamax
     &,ra_rmax
     &,c_y_r0
     &,sm_r1
     &,sh_r1
     &,ss_r1
     &,slq2_r1
     &,b1,visc_cbu_limit,diff_cbt_limit
     &,theta_rcrp,theta_rcrn


!! end public

      integer, allocatable, dimension (:,:) ::
     & klstsv,        ! copy of k-index for momtum
     & jerlov         ! jerlov water type 1-5

!!      common /kpp_integ_arrays/
!!     .  klstsv,jerlov

      real, allocatable, dimension (:,:,:) :: 
     & t1sav,         ! upper sublayer temperature
     & s1sav,         ! upper sublayer salinity
     & tmlb,          ! temp in lyr. containing mlb.
     & smlb           ! saln in lyr. containing mlb

      real, allocatable, dimension (:,:) ::
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

      real, allocatable, dimension (:,:,:) ::
     & akpar          ! photosynthetically available radiation coefficent

c --- kpp variables
      real, allocatable, dimension (:,:,:) ::
     & zgrid          !  grid levels in meters
     &,vcty           !  vert. viscosity coefficient
     &,difs           !  vert. scalar diffusivity
     &,dift           !  vert. temperature diffusivity
     &,ghats          !  nonlocal transport

      real, allocatable, dimension (:,:) ::
     & buoflx,        ! surface buoyancy flux
     & bhtflx,        ! surface buoyancy flux from heat
     & mixflx,        ! mixed layer thermal energy flux
     & sswflx         ! short-wave rad'n flux at surface

!!      common /kpp_real_arrays/
!!     .  t1sav,s1sav,tmlb,smlb,hekman,hmonob,dpbl,dpbbl,
!!     .  dpmold,tmix,smix,thmix,umix,vmix,betard,betabl,redfac,akpar,
!!     .  zgrid,vcty,difs,dift,ghats,buoflx,bhtflx,mixflx,sswflx
c
      integer,parameter :: lookup=762
      integer, save ::
     & irimax(-lookup:lookup)
     &,nb

!!      common/giss_int2/
!!     &        irimax,nb
!!      save  /giss_int2/
c
      real, save ::
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

!!      common/giss_real2/
!!     &        ribtbl,ridb,slq2b,dri,smb,shb,ssb,back_ra_r,
!!     &        sisamax,ra_rmax,c_y_r0,sm_r1,sh_r1,ss_r1,slq2_r1,
!!     &        b1,visc_cbu_limit,diff_cbt_limit,
!!     &        theta_rcrp,theta_rcrn
!!      save  /giss_real2/

      contains

!!! switch to this when doing parallelization
      subroutine alloc_kprf_arrays_local

      allocate(
     & klstsv(I_0H:I_1H,J_0H:J_1H),
     & jerlov(I_0H:I_1H,J_0H:J_1H) )

      allocate(
     & t1sav(I_0H:I_1H,J_0H:J_1H,2),
     & s1sav(I_0H:I_1H,J_0H:J_1H,2),
     & tmlb(I_0H:I_1H,J_0H:J_1H,2),
     & smlb(I_0H:I_1H,J_0H:J_1H,2) )

      allocate(
     & hekman(I_0H:I_1H,J_0H:J_1H),
     & hmonob(I_0H:I_1H,J_0H:J_1H),
     & dpbl(I_0H:I_1H,J_0H:J_1H),
     & dpbbl(I_0H:I_1H,J_0H:J_1H),
     & dpmold(I_0H:I_1H,J_0H:J_1H),
     & tmix(I_0H:I_1H,J_0H:J_1H),
     & smix(I_0H:I_1H,J_0H:J_1H),
     & thmix(I_0H:I_1H,J_0H:J_1H),
     & umix(I_0H:I_1H,J_0H:J_1H),  vmix(I_0H:I_1H,J_0H:J_1H) )

      allocate(
     & akpar(I_0H:I_1H,J_0H:J_1H,4) )

      allocate(
     & zgrid(I_0H:I_1H,J_0H:J_1H,kdm+1)
     &,vcty(I_0H:I_1H,J_0H:J_1H,kdm+1)
     &,difs(I_0H:I_1H,J_0H:J_1H,kdm+1)
     &,dift(I_0H:I_1H,J_0H:J_1H,kdm+1)
     &,ghats(I_0H:I_1H,J_0H:J_1H,kdm+1) )

      allocate(
     & buoflx(I_0H:I_1H,J_0H:J_1H),
     & bhtflx(I_0H:I_1H,J_0H:J_1H),
     & mixflx(I_0H:I_1H,J_0H:J_1H),
     & sswflx(I_0H:I_1H,J_0H:J_1H) )

      end subroutine alloc_kprf_arrays_local

      subroutine alloc_kprf_arrays

      allocate(
     & klstsv(IDM,JDM),
     & jerlov(IDM,JDM) )

      allocate(
     & t1sav(IDM,JDM,2),
     & s1sav(IDM,JDM,2),
     & tmlb(IDM,JDM,2),
     & smlb(IDM,JDM,2) )

      allocate(
     & hekman(IDM,JDM),
     & hmonob(IDM,JDM),
     & dpbl(IDM,JDM),
     & dpbbl(IDM,JDM),
     & dpmold(IDM,JDM),
     & tmix(IDM,JDM),
     & smix(IDM,JDM),
     & thmix(IDM,JDM),
     & umix(IDM,JDM),  vmix(IDM,JDM) )

      allocate(
     & akpar(IDM,JDM,4) )

      allocate(
     & zgrid(IDM,JDM,kdm+1)
     &,vcty(IDM,JDM,kdm+1)
     &,difs(IDM,JDM,kdm+1)
     &,dift(IDM,JDM,kdm+1)
     &,ghats(IDM,JDM,kdm+1) )

      allocate(
     & buoflx(IDM,JDM),
     & bhtflx(IDM,JDM),
     & mixflx(IDM,JDM),
     & sswflx(IDM,JDM) )

      end subroutine alloc_kprf_arrays



      end module kprf_arrays
