#include "rundeck_opts.h"

      MODULE TRACER_SOURCES

#ifdef TRACERS_ON

      USE TRACER_COM, only: ntm

      IMPLICIT NONE
      SAVE

!@param Laircr the number of layers of aircraft data read from file
!@param Lsulf the number of layers of sulfate SA data read from file
      INTEGER, PARAMETER :: 
     &                      Laircr      =25,
     &                      Lsulf       =23 ! not LM
#ifdef SHINDELL_STRAT_EXTRA
      REAL*8, PARAMETER ::  GLTic = 1.d-9 ! pppv
#endif   

!@param aircraft_Tyr1, aircraft_Tyr2 the starting and ending years
!@+     for transient tracer aircraft emissions (= means non transient)
      integer :: aircraft_Tyr1=0,aircraft_Tyr2=0

!@var airtracer 3D source of tracer from aircraft (on model levels)
      real*8, dimension(:,:,:), allocatable :: airtracer

#ifdef INTERACTIVE_WETLANDS_CH4
!@dbparam nn_or_zon approach to use for expanding wetlands 1=
!@+ zonal average, 0=nearest neighbor average
!@dbparam int_wet_dist to turn on/off interacive SPATIAL wetlands
!@dbparam ice_age if not = 0., allows no wetl emis for latitudes
!@+ poleward of +/- ice_age in degrees
!@dbparam ns_wet the number of source that is the wetlands src 
!@dbparam exclude_us_eu to exclude (=1) the U.S. and E.U. from 
!@+ interactive wetland distributiont
!@dbparam topo_lim upper limit on topography for int wetl dist
!@dbparam sat_lim lower limit on surf air temp for int wetl dist
!@dbparam gw_llim lower limit on ground wetness for int wetl dist
!@dbparam gw_ulim upper limit on ground wetness for int wetl dist
!@dbparam SW_lim lower limit on SW down flux for int wetl dist
!@param nra_ch4 number of running averages needed for int-wetlands
!@param nra_ncep # of ncep running averages needed for int-wetlands
!@param nday_ncep number of days in ncep running averages
!@param nday_ch4 number of days in running average
!@var maxHR_ch4 maximum number of sub-daily accumulations
!@var by_nday_ncep  1/real(nday_ncep)
!@var by_nday_ncep  1/real(nday_ch4)
!@var day_ncep daily NCEP temperature & precipitation 
!@var DRA_ch4 daily running average of model (prec,temp)
!@var avg_ncep running average (prec,temp) over nday_ncep days
!@var avg_model equivalent of avg_ncep, but based on model variables
!@var sum_ncep temp arrays for computing running average
!@var PRS_ch4 period running sum of model (prec,temp)
!@var HRA_ch4 hourly running average of model (prec,temp)
!@var iday_ncep current day (counter) of averaging period (prec,temp)
!@var iHch4 "hourly" index for averages of model (prec,temp)
!@var iDch4 "daily"  index for averages of model (prec,temp)
!@var i0ch4 ponter to current index in running sum of mode (prec,temp)
!@var first_ncep whether in the first ncep averaging period (prec,temp)
!@var first_mod whether in the first model averaging per.   (prec,temp)
!@var PTBA variable to hold the pressure, temperature, beta, and alpha
!@var PTBA1, PTBA2 for interpolations of PTBA
      integer, parameter :: nra_ch4 = 5, nra_ncep=2, n__prec=1,
     &  n__temp=2, n__SW=3, n__SAT=4, n__gwet=5, max_days=28, nncep=4
      integer :: int_wet_dist=0,exclude_us_eu=1,nn_or_zon=0,ns_wet=11
      real*8 :: topo_lim = 205.d0, sat_lim=-9.d0, 
     & gw_ulim=100.d0, gw_llim=18.d0, SW_lim=27.d0, ice_age=0.d0
      integer, parameter, dimension(nra_ch4) :: 
     &                                  nday_ch4=(/28,14,14,14,28/)
      integer, parameter, dimension(nra_ncep):: nday_ncep=(/28,14/)
      real*8, dimension(nra_ch4)             :: by_nday_ch4
      real*8, dimension(nra_ncep)            :: by_nday_ncep
      integer, dimension(nncep)      :: ncep_units,jmon_nc,jdlnc=0
      logical, dimension(nncep)      :: nc_first=.true.
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:):: day_ncep
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:):: DRA_ch4
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:)  :: avg_model,PRS_ch4
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:)  :: avg_ncep,sum_ncep
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:):: HRA_ch4
      integer, dimension(nra_ncep)           :: iday_ncep, i0_ncep,
     &                                          first_ncep
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: iHch4,iDch4,i0ch4,
     &                                          first_mod
      real*8, allocatable, dimension(:,:,:) ::  PTBA,PTBA1,PTBA2
      real*8, allocatable, dimension(:,:) :: add_wet_src
      integer :: maxHR_ch4
#endif
#endif
      END MODULE TRACER_SOURCES


      subroutine alloc_tracer_sources(grid)
!@SUM  To alllocate arrays whose sizes now need to be determined
!@+    at run-time
!@auth G.Faluvegi
      use domain_decomp_atm, only : dist_grid, getDomainBounds
      use domain_decomp_atm, only : write_parallel
      USE Dictionary_mod, only : get_param, is_set_param
      use resolution, only : lm
      use model_com, only : DTsrc
      use TimeConstants_mod, only: SECONDS_PER_HOUR
      use tracer_com, only : NTM
      use tracer_sources
      use fluxes, only : NIsurf
      IMPLICIT NONE

      type (dist_grid), intent(in) :: grid
      integer :: ier, J_1H, J_0H, I_1H, I_0H
      logical :: init = .false.
      real*8 :: DTsrc_LOCAL
      integer :: NIsurf_LOCAL

      if(init)return
      init=.true.
    
      call getDomainBounds( grid , J_STRT_HALO=J_0H, J_STOP_HALO=J_1H )
      I_0H = grid%I_STRT_HALO
      I_1H = grid%I_STOP_HALO

#ifdef INTERACTIVE_WETLANDS_CH4
      ! here I want to define how many surface calls expected in
      ! each day, but real DTsrc and NIsurf are not available yet
      ! so use local copy from the database:

      DTsrc_LOCAL = DTsrc
      if(is_set_param("DTsrc"))call get_param("DTsrc",DTsrc_LOCAL)
      NIsurf_LOCAL = NIsurf
      if(is_set_param("NIsurf"))call get_param("NIsurf",NIsurf_LOCAL)

      maxHR_ch4=24*NIsurf_LOCAL*NINT(SECONDS_PER_HOUR/DTsrc_LOCAL)
#endif
 
      allocate(airtracer(I_0H:I_1H,J_0H:J_1H,LM))
      airtracer = 0.

#ifdef INTERACTIVE_WETLANDS_CH4
      allocate( first_mod(I_0H:I_1H,J_0H:J_1H,nra_ch4) )
      allocate( iHch4(I_0H:I_1H,J_0H:J_1H,nra_ch4) )
      allocate( iDch4(I_0H:I_1H,J_0H:J_1H,nra_ch4) )
      allocate( i0ch4(I_0H:I_1H,J_0H:J_1H,nra_ch4) )
      allocate( day_ncep(I_0H:I_1H,J_0H:J_1H,max_days,nra_ncep) )
      allocate( DRA_ch4(I_0H:I_1H,J_0H:J_1H,max_days,nra_ch4) )
      allocate( avg_model(I_0H:I_1H,J_0H:J_1H,nra_ch4) )
      allocate( PRS_ch4(I_0H:I_1H,J_0H:J_1H,nra_ch4) )
      allocate( avg_ncep(I_0H:I_1H,J_0H:J_1H,nra_ncep) )
      allocate( sum_ncep(I_0H:I_1H,J_0H:J_1H,nra_ncep) )
      allocate( HRA_ch4(I_0H:I_1H,J_0H:J_1H,maxHR_ch4,nra_ch4) )
      allocate( PTBA(I_0H:I_1H,J_0H:J_1H,nncep) )
      allocate( PTBA1(I_0H:I_1H,J_0H:J_1H,nncep) )
      allocate( PTBA2(I_0H:I_1H,J_0H:J_1H,nncep) )
      allocate( add_wet_src(I_0H:I_1H,J_0H:J_1H) )
#endif      
      
      return
      end subroutine alloc_tracer_sources 
            
      
#ifdef SHINDELL_STRAT_EXTRA
      subroutine overwrite_GLT
!@sum L=1 overwriting of generic linear tracer    
!@vers 2013/03/26
!@auth Greg Faluvegi
C****
C**** Right now, there is just one L=1 source that changes 
C**** linearly in time (at 1% increase per year)
      USE RESOLUTION, only : im,jm
      USE MODEL_COM, only: itime,itimei,DTsrc
      USE GEOM, only: axyp,IMAXJ  
      USE ATM_COM, only: MA
      use OldTracer_mod, only: trname, vol2mass, itime_tr0
      USE TRACER_COM, only: trm,n_GLT
      USE TRACER_SOURCES, only: GLTic
      USE FLUXES, only : tr3Dsource
      USE DOMAIN_DECOMP_ATM, ONLY : getDomainBounds,grid,write_parallel
      
      IMPLICIT NONE
      
      INTEGER :: J_0, J_1, I_0, I_1
           
!@var by_s_in_yr recip. of # seconds in a year
!@var new_mr mixing ratio to overwrite in L=1 this time step
!@var new_mass mass to overwrite in L=1 this time step

      REAL*8 bydtsrc, by_s_in_yr, new_mr, new_mass
      integer i,j

      call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

      bydtsrc=1.d0/DTsrc
      !by_s_in_yr = 1.d0/(365.d0*24.d0*60.d0*60.d0)
      by_s_in_yr = 1.d0/SECONDS_PER_YEAR

C initial source is an overwriting of GLTic pppv, then add
C 1% every year, linearly in time. (note: vol2mass should
C just be 1 for this tracer, but kept it in here, in case
C we change that.)
      new_mr = GLTic * (1.d0 +
     &(Itime-ItimeI-itime_tr0(n_GLT))*DTsrc*by_s_in_yr*1.d-2) !pppv
      do j=J_0,J_1; do i=I_0,imaxj(j)
        new_mass = new_mr*vol2mass(n_GLT)*MA(1,i,j)*AXYP(i,j) ! kg
        tr3Dsource(i,j,1,1,n_GLT)=(new_mass-trm(i,j,1,n_GLT))*bydtsrc
        !i.e. tr3Dsource in kg/s 
      end do   ; end do

      return
      end subroutine overwrite_GLT
#endif


      SUBROUTINE get_aircraft_tracer(xyear,xday,phi,need_read)
!@sum  get_aircraft_tracer to define the 3D source of tracers from aircraft
!@auth Drew Shindell? / Greg Faluvegi / Jean Learner
!@ver  2.0 (based on DB396Tds3M23 -- adapted for AR5 emissions)
      USE RESOLUTION, only : im,jm
      USE RESOLUTION, only : lm
      use model_com, only: itime
      use domain_decomp_atm, only: GRID
      use domain_decomp_atm, only: getDomainBounds, write_parallel
      use constant, only: bygrav
      use filemanager, only: openunit,closeunit
      use fluxes, only: tr3Dsource
      use geom, only: axyp
      use OldTracer_mod, only: itime_tr0,trname
#ifdef TRACERS_SPECIAL_Shindell
      use TRACER_COM, only: n_NOx
#endif
#ifdef TRACERS_AEROSOLS_Koch
      use TRACER_COM, only: n_BCIA
#endif
#ifdef TRACERS_TOMAS
      use TRACER_COM, only: IDTECOB
#endif
!#ifdef TRACERS_AMP
!          use TRACER_COM, only: n_M_BC1_BC
!#endif
      use TRACER_COM, only:
     *                      nAircraft
      use tracer_sources, only: Laircr,aircraft_Tyr1,aircraft_Tyr2
     &     ,airtracer
      IMPLICIT NONE
 
      integer, intent(IN) :: xyear,xday
      real*8, dimension(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                  GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM),
     &     intent(IN) :: phi
      logical, intent(IN) :: need_read

      character(len=300) :: out_line
      integer, parameter :: nanns=0
#if (defined TRACERS_SPECIAL_Shindell) && (defined TRACERS_AEROSOLS_Koch)
!#if ((defined TRACERS_SPECIAL_Shindell) && (defined TRACERS_AEROSOLS_Koch)) ||\
!    ((defined TRACERS_SPECIAL_Shindell) && (defined TRACERS_AMP))
      integer, parameter :: nmons=2
#elif (defined TRACERS_SPECIAL_Shindell) && (defined TRACERS_TOMAS)
      integer, parameter :: nmons=2
#else
      integer, parameter :: nmons=1
#endif
      integer, dimension(nmons) :: mon_units, imon
      integer l,i,j,k,ll
      character*13, dimension(nmons) :: 
#if (defined TRACERS_SPECIAL_Shindell) && (defined TRACERS_AEROSOLS_Koch)
     *  mon_files=(/'NOx_AIRC ','BCIA_AIRC'/)
!#elif (defined TRACERS_SPECIAL_Shindell) && (defined TRACERS_AMP)
!     *  mon_files=(/'NOx_AIRC','M_BC1_BC_AIRC'/)
#elif (defined TRACERS_SPECIAL_Shindell) && (defined TRACERS_TOMAS)
     *  mon_files=(/'NOx_AIRC     ','AECOB_01_AIRC'/)
#elif (defined TRACERS_SPECIAL_Shindell)
     *  mon_files=(/'NOx_AIRC'/)
#elif (defined TRACERS_AEROSOLS_Koch)
     *  mon_files=(/'BCIA_AIRC'/)
#elif (defined TRACERS_TOMAS)
     *  mon_files=(/'AECOB_01_AIRC'/)
#else
     *  mon_files=(/'M_BC1_BC_AIRC'/)
#endif

      integer, dimension(nmons) :: mon_tracers ! define them later
#if (defined TRACERS_SPECIAL_Shindell) && (defined TRACERS_AEROSOLS_Koch)
!#if ((defined TRACERS_SPECIAL_Shindell) && (defined TRACERS_AEROSOLS_Koch)) ||\
!    ((defined TRACERS_SPECIAL_Shindell) && (defined TRACERS_AMP))
      logical, dimension(nmons) :: mon_bins=(/.true.,.true./) ! binary file?
#elif (defined TRACERS_SPECIAL_Shindell) && (defined TRACERS_TOMAS)
      logical, dimension(nmons) :: mon_bins=(/.true.,.true./) ! binary file?
#else /* this is for both TRACERS_AEROSOLS_Koch and TRACERS_AMP and TRACERS_TOMAS*/
      logical, dimension(nmons) :: mon_bins=(/.true./) ! binary file?
#endif
      real*8, dimension(GRID%I_STRT_HALO:GRID%I_STOP_HALO
     *     ,GRID%J_STRT_HALO:GRID%J_STOP_HALO,Laircr,nmons):: src
!@var zmod approx. geometric height at model layer(m), phi/grav
      real*8, dimension(LM)                :: zmod
!@var zairL heights of AR5 aircraft emissions (km)
      real*4, parameter, dimension(Laircr) :: zairL = ! alt in km:
     & (/0.305, 0.915, 1.525, 2.135, 2.745, 3.355, 3.965, 4.575, 5.185,
     & 5.795, 6.405, 7.015, 7.625, 8.235001, 8.845, 9.455001, 10.065,
     & 10.675, 11.285, 11.895, 12.505, 13.115, 13.725, 14.335, 14.945/)
      integer :: J_1, J_0, I_0, I_1
      logical :: trans_emis=.false.
      integer :: yr1=0, yr2=0
 
! Aircraft tracer source input is monthly, on 25 levels.
! Read it in here and interpolated each day.

#if (defined TRACERS_SPECIAL_Shindell) && (defined TRACERS_AEROSOLS_Koch)
      mon_tracers(1)=n_NOx
      mon_tracers(2)=n_BCIA
#elif (defined TRACERS_SPECIAL_Shindell) && (defined TRACERS_TOMAS)
      mon_tracers(1)=n_NOx
      mon_tracers(2)=IDTECOB
!#elif (defined TRACERS_SPECIAL_Shindell) && (defined TRACERS_AMP)
!      mon_tracers(1)=n_NOx
!      mon_tracers(2)=n_M_BC1_BC
#elif (defined TRACERS_SPECIAL_Shindell)
      mon_tracers(1)=n_NOx
#elif (defined TRACERS_AEROSOLS_Koch)
      mon_tracers(1)=n_BCIA
#elif (defined TRACERS_TOMAS)
      mon_tracers(1)=IDTECOB
#else
      mon_tracers(1)=n_M_BC1_BC
#endif
      do k=1,nmons
        if (mon_tracers(k) == 0) then
          call stop_model("mon_tracers(k) not defined",255)
        endif
        if (itime < itime_tr0(mon_tracers(k))) cycle
        call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1)
        call getDomainBounds(grid, I_STRT=I_0, I_STOP=I_1)

! Monthly sources are interpolated to the current day
! Units are KG(N)/m2/s, so no conversion is necessary:
        if(aircraft_Tyr1==aircraft_Tyr2)then
          trans_emis=.false.; yr1=0; yr2=0
        else
          trans_emis=.true.; yr1=aircraft_Tyr1; yr2=aircraft_Tyr2
        endif
        if (trans_emis .and. xyear < 1900) return !<-- hardcode for NO AIRCRAFT before 1900

        if(need_read) then

        call openunit(mon_files(k),mon_units(k),mon_bins(k))
        call read_monthly_3Dsources(Laircr,mon_units(k),
     &   src(:,:,:,k),trans_emis,yr1,yr2,xyear,xday)
        call closeunit(mon_units(k))

! Place aircraft sources onto model levels:
        airtracer = 0.
        do j=J_0,J_1
          do i=I_0,I_1
            zmod(:)=phi(i,j,:)*bygrav*1.d-3 ! km
            do LL=1,Laircr
              if(src(i,j,LL,k) > 0.)then
                loop_l: do L=1,LM
                  if(zairL(LL) <= zmod(L)) then
      airtracer(i,j,l) = airtracer(i,j,l) + src(i,j,LL,k)*axyp(i,j)
                    exit loop_l
                  endif
                if(L==LM)call stop_model("aircraft level problem",255)
                enddo loop_l
              endif  ! is there a source?
            enddo   ! LL
          enddo    ! I
        enddo     ! J

        endif                     ! need_read?

        tr3Dsource(I_0:I_1,J_0:J_1,:,nAircraft,mon_tracers(k)) =
     &    airtracer(I_0:I_1,J_0:J_1,:)
      enddo ! k

      return
      end subroutine get_aircraft_tracer
 

      subroutine check_aircraft_sectors(n_NOx)
!@sum check_aircraft_sectors checks parameters for user-
!@+ set sector for NOx aircraft source.
!@auth Greg Faluvegi
      use TracerSource_mod, only: TracerSource3D
      use Tracer_mod, only: Tracer
      use tracer_com, only: nAircraft, tracers,
     & sect_name,num_sectors,
     & n_max_sect,ef_fact,num_regions,ef_fact,ef_fact3d
      use Dictionary_mod, only: sync_param
      IMPLICIT NONE
      integer, intent(in) :: n_NOx
      integer :: i,j,ns,nsect,nn
      character*124 :: tr_sectors_are
      character*32 :: pname

      type (TracerSource3D), pointer :: source
      class (Tracer), pointer :: pTracer

      tr_sectors_are = ' '
      pTracer => tracers%getReference('NOx')
      source => pTracer%sources3D(nAircraft)

      pname='NOx_AIRC_sect'
      call sync_param(pname,tr_sectors_are)
      source%num_tr_sectors = 0

      i=1
      do while(i < len(tr_sectors_are))
        j=index(tr_sectors_are(i:len(tr_sectors_are))," ")
        if (j > 1) then
          source%num_tr_sectors = source%num_tr_sectors + 1
          i=i+j
        else
          i=i+1
        end if
      enddo
      ns=source%num_tr_sectors
      if(ns > n_max_sect)
     &call stop_model("num_tr_sectors3D problem",255)
      if(ns > 0)then
        read(tr_sectors_are,*) source%tr_sect_name(1:ns)
        do nsect=1,ns
          source%tr_sect_index(nsect) = 0
          loop_nn: do nn=1,num_sectors
            if(trim(source%tr_sect_name(nsect)) ==
     &         trim(sect_name(nn))) then
              source%tr_sect_index(nsect) = nn
              ef_fact3d(nn,1:num_regions)=
     &        ef_fact(nn,1:num_regions)
              exit loop_nn
            endif
          enddo loop_nn
        enddo
      endif

      return
      end subroutine check_aircraft_sectors
 

      SUBROUTINE read_aero(field,fn)
!@sum read_aero, read in monthly mean fields of SO2 concentration or
!@+ DMS concentration or sulfate surface area for use in chemistry.
!@+ This source is not applied directly to the tracer mass like
!@+ other 3D sources are but used for e.g. HOx, N2O5 chemistry.
!@+ call it like so:
!@+    call read_aero(dms_offline,'DMS_FIELD')
!@+    call read_aero(sulfate,'SULFATE_SA')
!@+    call read_aero(so2_offline,'SO2_FIELD')
!@auth Drew Shindell / Greg Faluvegi
      USE RESOLUTION, only : ptop,psf
      USE RESOLUTION, only : lm
      USE RESOLUTION, only : im,jm
      use model_com, only: modelEclock
      USE DYNAMICS, only : sig
      USE DOMAIN_DECOMP_ATM, only: GRID
      USE DOMAIN_DECOMP_ATM, only: getDomainBounds, write_parallel
      USE FILEMANAGER, only: openunit,closeunit
      use TRACER_SOURCES, only: Lsulf
 
      IMPLICIT NONE
      real*8, dimension(GRID%I_STRT_HALO:GRID%I_STOP_HALO
     *     ,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM) ::field
      CHARACTER(LEN=*),  INTENT(IN) :: fn

!@var nmons: number of monthly input files
!@param Psulf pressure levels of the input file
      real*8, parameter, dimension(Lsulf) :: Psulf = (/
     & 0.9720D+03,0.9445D+03,0.9065D+03,
     & 0.8515D+03,0.7645D+03,0.6400D+03,0.4975D+03,0.3695D+03,
     & 0.2795D+03,0.2185D+03,0.1710D+03,0.1335D+03,0.1016D+03,
     & 0.7120D+02,0.4390D+02,0.2470D+02,0.1390D+02,0.7315D+01,
     & 0.3045D+01,0.9605D+00,0.3030D+00,0.8810D-01,0.1663D-01/)
      integer, parameter :: ncalls=3
      integer, dimension(ncalls):: mon_units
      integer i,j,iu,k,l,nc
      character*80 :: title
      character(len=300) :: out_line
      logical, dimension(ncalls) :: mon_bins=(/.true.,.true.,.true./)
      REAL*8, DIMENSION(LM)    :: pres,srcLout
      REAL*8, DIMENSION(Lsulf) :: srcLin
      REAL*8, DIMENSION(GRID%I_STRT_HALO:GRID%I_STOP_HALO
     *     ,GRID%J_STRT_HALO:GRID%J_STOP_HALO,Lsulf,ncalls):: src
      logical :: trans_emis=.false.
      INTEGER :: J_1, J_0, J_0H, J_1H, I_0, I_1

      call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1)
      call getDomainBounds(grid, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

      select case(trim(fn))
      case('DMS_FIELD') ; nc=1
      case('SO2_FIELD') ; nc=2
      case('SULFATE_SA'); nc=3
      case default
        call stop_model('please address filename in read_aero',255)
      end select
 
C**** Input file is monthly, on LM levels.
C     Read it in here and interpolated each day.
C     Interpolation in the vertical.
  
      call openunit(trim(fn),mon_units(nc),mon_bins(nc))
      call read_monthly_3Dsources(Lsulf,mon_units(nc),
     &   src(:,:,:,nc),trans_emis,0,0,modelEclock%year(),
     &   modelEclock%dayOfYear())
      call closeunit(mon_units(nc))
C====
C====   Place field onto model levels
C====              
      PRES(:)=SIG(:)*(PSF-PTOP)+PTOP
      DO J=J_0,J_1; DO I=I_0,I_1
        srcLin(1:Lsulf)=src(I,J,1:Lsulf,nc)
        call LOGPINT(Lsulf,Psulf,srcLin,LM,PRES,srcLout,.true.)
        field(I,J,1:LM)=srcLout(1:LM)
      END DO   ; END DO    
  
      return
      END SUBROUTINE read_aero
 

      SUBROUTINE read_monthly_3Dsources
     & (Ldim,iu,data1,trans_emis,yr1,yr2,xyear,xday)
!@sum Read in monthly sources and interpolate to current day
!@auth Jean Lerner and others / Greg Faluvegi
      USE RESOLUTION, only : im,jm
      USE MODEL_COM, only: idofm=>JDmidOfM
      USE FILEMANAGER, only : NAMEUNIT
      USE DOMAIN_DECOMP_ATM, only : GRID,getDomainBounds,READT_PARALLEL
     &     ,REWIND_PARALLEL
     &     ,write_parallel,backspace_parallel,am_i_root
      implicit none
!@var Ldim how many vertical levels in the read-in file?
!@var L dummy vertical loop variable
      integer :: Ldim,L,imon,iu,ipos,k,nn,k2,kstep=10
      character(len=300) :: out_line
      real*8 :: frac, alpha
      real*8, DIMENSION(GRID%I_STRT_HALO:GRID%I_STOP_HALO
     *     ,GRID%J_STRT_HALO:GRID%J_STOP_HALO) ::A2D,B2D,dummy
      real*8, DIMENSION(GRID%I_STRT_HALO:GRID%I_STOP_HALO
     *     ,GRID%J_STRT_HALO:GRID%J_STOP_HALO,Ldim) ::tlca,tlcb,data1
     *     ,sfc_a,sfc_b
      logical, intent(in):: trans_emis
      integer, intent(in):: yr1,yr2,xyear,xday
     
      integer :: J_0, J_1, I_0, I_1

      call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1)     
      call getDomainBounds(grid, I_STRT=I_0, I_STOP=I_1)     

C No doubt this code can be combined/compressed, but I am going to
C do the transient and non-transient cases separately for the moment:

! -------------- non-transient emissions ----------------------------!
      if(.not.trans_emis) then
C
      imon=1                ! imon=January
      if (xday <= 16)  then ! DAY in Jan 1-15, first month is Dec
        if(am_i_root())write(6,*) 'Not using this first record:'
        call readt_parallel(grid,iu,nameunit(iu),dummy,Ldim*11)
        do L=1,Ldim
          call readt_parallel(grid,iu,nameunit(iu),A2D,1)
          tlca(I_0:I_1,J_0:J_1,L)=A2D(I_0:I_1,J_0:J_1)
        enddo  
        call rewind_parallel(iu)
      else              ! DAY is in Jan 16 to Dec 16, get first month
        do while(xday > idofm(imon) .AND. imon <= 12)
          imon=imon+1
        enddo
        if(imon/=2)then ! avoids advancing records at start of file
          if(am_i_root())write(6,*) 'Not using this first record:'
          call readt_parallel(grid,iu,nameunit(iu),dummy,Ldim*(imon-2))
        end if
        do L=1,Ldim
          call readt_parallel(grid,iu,nameunit(iu),A2D,1)
          tlca(I_0:I_1,J_0:J_1,L)=A2D(I_0:I_1,J_0:J_1)
        enddo   
        if(imon==13) call rewind_parallel(iu)
      end if
      do L=1,Ldim
        call readt_parallel(grid,iu,nameunit(iu),B2D,1)
        tlcb(I_0:I_1,J_0:J_1,L)=B2D(I_0:I_1,J_0:J_1)
      enddo 
c**** Interpolate two months of data to current day
      frac = float(idofm(imon)-xday)/(idofm(imon)-idofm(imon-1))
      data1(I_0:I_1,J_0:J_1,:) =
     & tlca(I_0:I_1,J_0:J_1,:)*frac + tlcb(I_0:I_1,J_0:J_1,:)*(1.-frac)
      write(out_line,*) '3D source monthly factor=',frac
      call write_parallel(trim(out_line))

! --------------- transient emissions -------------------------------!
      else
        ! 3D source files as of now have no meta-data so assume
        ! that transient time slices are decadal:
        kstep=10
        ipos=1
        k2=yr1
        alpha=0.d0 ! before start year, use start year value
        if(xyear>yr2.or.(xyear==yr2.and.xday>=183))then
          alpha=1.d0 ! after end year, use end year value
          ipos=(yr2-yr1)/kstep
          k2=yr2-kstep
        endif
        do k=yr1,yr2-kstep,kstep
          if(xyear>k .or. (xyear==k.and.xday>=183)) then
            if(xyear<k+kstep .or. (xyear==k+kstep.and.xday<183))then
              ipos=1+(k-yr1)/kstep ! (integer artithmatic)
              alpha=real(xyear-k)/real(kstep)
              k2=k
              exit
            endif
          endif
        enddo
!
! read the two necessary months from the first decade:
!
      imon=1                ! imon=January
      if (xday <= 16)  then ! DAY in Jan 1-15, first month is Dec
        if(am_i_root())write(6,*) 'Not using this first record:'
        call readt_parallel
     &  (grid,iu,nameunit(iu),dummy,(ipos-1)*12*Ldim+Ldim*11)
        do L=1,Ldim
          call readt_parallel(grid,iu,nameunit(iu),A2D,1)
          tlca(I_0:I_1,J_0:J_1,L)=A2D(I_0:I_1,J_0:J_1)
        enddo
        do nn=1,12*Ldim; call backspace_parallel(iu); enddo
      else              ! DAY is in Jan 16 to Dec 16, get first month
        do while(xday > idofm(imon) .AND. imon <= 12)
          imon=imon+1
        enddo
        if(imon/=2 .or. ipos/=1)then ! avoids advancing records at start of file
          if(am_i_root())write(6,*) 'Not using this first record:' 
          call readt_parallel
     &    (grid,iu,nameunit(iu),dummy,(ipos-1)*12*Ldim+Ldim*(imon-2))
        end if
        do L=1,Ldim
          call readt_parallel(grid,iu,nameunit(iu),A2D,1)
          tlca(I_0:I_1,J_0:J_1,L)=A2D(I_0:I_1,J_0:J_1)
        enddo
        if(imon==13)then
          do nn=1,12*Ldim; call backspace_parallel(iu); enddo
        endif
      end if
CCCCC write(6,*) 'Not using this first record:'
CCCCC call readt_parallel(grid,iu,nameunit(iu),dummy,Ldim*(imon-1))
      do L=1,Ldim
        call readt_parallel(grid,iu,nameunit(iu),B2D,1)
        tlcb(I_0:I_1,J_0:J_1,L)=B2D(I_0:I_1,J_0:J_1)
      enddo
      frac = float(idofm(imon)-xday)/(idofm(imon)-idofm(imon-1))
      sfc_a(I_0:I_1,J_0:J_1,:) =
     & tlca(I_0:I_1,J_0:J_1,:)*frac + tlcb(I_0:I_1,J_0:J_1,:)*(1.-frac)
      call rewind_parallel( iu )

      ipos=ipos+1
      imon=1                ! imon=January
      if (xday <= 16)  then ! DAY in Jan 1-15, first month is Dec
        if(am_i_root())write(6,*) 'Not using this first record:'
        call readt_parallel
     &  (grid,iu,nameunit(iu),dummy,(ipos-1)*12*Ldim+Ldim*11)
        do L=1,Ldim
          call readt_parallel(grid,iu,nameunit(iu),A2D,1)
          tlca(I_0:I_1,J_0:J_1,L)=A2D(I_0:I_1,J_0:J_1)
        enddo
        do nn=1,12*Ldim; call backspace_parallel(iu); enddo
      else              ! DAY is in Jan 16 to Dec 16, get first month
        do while(xday > idofm(imon) .AND. imon <= 12)
          imon=imon+1
        enddo
        if(am_i_root())write(6,*) 'Not using this first record:'
        call readt_parallel
     &  (grid,iu,nameunit(iu),dummy,(ipos-1)*12*Ldim+Ldim*(imon-2))
        do L=1,Ldim
          call readt_parallel(grid,iu,nameunit(iu),A2D,1)
          tlca(I_0:I_1,J_0:J_1,L)=A2D(I_0:I_1,J_0:J_1)
        enddo
        if(imon==13)then
          do nn=1,12*Ldim; call backspace_parallel(iu); enddo
        endif
      end if
CCCCCCwrite(6,*) 'Not using this first record:'
CCCCCCcall readt_parallel(grid,iu,nameunit(iu),dummy,Ldim*(imon-1))
      do L=1,Ldim
        call readt_parallel(grid,iu,nameunit(iu),B2D,1)
        tlcb(I_0:I_1,J_0:J_1,L)=B2D(I_0:I_1,J_0:J_1)
      enddo
      frac = float(idofm(imon)-xday)/(idofm(imon)-idofm(imon-1))
      sfc_b(I_0:I_1,J_0:J_1,:) =
     & tlca(I_0:I_1,J_0:J_1,:)*frac + tlcb(I_0:I_1,J_0:J_1,:)*(1.-frac)

! now interpolate between the two time periods:

      data1(I_0:I_1,J_0:J_1,:) = sfc_a(I_0:I_1,J_0:J_1,:)*(1.d0-alpha) 
     & + sfc_b(I_0:I_1,J_0:J_1,:)*alpha

      write(out_line,*) '3D source at',
     &100.d0*alpha,' % of period this day ',k2,' to this day ',k2+kstep,
     &' and monthly fraction= ',frac 
      call write_parallel(trim(out_line))

      endif ! transient or not

      return
      end subroutine read_monthly_3Dsources


      subroutine get_CH4_IC(icall)
!@sum get_CH4_IC to generate initial conditions for methane.
!@vers 2013/03/26
!@auth Greg Faluvegi/Drew Shindell
      USE RESOLUTION, only : ls1
      USE RESOLUTION, only : im,jm,lm
      USE MODEL_COM, only  : DTsrc
      USE DOMAIN_DECOMP_ATM, only : GRID, getDomainBounds,
     *     write_parallel,am_i_root
      USE GEOM, only       : axyp,lat2d_dg
      USE ATM_COM, only: MA
      USE CONSTANT, only: mair
      use OldTracer_mod, only: vol2mass
      USE TRACER_COM, only : trm, n_CH4, nOverwrite
      USE FLUXES, only: tr3Dsource
      USE TRCHEM_Shindell_COM, only: CH4altT, CH4altX, ch4_init_sh,
     *     ch4_init_nh,fix_CH4_chemistry
 
      IMPLICIT NONE
      integer, intent(in) :: icall
 
!@var CH4INIT temp variable for ch4 initial conditions
!@var I,J,L dummy loop variables
!@param bymair 1/molecular wt. of air = 1/mair
!@var icall =1 (during run) =0 (first time)
      REAL*8, PARAMETER :: bymair = 1.d0/mair
      REAL*8 CH4INIT,bydtsrc
      INTEGER I, J, L
      integer :: J_0, J_1, I_0, I_1

      call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1)
      call getDomainBounds(grid, I_STRT=I_0, I_STOP=I_1)

      bydtsrc=1.d0/DTsrc
C     First, the troposphere:
      DO J=J_0,J_1
      DO I=I_0,I_1
C       Initial latitudinal gradient for CH4:
        IF(LAT2D_DG(I,J) < 0.) THEN ! Southern Hemisphere
          CH4INIT=ch4_init_sh*vol2mass(n_CH4)*1.d-6
        ELSE                        ! Northern Hemisphere
          CH4INIT=ch4_init_nh*vol2mass(n_CH4)*1.d-6
        ENDIF
        select case(icall)
        case(0) ! initial conditions
          DO L=1,LS1-1
            trm(i,j,l,n_CH4) = MA(L,I,J)*CH4INIT*AXYP(I,J)
          END DO
        case(1) ! overwriting
          DO L=1,LS1-1
            tr3Dsource(i,j,l,nOverwrite,n_CH4) = (MA(L,I,J)*
     &           CH4INIT*AXYP(I,J)-trm(i,j,l,n_CH4))*bydtsrc
          END DO
        end select
      END DO
      END DO
 
C     Now, the stratosphere:
      do L=LS1,LM
      do j=J_0,J_1
      do i=I_0,I_1
 
c     Define stratospheric ch4 based on HALOE obs for tropics
c     and extratropics and scale by the ratio of initial troposphere
c     mixing ratios to 1.79 (observed):
        IF(LAT2D_DG(I,J) < 0.) THEN ! Southern Hemisphere
          CH4INIT=ch4_init_sh/1.79d0*vol2mass(n_CH4)*1.E-6
        ELSE                        ! Northern Hemisphere
          CH4INIT=ch4_init_nh/1.79d0*vol2mass(n_CH4)*1.E-6
        ENDIF
        IF(ABS(LAT2D_DG(I,J)) > 30.) THEN ! extratropics
          CH4INIT=CH4INIT*CH4altX(L)
        ELSE                              ! tropics
          CH4INIT=CH4INIT*CH4altT(L)
        END IF
        select case(icall)
        case(0) ! initial conditions
          trm(i,j,l,n_CH4) = MA(L,I,J)*CH4INIT*AXYP(I,J)
        case(1) ! overwriting
          tr3Dsource(i,j,l,nOverwrite,n_CH4) = (MA(L,I,J)*
     &    CH4INIT*AXYP(I,J)-trm(i,j,l,n_CH4))*bydtsrc
        end select
      end do   ! i
      end do   ! j
      end do   ! l
 
      RETURN
      end subroutine get_CH4_IC
   
      
      SUBROUTINE get_sza(I,J,tempsza)
!@sum get_sza calculates the solar angle.  The intention is
!@+   that this routine will only be used when the COSZ1 from the 
!@+   radiation code is < 0, i.e. SZA > 90 deg.
!@auth Greg Faluvegi, based on the Harvard CTM routine SCALERAD
      USE RESOLUTION, only : im,jm
      USE MODEL_COM, only: nday,Itime,DTsrc 
      use TimeConstants_mod, only: SECONDS_PER_HOUR, SECONDS_PER_DAY,
     &                             INT_DAYS_PER_YEAR
      USE CONSTANT, only: PI, radian
      USE TRCHEM_Shindell_COM, only: byradian
      use geom, only : lon2d_dg,lat2d_dg
      use domain_decomp_atm, only : grid,hasSouthPole, hasNorthPole
  
      IMPLICIT NONE
      REAL*8, INTENT(OUT) :: tempsza
      INTEGER, INTENT(IN) :: I,J

!@param ANG1 ?
!@var DX degree width of a model grid cell (360/IM)
!@var DY degree width of a model grid cell ~(180/JM)
!@var tempsza the solar zenith angle, returned in degrees
!@var P1,P2,P3 ? angles needed to compute COS(SZA) in degrees    
!@var VLAT,VLOM latitude and longitude in degrees
!@var current julian time in seconds
!@var FACT,temp are temp variables
      REAL*8, PARAMETER ::  ANG1 = 90.d0/91.3125d0
      REAL*8 P1,P2,P3,VLAT,VLON,TIMEC,FACT,temp
      INTEGER DY

      DY   = NINT(180./REAL(JM)) ! only used for latlon grid
      vlat = lat2d_dg(i,j)
      if (J == 1  .and. hasSouthPole(grid))
     &     vlat= -90. + 0.5*REAL(DY)
      if (j == JM .and. hasNorthPole(grid))
     &     vlat=  90. - 0.5*REAL(DY)
      vlon = -lon2d_dg(i,j) ! i=1 in lon2d_dg is -180 + dlon/2
      TIMEC = (mod(itime,INT_DAYS_PER_YEAR*nday) + nday + 0.5)*DTsrc  
      P1 = 15.*(TIMEC/SECONDS_PER_HOUR - VLON/15. - 12.)
      FACT = (TIMEC/SECONDS_PER_DAY - 81.1875)*ANG1
      P2 = 23.5*SIN(FACT*radian)
      P3 = VLAT
      temp = (SIN(P3*radian)*SIN(P2*radian)) +
     &(COS(P1*radian)*COS(P2*radian)*COS(P3*radian))
      tempsza = acos(temp)*byradian

      RETURN
      END SUBROUTINE get_sza
 
      subroutine interpolateAltitude()
      USE TRCHEM_Shindell_COM,only:LCOalt,PCOalt,
     &     ClOXaltIN, ClOxalt,
     &     BrOXaltIN, BrOxalt,
     &     HClAltIN, Hclalt,
     &     ClONO2altIN, ClONO2alt
      USE ATM_COM, only: pmidl00
      USE RESOLUTION, only : LM

C     Interpolate ClOx altitude-dependence to model resolution:
      CALL LOGPINT(LCOalt,PCOalt,ClOxaltIN,LM,PMIDL00,ClOxalt,.true.)
C     Interpolate BrOx altitude-dependence to model resolution:
      CALL LOGPINT(LCOalt,PCOalt,BrOxaltIN,LM,PMIDL00,BrOxalt,.true.)
C     Interpolate HCl altitude-dependence to model resolution:
      CALL LOGPINT(LCOalt,PCOalt,HClaltIN,LM,PMIDL00,HClalt,.true.)
C     Interpolate ClONO2 altitude-dependence to model resolution:
      CALL
     &    LOGPINT(LCOalt,PCOalt,ClONO2altIN,LM,PMIDL00,ClONO2alt,.true.)
      end subroutine interpolateAltitude
 
      SUBROUTINE LOGPINT(LIN,PIN,AIN,LOUT,POUT,AOUT,min_zero)
!@sum LOGPINT does vertical interpolation of column variable,
!@+   linearly in ln(P).
!@auth Greg Faluvegi
 
      IMPLICIT NONE

!@var LIN number of levels for initial input variable
!@var LOUT number of levels for output variable
!@var PIN pressures at LIN levels
!@var POUT pressures at LOUT levels
!@var AIN initial input column variable
!@var AOUT output (interpolated) column variable
!@var LNPIN natural log of PIN
!@var LNPOUT natural log of POUT
!@var min_zero if true, don't allow negatives
!@var slope slope of line used for extrapolations
 
      LOGICAL, INTENT(IN)                 :: min_zero
      INTEGER, INTENT(IN)                 :: LIN, LOUT
      REAL*8, INTENT(IN), DIMENSION(LIN)  :: PIN, AIN
      REAL*8, INTENT(IN),DIMENSION(LOUT)  :: POUT
      REAL*8, INTENT(OUT),DIMENSION(LOUT) :: AOUT 
      REAL*8, DIMENSION(LIN)              :: LNPIN
      REAL*8, DIMENSION(LOUT)             :: LNPOUT       
      INTEGER L1,L2
      REAL*8 slope
C      
      LNPIN(:) = LOG(PIN(:))                   ! take natural log
      LNPOUT(:)= LOG(POUT(:))                  ! of pressures
C      
      DO L1=1,LOUT
       IF (LNPOUT(L1)>LNPIN(1)) THEN        ! extrapolate
         slope=(AIN(2)-AIN(1))/(LNPIN(2)-LNPIN(1))
         AOUT(L1)=AIN(1)-slope*(LNPIN(1)-LNPOUT(L1))
       ELSE IF (LNPOUT(L1) < LNPIN(LIN)) THEN ! extrapolate
         slope=(AIN(LIN)-AIN(LIN-1))/(LNPIN(LIN)-LNPIN(LIN-1))
         AOUT(L1)=AIN(LIN)+slope*(LNPOUT(L1)-LNPIN(LIN))
       ELSE                                    ! interpolate
        DO L2=1,LIN-1
         IF(LNPOUT(L1) == LNPIN(L2)) THEN
           AOUT(L1)=AIN(L2)
         ELSE IF(LNPOUT(L1) == LNPIN(L2+1)) THEN
           AOUT(L1)=AIN(L2+1)
         ELSE IF(LNPOUT(L1) < LNPIN(L2) .and.
     &   LNPOUT(L1) > LNPIN(L2+1)) THEN
           AOUT(L1)=(AIN(L2)*(LNPIN(L2+1)-LNPOUT(L1)) 
     &     +AIN(L2+1)*(LNPOUT(L1)-LNPIN(L2)))/(LNPIN(L2+1)-LNPIN(L2))
         END IF
        END DO
       END IF
      END DO          
C      
C If necessary: limit interpolated array to positive numbers:
C
      IF(min_zero)AOUT(:)=MAX(0.d0,AOUT(:))
C           
      RETURN
      END SUBROUTINE LOGPINT
 
 
#ifdef TRACERS_ON
      SUBROUTINE special_layers_init
!@sum special_layers_init determined special altitude levels that
!@+ are important to Drew Shindell's tracer code.  This is done to
!@+ so that hardcoded levels are not needed when switching
!@+ vertical resolutions.
!@auth Greg Faluvegi
      USE RESOLUTION, only : ptop,psf
      USE RESOLUTION, only : LM
      USE DYNAMICS, only : sig
      USE TRCHEM_Shindell_COM, only: L75P,L75M,F75P,F75M, ! FACT1
     &                           L569P,L569M,F569P,F569M  ! CH4
 
      IMPLICIT NONE

!@var natural log of nominal pressure for verticle interpolations
!@var log75 natural log of 75 hPa
!@var log569 natural log of 569 hPa
      REAL*8 log75,log569
      REAL*8, DIMENSION(LM) :: LOGP
      INTEGER L
       
      LOGP(:)=LOG(SIG(:)*(PSF-PTOP)+PTOP)
      log75=LOG(75.d0)
      log569=LOG(569.d0)
  
      DO L=1,LM-1
        IF(LOGP(L) > log75 .and. LOGP(L+1) < log75) THEN
          L75P=L+1 ! these are for FACT1 variable in strat overwrite
          L75M=L   ! hence effects several tracers
          F75P=(log75-LOGP(L75M))/(LOGP(L75P)-LOGP(L75M))
          F75M=(LOGP(L75P)-log75)/(LOGP(L75P)-LOGP(L75M))
        END IF
        IF(LOGP(L) > log569 .and. LOGP(L+1) < log569) THEN
          L569P=L+1 ! these are for the CH4 strat overwrite
          L569M=L
          F569P=(log569-LOGP(L569M))/(LOGP(L569P)-LOGP(L569M))
          F569M=(LOGP(L569P)-log569)/(LOGP(L569P)-LOGP(L569M))
        END IF
      END DO
       
      RETURN
      END SUBROUTINE special_layers_init
#endif


#ifdef INTERACTIVE_WETLANDS_CH4
      subroutine running_average(var_in,I,J,nicall,m)
!@sum running_average keeps a running average of the model variables
!@+ for use with the interactive wetlands CH4. Currently: 1st layer
!@+ ground temperature, precipitaion, downward SW rad flux, surface
!@+ air temperature, and ground wetness. 
!@+ I suppose I could generalized this in the future.
!@auth Greg Faluvegi
C
C**** Global variables:
c
      USE MODEL_COM, only: DTsrc
      use TimeConstants_mod, only: SECONDS_PER_HOUR, HOURS_PER_DAY
      USE TRACER_SOURCES, only: iH=>iHch4,iD=>iDch4,i0=>i0ch4,
     & first_mod,HRA=>HRA_ch4,DRA=>DRA_ch4,PRS=>PRS_ch4,
     & nday_ch4,by_nday_ch4,nra_ch4,maxHR_ch4,avg_model
C
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments
C     
!@var var_in model variable which will be used in the average
!@var nicall number of times routine is called per main timesetep
!@var m the index of the variable in question
!@var temp just for holding current day average for use in avg_model
!@var nmax number of accumulations in one day
!@var bynmax reciprocal of nmax
      real*8, intent(IN) :: var_in, nicall
      real*8 temp, bynmax
      integer, intent(IN):: m, I, J
      integer n, nmax
      
      if(m > nra_ch4.or.m < 1)call stop_model('nra_ch4 problem',255)
      if(iH(I,J,m) < 0.or.iH(I,J,m) > maxHR_ch4) then
        write(6,*) 'IJM iH maxHR_ch4=',I,J,m,iH(I,J,m),maxHR_ch4
        call stop_model('iH or maxHR_ch4 problem',255)
      endif
      nmax=NINT(HOURS_PER_DAY*nicall*SECONDS_PER_HOUR/DTsrc)
      bynmax=1.d0/real(nmax)
      by_nday_ch4(:)=1.d0/real(nday_ch4(:))
      iH(I,J,m) = iH(I,J,m) + 1
      HRA(I,J,iH(I,J,m),m) = var_in
      ! do no more, unless it is the end of the day:
      
      if(iH(I,J,m) == nmax)then ! end of "day":
        iH(I,J,m) = 0
        if(first_mod(I,J,m) == 1)then ! first averaging period only
          iD(I,J,m) = ID(I,J,m) + 1 
          do n=1,nmax
            DRA(I,J,iD(I,J,m),m) = DRA(I,J,iD(I,J,m),m) + HRA(I,J,n,m)
          end do
          DRA(I,J,iD(I,J,m),m) = DRA(I,J,iD(I,J,m),m)*bynmax
          if(iD(I,J,m) == nday_ch4(m))then !end first period
            PRS(I,J,m) = 0.d0
            do n=1,nday_ch4(m)
              PRS(I,J,m) = PRS(I,J,m) + DRA(I,J,n,m)
            end do
            avg_model(I,J,m)= PRS(I,J,m) * by_nday_ch4(m)
            first_mod(I,J,m)=0
            iD(I,J,m)=0
            i0(I,J,m)=0
          end if
        else ! not first averaging period: update the running average
          i0(I,J,m) = i0(I,J,m) + 1 ! move pointer
          if(i0(I,J,m)  ==  nday_ch4(m)+1) i0(I,J,m)=1
          temp=0.d0
          do n=1,nmax
            temp = temp + HRA(I,J,n,m) 
          end do
          temp = temp * bynmax ! i.e. today's average
          PRS(I,J,m) = PRS(I,J,m) - DRA(I,J,i0(I,J,m),m)
          DRA(I,J,i0(I,J,m),m) = temp
          PRS(I,J,m) = PRS(I,J,m) + DRA(I,J,i0(I,J,m),m)
          avg_model(I,J,m)= PRS(I,J,m) * by_nday_ch4(m)
        end if           
      end if

  
      END SUBROUTINE running_average
#endif

#ifdef INTERACTIVE_WETLANDS_CH4
      subroutine read_ncep_for_wetlands(end_of_day)
!@sum reads NCEP precip and temperature data and the coefficients
!@+ used to parameterize CH4 wetlands emissions from these. Keeps
!@+ running average of these. Calculated the portion to add to
!@+ the CH4 source to be used later in subroutine alter_wetlands_source.
!@auth Greg Faluvegi based on Jean Lerner
      USE RESOLUTION, only : im,jm
      use model_com, only: modelEclock
      USE DOMAIN_DECOMP_ATM, only: GRID, getDomainBounds, 
     &     am_i_root, write_parallel
      USE FILEMANAGER, only: openunit,closeunit
      USE TRCHEM_Shindell_COM, only: fix_CH4_chemistry
      USE TRACER_SOURCES, only: nday_ncep,by_nday_ncep,first_ncep,
     &iday_ncep,day_ncep,i0_ncep,avg_ncep,sum_ncep,nra_ncep,max_days,
     &PTBA,PTBA1,PTBA2,nncep,ncep_units,jmon_nc,jdlnc,nc_first
 
      implicit none
      
      integer :: n,m,i,j,k
      logical, intent(in) :: end_of_day
      character(len=300) :: out_line
      real*8 :: frac
      character*10, dimension(nncep) :: ncep_files =
     & (/'PREC_NCEP ','TEMP_NCEP ','BETA_NCEP ','ALPHA_NCEP'/)
      logical, dimension(nncep) :: ncep_bins =
     & (/.true.,.true.,.true.,.true./)
      real*8,dimension(GRID%I_STRT_HALO:GRID%I_STOP_HALO
     *     ,GRID%J_STRT_HALO:GRID%J_STOP_HALO,nra_ncep)::day_ncep_tmp

      INTEGER :: J_1, J_0, J_0H, J_1H, I_0, I_1

      call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1)
      call getDomainBounds(grid, I_STRT=I_0, I_STOP=I_1)

      if (fix_CH4_chemistry == 1) return ! don't bother if no CH4 chem

      by_nday_ncep(:)=1.d0/real(nday_ncep(:))

      do n=1,nra_ncep
       if(nday_ncep(n) > max_days .or. nday_ncep(n) < 1)
     & call stop_model('nday_ncep out of range',255)
      end do

! read the monthly data and interpolate to current day:   
      
      do k = 1,nncep
CCCCC   if(nc_first(k))then
          call openunit(ncep_files(k),ncep_units(k),ncep_bins(k))
CCCCC     nc_first(k)=.false. 
CCCCC   endif
        call read_mon_src_3(ncep_units(k),PTBA(:,:,k),frac)
CCCCC   jdlnc(k) = jday ! not used at the moment...
        call closeunit(ncep_units(k))
      end do
      write(out_line,*)
     &'NCEP for wetlands interpolated to current day',frac
      call write_parallel(trim(out_line))

! Then update running averages of NCEP Precip(m=1) & Temp(m=2)

      IF(.not. end_of_day) return ! avoid redundant accumulation on restarts

      do m=1,nra_ncep
        if(m > nra_ncep)call stop_model('check on ncep index m',255)
        if(first_ncep(m) == 1)then !accumulate first nday_ncep(m) days
          iday_ncep(m) = iday_ncep(m) + 1
          day_ncep(I_0:I_1,J_0:J_1,iday_ncep(m),m)=
     &    PTBA(I_0:I_1,J_0:J_1,m)
          if(iday_ncep(m) == nday_ncep(m))then !end of averaging period
            sum_ncep(I_0:I_1,J_0:J_1,m)=0.d0
            do n=1,nday_ncep(m)
              sum_ncep(I_0:I_1,J_0:J_1,m)=
     &        sum_ncep(I_0:I_1,J_0:J_1,m)+day_ncep(I_0:I_1,J_0:J_1,n,m)
            end do
            first_ncep(m)= 0
            iday_ncep(m) = 0
            i0_ncep(m)   = 0
            avg_ncep(I_0:I_1,J_0:J_1,m)=
     &      sum_ncep(I_0:I_1,J_0:J_1,m)*by_nday_ncep(m)
          end if
        else                     ! no longer first averaging period
          i0_ncep(m) = i0_ncep(m) + 1
          if(i0_ncep(m) > nday_ncep(m)) i0_ncep(m) = 1
          day_ncep_tmp(I_0:I_1,J_0:J_1,m) = PTBA(I_0:I_1,J_0:J_1,m)
          sum_ncep(I_0:I_1,J_0:J_1,m)= sum_ncep(I_0:I_1,J_0:J_1,m) - 
     &    day_ncep(I_0:I_1,J_0:J_1,i0_ncep(m),m)
          day_ncep(I_0:I_1,J_0:J_1,i0_ncep(m),m) = 
     &    day_ncep_tmp(I_0:I_1,J_0:J_1,m)
          sum_ncep(I_0:I_1,J_0:J_1,m)= sum_ncep(I_0:I_1,J_0:J_1,m) + 
     &    day_ncep(I_0:I_1,J_0:J_1,i0_ncep(m),m)
          avg_ncep(I_0:I_1,J_0:J_1,m) = sum_ncep(I_0:I_1,J_0:J_1,m)*
     &    by_nday_ncep(m)
        endif
      end do ! m

      return
      end subroutine read_ncep_for_wetlands


      subroutine alter_wetlands_source(n,ns_wet)
!@sum alter_wetlands_source changes the magnitude of the CH4 wetlands
!@+ (+ Tundra) source based on B.Walter's parameterization comparing
!@+ 1st layer ground temperature from 1 week ago and precipitation 
!@+ from 2 weeks ago, such that:
!@+  CH4emis = CH4emis + (alpha*TempAnom + beta*PrecAnom)
!@+ It then optionally allows for change in the wetlands source IJ
!@+ distribution based on criteria for topography, surf. air temp.,
!@+ downward SW radiation, land fraction, and ground wetness.
!@+ You also have the option to exclude the U.S. and E.U. from 
!@+ wetland distribution changes.
!@auth Greg Faluvegi
      USE RESOLUTION, only : im,jm
      USE FLUXES, only : fland
      USE MODEL_COM, only: itime, modelEclock
      USE DOMAIN_DECOMP_ATM, only: GRID, getDomainBounds,
     &   broadcast, write_parallel
#ifdef CUBED_SPHERE
      USE DD2D_UTILS, only : 
#else
      USE DOMAIN_DECOMP_1D, only : 
#endif
     &     pack_data
      use OldTracer_mod, only: itime_tr0,trname
      USE TRACER_COM, only: sfc_src,ntsurfsrcmax
      use TRACER_SOURCES, only: PTBA,nncep,first_ncep,avg_ncep,
     &   avg_model,nra_ncep,int_wet_dist,topo_lim,sat_lim,
     &   gw_ulim,gw_llim,SW_lim,exclude_us_eu,nra_ch4,first_mod,
     &   n__temp,n__sw,n__gwet,n__SAT,nn_or_zon,ice_age,add_wet_src
      use GEOM, only : lat2d_dg, lon2d_dg, imaxj
      use ghy_com, only : top_dev_ij,fearth
      USE TRCHEM_Shindell_COM, only: fix_CH4_chemistry

      implicit none
      
      integer, intent(in) :: n,ns_wet
      integer i,j,nt,iu,k
      character(len=300) :: out_line
      integer m,ii,ix,jj
      real*8 :: zm,zmcount
#ifdef CUBED_SPHERE
      real*8 :: src_glob(IM,IM,6),src_flatglob(1-im:2*im,1-im:2*im)
#else
      real*8, dimension(IM,JM) :: src_glob
#endif
      INTEGER :: J_1, J_0, J_0H, J_1H, I_0, I_1, J_0S,J_1S

      call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1,
     &     J_STRT_SKP=J_0S, J_STOP_SKP=J_1S)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

      if(ns_wet < 0 .or. ns_wet > ntsurfsrcmax)call stop_model
     & ('problem with ns_wet parameter',255)

! Don't alter the wetlands if:
      ! -->  the tracer is not supposed to be on yet:
      if (itime < itime_tr0(n)) return
      ! --> the methane is supposed to be fixed value:
      if (fix_CH4_chemistry == 1) return
      ! --> the first averaging period isn't done for ncep vars:
      do m=1,nra_ncep; if(first_ncep(m)==1) RETURN; end do

! Otherwise calculate the magnitude adjustments (not at poles or
! over water (/ice?) though. And not until the first averaging
! period is through (first_mod criterion):
      do m=1,nra_ncep
        loop_j: do j=J_0S,J_1S ! skipping poles
          loop_i: do i=I_0,imaxj(j)
            add_wet_src(i,j)=0.d0
            if(fearth(i,j) <= 0.) cycle loop_i
            if(first_mod(i,j,m)/=1.and.sfc_src(i,j,n,ns_wet)/=0.)
     &      add_wet_src(i,j) = add_wet_src(i,j) + PTBA(i,j,m+nra_ncep)*
     &      (avg_model(i,j,m) - avg_ncep(i,j,m))
          end do loop_i
        end do loop_j
      end do

! Now, determine the distribution (spatial) adjustments:
      if(int_wet_dist > 0)then ! if option is on

! Determine if there are new wetlands for given point,
! or remove existing wetlands:

       ! for nearest-neighbor all processors must know global source:
       if(nn_or_zon==0)then 
         call pack_data(grid,sfc_src(:,:,n,ns_wet),src_glob)
         call broadcast(grid,src_glob)
#ifdef CUBED_SPHERE
         ! Reorient the global array to the i-j index space of this
         ! processor. Implicit assumption (for now): wetlands exist on
         ! at least 2 faces; no data needed from the opposing face.
         call flatten_cube(src_glob,src_flatglob,im,grid%tile)
#endif
       else
#ifdef CUBED_SPHERE
         call stop_model('alter_wetlands_source: zonal mean '//
     &        'to be implemented using zonalmean_ij2ij',255)
#endif
       end if

       do j=J_0S,J_1S ! skipping poles
       do i=I_0,imaxj(j)
         if(fearth(i,j) <= 0.) cycle ! and boxes with no land

         ! if either {(point is not in U.S. or E.U.) .or. (it is but 
         ! the exclusion of U.S. and E.U. wetlands is OFF)} AND
         ! (the point has some land):

         if((exclude_us_eu == 0 .OR. .NOT.((lon2d_dg(i,j)
     &   >= -122.5.and.lon2d_dg(i,j) <= -72.5.and.lat2d_dg(i,j) >= 34.0
     &   .and.lat2d_dg(i,j) <= 46.5).or.(lon2d_dg(i,j) >= -12.5.and.
     &   lon2d_dg(i,j) <= 17.5.and.lat2d_dg(i,j) >= 37.5
     &        .and.lat2d_dg(i,j)<= 62.0))) .AND. (fland(i,j) > 0.))then

           ! if the topography slope, surface are temp., SW radiation,
           ! and ground wetness are within certain limits:

           if(top_dev_ij(i,j) < topo_lim .AND. avg_model(i,j,n__SAT)
     &     > sat_lim .AND. avg_model(i,j,n__SW) > SW_lim .AND.
     &     avg_model(i,j,n__gwet) > gw_llim .AND.
     &     avg_model(i,j,n__gwet) < gw_ulim) then

             ! if no wetlands there yet:
             
             if(sfc_src(i,j,n,ns_wet) == 0.)then
               zm=0.d0; zmcount=0.d0
               select case (nn_or_zon)
               case(1) ! use zonal average (of existing wetlands)
                 do ii=I_0, I_1  !EXCEPTIONAL CASE? 
                  if(sfc_src(ii,j,n,ns_wet) > 0.d0)then
                    zm=zm+sfc_src(ii,j,n,ns_wet)
                    zmcount=zmcount+1.d0
                  endif
                 enddo
               case(0) ! use nearest neighbor approach
                 ix=0
                 do while(zmcount == 0.)
                   ix=ix+1
#ifdef CUBED_SPHERE
                   if(ix>2*im)
     &                  call stop_model('ix>2*im int wetl dist',255)
                   do jj=j-ix,j+ix
                     if(jj<1-im .or. jj>2*im) cycle
                     do ii=i-ix,i+ix,1+(2*ix-1)*(1-abs((jj-j)/ix))
                       if(ii<1-im .or. ii>2*im) cycle
                       if(src_flatglob(ii,jj).le.0.) cycle
                       zm=zm+src_flatglob(ii,jj)
                       zmcount=zmcount+1.d0
                     enddo
                   enddo
#else
                   if(ix>im)call stop_model('ix>im int wetl dist',255)
                   do ii=i-ix,i+ix ;do jj=j-ix,j+ix
                     if(ii>0.and.ii<=im.and.jj>0.and.jj<=jm .and.
     &               src_glob(ii,jj) > 0.)then
                       zm=zm+src_glob(ii,jj)
                       zmcount=zmcount+1.d0
                     endif
                   enddo           ;enddo
#endif
                 enddo
               case default
                 call stop_model('problem with nn_or_zon',255)
               end select
               if(zmcount <= 0.)then
                 write(out_line,*)'zmcount for wetl src <= 0 @ IJ=',I,J
                 call write_parallel(trim(out_line),unit=6,crit=.true.)
                 add_wet_src(i,j)=-1.d0*sfc_src(i,j,n,ns_wet)
               else
                 add_wet_src(i,j)=zm/zmcount
               end if
             end if       ! end of no prescribed wetlands there
           else           ! no wetlands should be here
             add_wet_src(i,j)=-1.d0*sfc_src(i,j,n,ns_wet)
           end if         ! end wetlands criteria
         end if           ! end U.S./E.U./land criteria
       end do             ! i loop
       end do             ! j loop
      end if              ! Consider interactive welands distrbution?

! Limit wetlands source to be positive:
      do j=J_0,J_1  
        do i=I_0,imaxj(j)
          if((sfc_src(i,j,n,ns_wet)+add_wet_src(i,j))<0.)
     &    add_wet_src(i,j)=-1.d0*sfc_src(i,j,n,ns_wet)
        enddo
      enddo

! Optionally disallow emissions over glacier latitudes:
      if(ice_age /= 0.) then
        do j=J_0,J_1 
        do i=I_0,imaxj(j)
          if(abs(lat2d_dg(i,j))>abs(ice_age))
     &    add_wet_src(i,j)=-1.d0*sfc_src(i,j,n,ns_wet)
        enddo
        enddo
      endif

      return
      end subroutine alter_wetlands_source

      subroutine flatten_cube(arr6,flat5,n,face)
!@sum Places data from 5 of the 6 faces of a cube into the i-j index
!@+   space of one of the faces. Referenced to that face, the cube is
!@+   comprised of
!@+   (1) that face
!@+   (2) the 4 adjacent faces in the directions left, right, down, up
!@+   (3) the opposing face
!@+   Data from the local/adjacent faces are stored in
!@+   flat5(1-n:2*n,1-n:2*n) as follows:
!@+    local face: i, j = 1  :n  , 1  :n  (each face is n by n points)
!@+   to the left: i, j = 1-n:0  , 1  :n
!@+         right: i, j = 1+n:2*n, 1  :n
!@+          down: i, j = 1  :n  , 1-n:0
!@+            up: i, j = 1  :n  , 1+n:2*n
!@+   flat5 is set to zero where both i and j are either <1 or >n
!@+   Data from the opposing side can also be placed into the flat array
!@+   in a symmetric manner, but this is not necessary for current
!@+   applications (filling missing data using neighboring values).
!@auth M. Kelley
      implicit none
!@ var arr6  : input global array with 3d-cube indexing
!@ var flat5 : output global array with flattened-cube indexing
!@ var n     : faces are n by n gridpoints
!@ var face  : the cube face to be used as reference for i-j indexing
      integer :: n,face
      real*8, dimension(n,n,6) :: arr6
      real*8, dimension(1-n:2*n,1-n:2*n) :: flat5
      integer :: facem1,facep1,facem2,facep2
      facem1 = 1+mod(6+(face-1)-1,6)
      facep1 = 1+mod(6+(face-1)+1,6)
      facem2 = 1+mod(6+(face-1)-2,6)
      facep2 = 1+mod(6+(face-1)+2,6)
      flat5 = 0.
      flat5(1:n,1:n) = arr6(:,:,face)
      if(mod(face,2).eq.0) then
        flat5(1-n:0,1:n)   = arr6(:,:,facem1)           ! l
        flat5(n+1:2*n,1:n) = rotcw(arr6(:,:,facep2),n)  ! r
        flat5(1:n,1-n:0) = rotccw(arr6(:,:,facem2),n)   ! d
        flat5(1:n,n+1:2*n) = arr6(:,:,facep1)           ! u
      else
        flat5(1-n:0,1:n) = rotcw(arr6(:,:,facem2),n)    ! l
        flat5(n+1:2*n,1:n) = arr6(:,:,facep1)           ! r
        flat5(1:n,1-n:0) = arr6(:,:,facem1)             ! d
        flat5(1:n,n+1:2*n) = rotccw(arr6(:,:,facep2),n) ! u
      endif
      return
      contains
      function rotcw(arr,m)
      integer :: m
      real*8, dimension(m,m) :: arr,rotcw
      rotcw = transpose(arr(m:1:-1,:))
      end function rotcw
      function rotccw(arr,m)
      integer :: m
      real*8, dimension(m,m) :: arr,rotccw
      rotccw = transpose(arr(:,m:1:-1))
      end function rotccw
      end subroutine flatten_cube


      SUBROUTINE read_mon_src_3(iu,data,frac)
!@sum Read in monthly sources and interpolate to current day
! I know... yet another copy of this kind of routine...
! I have them all combined in one, but there is a bug so I
! can't commit it yet...
!@auth Greg Faluvegi, Jean Lerner and others

      USE FILEMANAGER, only : NAMEUNIT
      USE DOMAIN_DECOMP_ATM, only : GRID, getDomainBounds, 
     & AM_I_ROOT,write_parallel,
     & READT_PARALLEL, REWIND_PARALLEL, BACKSPACE_PARALLEL
      USE RESOLUTION, only : im,jm
      use model_com, only: modelEclock
      USE MODEL_COM, only: idofm=>JDmidOfM

      implicit none

      real*8, DIMENSION(GRID%I_STRT_HALO:GRID%I_STOP_HALO
     *     ,GRID%J_STRT_HALO:GRID%J_STOP_HALO) ::tlca,tlcb,data
      real*8 :: frac,alpha
      integer ::  imon,iu,k
      character(len=300) :: out_line

      integer :: J_0, J_1, I_0, I_1

      call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1)
      call getDomainBounds(grid, I_STRT=I_0, I_STOP=I_1)

      imon=1
      if (modelEclock%dayOfYear() <= 16)  then ! JDAY in Jan 1-15, first month is Dec
        call readt_parallel(grid,iu,nameunit(iu),tlca,12)
        call rewind_parallel( iu )
      else            ! JDAY is in Jan 16 to Dec 16, get first month
        do while(modelEclock%dayOfYear() > idofm(imon) .AND. imon <= 12)
          imon=imon+1
        enddo
        call readt_parallel(grid,iu,nameunit(iu),tlca,imon-1)
        if (imon == 13)then
          call rewind_parallel( iu )
        endif
      end if
      call readt_parallel(grid,iu,nameunit(iu),tlcb,1)

c**** Interpolate two months of data to current day
      frac = float(idofm(imon)-modelEclock%dayOfYear()) / 
     & (idofm(imon)-idofm(imon-1))
      data(I_0:I_1,J_0:J_1)=tlca(I_0:I_1,J_0:J_1)*frac + 
     & tlcb(I_0:I_1,J_0:J_1)*(1.-frac)

      return
      end subroutine read_mon_src_3
#endif

