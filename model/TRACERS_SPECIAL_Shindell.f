#include "rundeck_opts.h"

      MODULE TRACER_SOURCES

#ifdef TRACERS_ON

      USE TRACER_COM

      IMPLICIT NONE
      SAVE

!@param Laircr the number of layers of aircraft data read from file
!@param Lsulf the number of layers of sulfate SA data read from file
      INTEGER, PARAMETER :: 
     &                      Laircr      =19,
     &                      Lsulf       =23 ! not LM
#ifdef GFED_3D_BIOMASS 
     &                     ,LbbGFED     =6
      real*8, allocatable, dimension(:,:,:,:) :: GFED_BB
#endif
#ifdef SHINDELL_STRAT_EXTRA
      REAL*8, PARAMETER ::  GLTic = 1.d-9 ! pppv
#endif   

#ifdef INTERACTIVE_WETLANDS_CH4
!@dbparam nn_or_zon approach to use for expanding wetlands 1=
!@+ zonal average, 0=nearest neighbor average
!@dbparam int_wet_dist to turn on/off interacive SPATIAL wetlands
!@dbparam ice_age if not = 0 allows no wetl emis for J >= ice_age
!@dbparam exclude_us_eu to exclude (=1) the U.S. and E.U. from 
!@+ interactive wetland distributiont
!@dbparam topo_lim upper limit on topography for int wetl dist
!@dbparam sat_lim lower limit on surf air temp for int wetl dist
!@dbparam gw_llim lower limit on ground wetness for int wetl dist
!@dbparam gw_ulim upper limit on ground wetness for int wetl dist
!@dbparam SW_lim lower limit on SW down flux for int wetl dist
!@param nra_ch4 number of running averages needed for int-wetlands
!@param nra_ncep # of ncep running averages needed for int-wetlands
!@param maxHR_ch4 maximum number of sub-daily accumulations
!@param fact_ncep for temp: nothing, for prec: 24 converts from model
!@+     units kg/m2/hr to mm/day
!@param nday_ncep number of days in ncep running averages
!@param nday_ch4 number of days in running average
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
      ! next line: 1800=dtsrc, 1=nisurf, but since those are not
      ! parameters, I don't know how to soft-code this. There is a 
      ! failsafe in the read_CH4_sources routine, however.
      integer, parameter :: nra_ch4 = 5, maxHR_ch4=24*1*3600/1800,
     & nra_ncep=2, n__prec=1, n__temp=2, n__SW=3, n__SAT=4, n__gwet=5,
     & max_days=28
      integer :: int_wet_dist=0,exclude_us_eu=1,nn_or_zon=0,ice_age=0
      real*8 :: topo_lim = 205.d0, sat_lim=-9.d0,
     & gw_ulim=100.d0, gw_llim=18.d0, SW_lim=27.d0
      real*8, parameter, dimension(nra_ncep) ::fact_ncep=(/24.d0,1.d0/)
      integer, parameter, dimension(nra_ch4) :: 
     &                                  nday_ch4=(/28,14,14,14,28/)
      integer, parameter, dimension(nra_ncep):: nday_ncep=(/28,14/)
      real*8, dimension(nra_ch4)             :: by_nday_ch4
      real*8, dimension(nra_ncep)            :: by_nday_ncep
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:):: day_ncep
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:):: DRA_ch4
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:)  :: avg_model,PRS_ch4
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:)  :: avg_ncep,sum_ncep
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:):: HRA_ch4
      integer, dimension(nra_ncep)           :: iday_ncep, i0_ncep,
     &                                          first_ncep
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: iHch4,iDch4,i0ch4,
     &                                          first_mod
#endif
#endif
      END MODULE TRACER_SOURCES



      subroutine alloc_tracer_sources(grid)
!@SUM  To alllocate arrays whose sizes now need to be determined
!@+    at run-time
!@auth G.Faluvegi
!@ver  1.0
      use domain_decomp, only : dist_grid, get, write_parallel
      use model_com, only     : im
      use tracer_com, only : ntm
#ifdef INTERACTIVE_WETLANDS_CH4
      use TRACER_SOURCES, only: maxHR_ch4,
     & day_ncep,DRA_ch4,avg_model,PRS_ch4,avg_ncep,sum_ncep,HRA_ch4,
     & iHch4,iDch4,i0ch4,first_mod,nra_ch4,nra_ncep,max_days
#endif
#ifdef GFED_3D_BIOMASS
      use TRACER_SOURCES, only: GFED_BB,LbbGFED
#endif

      IMPLICIT NONE

      type (dist_grid), intent(in) :: grid
      integer :: ier, J_1H, J_0H
      logical :: init = .false.

      if(init)return
      init=.true.
    
      call get( grid , J_STRT_HALO=J_0H, J_STOP_HALO=J_1H )
 
#ifdef GFED_3D_BIOMASS
      allocate( GFED_BB(IM,J_0H:J_1H,LbbGFED,ntm) )
#endif
#ifdef INTERACTIVE_WETLANDS_CH4
      allocate( first_mod(IM,J_0H:J_1H,nra_ch4) )
      allocate( iHch4(IM,J_0H:J_1H,nra_ch4) )
      allocate( iDch4(IM,J_0H:J_1H,nra_ch4) )
      allocate( i0ch4(IM,J_0H:J_1H,nra_ch4) )
      allocate( day_ncep(IM,J_0H:J_1H,max_days,nra_ncep) )
      allocate( DRA_ch4(IM,J_0H:J_1H,max_days,nra_ch4) )
      allocate( avg_model(IM,J_0H:J_1H,nra_ch4) )
      allocate( PRS_ch4(IM,J_0H:J_1H,nra_ch4) )
      allocate( avg_ncep(IM,J_0H:J_1H,nra_ncep) )
      allocate( sum_ncep(IM,J_0H:J_1H,nra_ncep) )
      allocate( HRA_ch4(IM,J_0H:J_1H,maxHR_ch4,nra_ch4) )
#endif      
      
      return
      end subroutine alloc_tracer_sources 
            
      
      MODULE LIGHTNING
!@sum  LIGHTNING_COM model variables lightning parameterization
!@auth Colin Price (modelEification by Greg Faluvegi)
!@ver  1.0 (taken from CB436Tds3M23)
#ifdef TRACERS_ON
      USE MODEL_COM, only : IM,JM,LM,LS1

      IMPLICIT NONE
      SAVE
      
!@var JN J at 30 N
!@var JS J at 30 S
!@var I dummy
      INTEGER, PARAMETER :: JS = JM/3 + 1, JN = 2*JM/3
      INTEGER I
      REAL*8, ALLOCATABLE, DIMENSION(:,:)  :: RNOx_lgt
      REAL*8, DIMENSION(LM) :: SRCLIGHT
      
!@var HGT Pickering vertical lightning distributions (1998)
      REAL*8 HGT(2,2,16)
      DATA (HGT(1,1,I),I=1,16)/8.2,1.9,2.1,1.6,1.1,1.6,3.0,5.8,
     &  7.6,9.6,10.5,12.3,11.8,12.5,8.1,2.3/
      DATA (HGT(1,2,I),I=1,16)/5.8,2.9,2.6,2.4,2.2,2.1,2.3,6.1,
     & 16.5,14.1,13.7,12.8,12.5,2.8,0.9,0.3/
      DATA (HGT(2,1,I),I=1,16)/20.1,2.3,0.8,1.5,3.4,5.3,3.6,3.8,
     &    5.4,6.6,8.3,9.6,12.8,10.0,6.2,0.3/
      DATA (HGT(2,2,I),I=1,16)/5.8,2.9,2.6,2.4,2.2,2.1,2.3,6.1,
     & 16.5,14.1,13.7,12.8,12.5,2.8,0.9,0.3/ 
 
      END MODULE LIGHTNING 
      
      
          
      subroutine alloc_lightning(grid)
!@SUM  To alllocate arrays whose sizes now need to be determined
!@+    at run-time
!@auth G.Faluvegi
!@ver  1.0
      use domain_decomp, only : dist_grid, get, write_parallel
      use LIGHTNING, only     : RNOx_lgt
      use model_com, only     : im

      IMPLICIT NONE

      type (dist_grid), intent(in) :: grid
      integer :: ier, J_1H, J_0H
      logical :: init = .false.

      if(init)return
      init=.true.
    
      call get( grid , J_STRT_HALO=J_0H, J_STOP_HALO=J_1H )
 
      allocate( RNOx_lgt(IM,J_0H:J_1H) )
      
      return
      end subroutine alloc_lightning
      
      
      
#ifdef SHINDELL_STRAT_EXTRA
      subroutine overwrite_GLT
!@sum L=1 overwriting of generic linear tracer    
!@auth Greg Faluvegi
C****
C**** Right now, there is just one L=1 source that changes 
C**** linearly in time (at 1% increase per year)

      USE MODEL_COM, only: itime,itime0,DTsrc,im,jm
      USE GEOM, only: dxyp,IMAXJ  
      USE DYNAMICS, only: am
      USE TRACER_COM, only: trname,trm,n_GLT
      USE TRACER_SOURCES, only: GLTic
      USE FLUXES, only : tr3Dsource
      USE TRCHEM_Shindell_COM,only: bymass2vol
      USE DOMAIN_DECOMP, ONLY : GET, grid, write_parallel
      
      IMPLICIT NONE
      
      INTEGER :: J_0, J_1
           
!@var by_s_in_yr recip. of # seconds in a year
!@var new_mr mixing ratio to overwrite in L=1 this time step
!@var new_mass mass to overwrite in L=1 this time step

      REAL*8 bydtsrc, by_s_in_yr, new_mr, new_mass
      integer i,j

      CALL GET(grid, J_STRT=J_0, J_STOP=J_1)

      bydtsrc=1.d0/DTsrc
      by_s_in_yr = 1.d0/(365.d0*24.d0*60.d0*60.d0)

C initial source is an overwriting of GLTic pppv, then add
C 1% every year, linearly in time. (note: bymass2vol should
C just be 1 for this tracer, but kept it in here, in case
C we change that.)
      new_mr = GLTic*(1.d0 + Itime*DTsrc*by_s_in_yr*1.d-2)! pppv
      do j=J_0,J_1; do i=1,imaxj(j)
        new_mass=new_mr*bymass2vol(n_GLT)*am(1,i,j)*DXYP(j) ! kg
        tr3Dsource(i,j,1,1,n_GLT)=(new_mass-trm(i,j,1,n_GLT))*bydtsrc
        !i.e. tr3Dsource in kg/s 
      end do   ; end do

      return
      end subroutine overwrite_GLT
#endif

#ifdef hell_freezes_over
! keeping this in for reference for now (regarding the
! wetlands portion, which must be reinserted somewhere):

      subroutine read_CH4_sources(nt,iact)
!@sum reads in CH4 surface sources and sinks
!@auth Jean Lerner
C****
C**** There are 3 monthly sources and 11 annual sources
C**** Annual sources are read in at start and re-start of run only
C**** Monthly sources are interpolated each day
      USE MODEL_COM, only: itime,JDperY,im,jm,jday,focean,NIsurf
     &                     ,DTsrc,fland,itimei
      USE DOMAIN_DECOMP, only: GRID, GET, readt_parallel, pack_data,
     &                         unpack_data, am_i_root, write_parallel
      USE LAKES_COM, only: flake
      USE CONSTANT, only: sday
      USE FILEMANAGER, only: openunit,closeunit, nameunit
      USE TRACER_COM, only: itime_tr0,trname
      use TRACER_SOURCES, only: ch4_src,nCH4src
#ifdef INTERACTIVE_WETLANDS_CH4      
     &       ,nday_ncep,by_nday_ncep,first_ncep,iday_ncep,day_ncep
     &       ,i0_ncep,avg_ncep,sum_ncep,fact_ncep,avg_model,nra_ncep
     &       ,max_days,int_wet_dist,topo_lim,sat_lim,gw_ulim,gw_llim
     &       ,SW_lim,exclude_us_eu,nra_ch4,first_mod,n__temp,n__sw
     &       ,n__gwet,n__SAT,nn_or_zon,nday_ch4,ice_age
#endif
      use GEOM, only : lat_dg, lon_dg, imaxj
      use ghy_com, only : top_dev_ij,fearth
      USE TRCHEM_Shindell_COM, only: PI_run,PIratio_indus,PIratio_bburn
     & ,fix_CH4_chemistry
      
      implicit none
      
!@var nanns,nmons: number of annual and monthly input files
!@var PIfact factor for altering source to preindustrial (or not)
!@var adj Factors that tune the total amount of individual sources         

      integer, parameter :: nanns=11,nmons=3,kwet=14
      integer, dimension(nanns-3) :: ann_units
      integer, dimension(nmons) :: mon_units, imon
      integer :: jdlast=0
      integer i,j,nt,iact,iu,k
      character*80 :: title
      character*12, dimension(nanns-3) :: ann_files =
     *  (/'CH4_ANIMALS ','CH4_COALMINE','CH4_GASLEAK ','CH4_GASVENT ',
     *    'CH4_CITYDUMP','CH4_SOIL_ABS','CH4_TERMITES','CH4_COALBURN'/)
      character*8, dimension(nmons) :: mon_files =
     *  (/'CH4_BURN','CH4_RICE','CH4_WETL'/)
      character(len=300) :: out_line
      logical, dimension(nanns-3) :: ann_bins=
     * (/.true.,.true.,.true.,.true.,.true.,.true.,.true.,.true./)
      logical :: ifirst=.true.
      logical, dimension(nmons) :: mon_bins=(/.true.,.true.,.true./)
      real*8 :: PIfact, frac 
      real*8, dimension(nanns+nmons) :: adj=
     *  (/1.3847,1.0285,3.904,1.659,1.233,1.194,0.999,
     *  7.2154, 3.5997e-5,17.330e-4,5.3558e-5,0.4369,0.7533,0.9818/)
      real*8, allocatable, dimension(:,:,:) :: tlca, tlcb  ! for monthly sources
      
#ifdef INTERACTIVE_WETLANDS_CH4
!@var PTBA variable to hold the pressure, temperature, beta, and alpha
!@var PTBA1, PTBA2 for interpolations of PTBA
!@var day_ncep_tmp temp array for computing running average
!@var zm,zmcount temp variables for computing averages
      integer, parameter :: nncep=4
      integer, dimension(nncep):: ncep_units,jmon
      integer n,m,ii,ix,jj
      character*10, dimension(nncep) :: ncep_files =
     & (/'PREC_NCEP ','TEMP_NCEP ','BETA_NCEP ','ALPHA_NCEP'/)
      logical, dimension(nncep) :: ncep_bins = 
     & (/.true.,.true.,.true.,.true./)
      real*8, allocatable, dimension(:,:,:) :: PTBA1, PTBA2
      real*8,dimension(im,GRID%J_STRT_HALO:GRID%J_STOP_HALO,nncep) :: 
     &                                                            PTBA
      real*8,dimension(im,GRID%J_STRT_HALO:GRID%J_STOP_HALO,nra_ncep)::
     &                                                    day_ncep_tmp
      real*8 :: zm,zmcount
      real*8, dimension(IM,JM) :: src_glob
#endif
      save ifirst,jdlast,tlca,tlcb,mon_units,imon
#ifdef INTERACTIVE_WETLANDS_CH4
     & ,jmon,PTBA1,PTBA2,ncep_units
#endif

      INTEGER :: J_1, J_0, J_0H, J_1H

      CALL GET(grid, J_STRT=J_0, J_STOP=J_1)

      if (itime < itime_tr0(nt)) return
      if (fix_CH4_chemistry == 1) return

#ifdef INTERACTIVE_WETLANDS_CH4
      by_nday_ncep(:)=1.d0/real(nday_ncep(:))
C     some sanity checks:
      if(nisurf /= 1) call stop_model
     &('nisurf no longer=1, please alter maxHR_ch4',255)
      if(DTsrc /= 1800.) call stop_model
     &('DTsrc no longer=1800, please alter maxHR_ch4',255)
      do n=1,nra_ch4
       if(nday_ch4(n) > max_days .or. nday_ch4(n) < 1)
     & call stop_model('nday_ch4 out of range',255)
      end do
      do n=1,nra_ncep
       if(nday_ncep(n) > max_days .or. nday_ncep(n) < 1)
     & call stop_model('nday_ncep out of range',255)
      end do
#endif
C****
C**** Annual Sources and sinks
C**** Apply adjustment factors to bring sources into balance
C**** Annual sources are in KG C/M2/Y
C**** Sources need to be kg/m^2 s; convert /year to /s
C****
      if (ifirst) then
        call GET(grid, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)
        allocate(tlca(im,j_0H:j_1H,nmons),tlcb(im,j_0H:j_1H,nmons))
#ifdef INTERACTIVE_WETLANDS_CH4
        allocate(ptba1(im,j_0H:j_1H,nncep),ptba2(im,j_0H:j_1H,nncep))
#endif
        do k = 1,nanns-3  
          call openunit(ann_files(k),ann_units(k),ann_bins(k))
          iu=ann_units(k)
          select case(PI_run)
          case(1) ! pre-industrial
            select case(k)
            case(1)        ; PIfact=0.18d0 ! animals
            case(6)        ; PIfact=0.38d0 ! soil abs
            case(2,3,4,5,8); PIfact=PIratio_indus
            case default; PIfact=1.d0
            end select
          case default
            PIfact=1.d0
          end select
          call readt_parallel(grid,iu,nameunit(iu),0,CH4_src(:,:,k),1)
          CH4_src(:,J_0:J_1,k)=CH4_src(:,J_0:J_1,k)*PIfact*adj(k)/
     &    (sday*JDperY)
C REDUCE SOME SOURCES MANUALLY:
          if(ann_files(k) == 'CH4_CITYDUMP')
     &    CH4_src(:,J_0:J_1,k)=CH4_src(:,J_0:J_1,k)*7.5d-1
          if(ann_files(k) == 'CH4_GASLEAK ')
     &    CH4_src(:,J_0:J_1,k)=CH4_src(:,J_0:J_1,k)*8.0d-1
          if(ann_files(k) == 'CH4_ANIMALS ')
     &    CH4_src(:,J_0:J_1,k)=CH4_src(:,J_0:J_1,k)*9.0d-1
C-----------------------------
          call closeunit(iu)
        end do
        ! 3 miscellaneous sources
        k=nanns-2
        CH4_src(:,J_0:J_1,k) = focean(:,J_0:J_1)*adj(k)/(sday*JDperY)
        k=nanns-1
        CH4_src(:,J_0:J_1,k) =  flake(:,J_0:J_1)*adj(k)/(sday*JDperY)
        k=nanns
        CH4_src(:,J_0:J_1,k) = fearth(:,J_0:J_1)*adj(k)/(sday*JDperY)
        ifirst = .false.
      endif
C****
C**** Monthly sources are interpolated to the current day
C****
C**** Also, Apply adjustment factors to bring sources into balance
C**** Monthly sources are in KG C/M2/S => src in kg/m^2 s
C****
      j = 0
      do k = nanns+1,nCH4src
        select case(PI_run)
        case(1) ! pre-industrial
          select case(k)
          case(nanns+1); PIfact=PIratio_bburn
          case(nanns+2); PIfact=PIratio_indus ! rice
#ifndef INTERACTIVE_WETLANDS_CH4
          case(nanns+3); PIfact=1.1d0         ! wetlands
#endif
          case default; PIfact=1.d0
          end select
        case default
          PIfact=1.d0
        end select
        j = j+1
        call openunit(mon_files(j),mon_units(j),mon_bins(j))
        call read_mon_src_2(mon_units(j),jdlast,
     *    tlca(:,:,j),tlcb(:,:,j),CH4_src(:,:,k),frac,imon(j))
        call closeunit(mon_units(j))
        CH4_src(:,J_0:J_1,k) = CH4_src(:,J_0:J_1,k)*adj(k)*PIfact
C REDUCE SOME SOURCES MANUALLY:
        if(mon_files(j) == 'CH4_RICE')
     &  CH4_src(:,J_0:J_1,k)=CH4_src(:,J_0:J_1,k)*7.5d-1
C-----------------------------
      end do
#ifdef INTERACTIVE_WETLANDS_CH4
       do k = 1,nncep
        call openunit(ncep_files(k),ncep_units(k),ncep_bins(k))
        call read_mon_src_2(ncep_units(k),jdlast,
     *    PTBA1(:,:,k),PTBA2(:,:,k),PTBA(:,:,k),frac,jmon(k))
        call closeunit(ncep_units(k))
      end do     
#endif
      jdlast = jday
      write(out_line,*)
     &trname(nt),'Sources interpolated to current day',frac
      call write_parallel(trim(out_line))   
C****
C**** Zonal adjustment for combined wetlands and tundra
C****
      do j=J_0,J_1
        if ((lat_dg(j,1) >= -30.) .and. (lat_dg(j,1) <= 30.))then
          CH4_src(:,j,kwet) = CH4_src(:,j,kwet)*1.761d0
        else
          CH4_src(:,j,kwet) = CH4_src(:,j,kwet)*0.6585d0
        endif
      end do
C****
C**** Also, increase the wetlands + tundra CH4 emissions:
C**** (limit positive):
      do j=J_0,J_1; do i=1,im
        CH4_src(i,j,kwet)=max(1.78d0*CH4_src(i,j,kwet),0.d0)
      end do      ; end do
C  
#ifdef INTERACTIVE_WETLANDS_CH4
C****
C**** Adjust the wetlands+tundra CH4 source based on 1st layer ground
C**** temperature from one week ago, precipitation from 2 weeks ago,
C**** and B. Walter's regression coefficients.  I.e.:
C**** 
C**** CH4emis = CH4emis + (alpha*TempAnom + beta*PrecAnom)
C****
C
C NCEP Precipitation(m=1) and Temperature(m=2) running averages:
C
      IF(Itime == ItimeI)first_ncep(:)=1      ! Initialize
      IF(IACT <= 0) return ! avoid redundant accumulation on restarts...
      do m=1,nra_ncep
        if(m > nra_ncep)call stop_model('check on ncep index m',255)
        if(first_ncep(m) == 1)then !accumulate first nday_ncep(m) days
          iday_ncep(m) = iday_ncep(m) + 1
          day_ncep(:,J_0:J_1,iday_ncep(m),m)=PTBA(:,J_0:J_1,m)
          if(iday_ncep(m) == nday_ncep(m))then !end of averaging period
            sum_ncep(:,J_0:J_1,m)=0.d0
            do n=1,nday_ncep(m)
              sum_ncep(:,J_0:J_1,m)=
     &        sum_ncep(:,J_0:J_1,m) + day_ncep(:,J_0:J_1,n,m)
            end do
            first_ncep(m)= 0
            iday_ncep(m) = 0
            i0_ncep(m)   = 0
            avg_ncep(:,J_0:J_1,m)=sum_ncep(:,J_0:J_1,m)*by_nday_ncep(m)
          end if
        else                     ! no longer first averaging period
          i0_ncep(m) = i0_ncep(m) + 1
          if(i0_ncep(m) > nday_ncep(m)) i0_ncep(m) = 1
C         UPDATE RUNNING AVERAGE:
          day_ncep_tmp(:,J_0:J_1,m) = PTBA(:,J_0:J_1,m)
          sum_ncep(:,J_0:J_1,m)=
     &    sum_ncep(:,J_0:J_1,m) - day_ncep(:,J_0:J_1,i0_ncep(m),m)
          day_ncep(:,J_0:J_1,i0_ncep(m),m) = day_ncep_tmp(:,J_0:J_1,m)
          sum_ncep(:,J_0:J_1,m)=
     &    sum_ncep(:,J_0:J_1,m) + day_ncep(:,J_0:J_1,i0_ncep(m),m)
          avg_ncep(:,J_0:J_1,m) = sum_ncep(:,J_0:J_1,m)*by_nday_ncep(m)
        endif
      end do ! m
C
C Don't alter source until enough statistics are built up:
      do m=1,nra_ncep; if(first_ncep(m)==1) RETURN; end do
      do j=J_0,J_1
        if(j == 1 .or. j==jm) cycle
        if(fearth(i,j) <= 0.) cycle
        do i=1,IMAXJ(j); do m=1,nra_ch4 
          if(first_mod(i,j,m) == 1) return
        end do         ; end do
      end do

C Otherwise, apply the magnitude adjustments:
      do m=1,nra_ncep
        loop_j: do j=J_0,J_1
          if(j == 1 .or. j==jm) cycle loop_j
          loop_i: do i=1,imaxj(j)
            if(fearth(i,j) <= 0.) cycle loop_i
            CH4_src(i,j,kwet) = CH4_src(i,j,kwet) +
     &      PTBA(i,j,m+nra_ncep) * (fact_ncep(m)*
     &      avg_model(i,j,m) - avg_ncep(i,j,m))
          end do loop_i
        end do loop_j
      end do
  
C Now, determine the distribution (spatial) adjustments:
      if(int_wet_dist > 0)then
C First, determine if there are new wetlands for given point:
! Conditions of these two IF statements: point is within the limits
! of topography, surface air temperature, downward SW radiation,
! and ground wetness, and either (1) the exclusion of U.S. and E.U.
! wetlands is turned off or (2) the point is outside those regions, 
! AND this box has some land in it.
       do j=J_0,J_1
       if(j == 1 .or. j==jm) cycle
       if(fearth(i,j) <= 0.) cycle
       do i=1,imaxj(j)
        if(exclude_us_eu == 0 .OR. .NOT.((lon_dg(i,1)
     &   >= -122.5.and.lon_dg(i,1) <= -72.5.and.lat_dg(j,1) >= 34.0
     &  .and.lat_dg(j,1) <= 46.5).or.(lon_dg(i,1) >= -12.5.and.
     &  lon_dg(i,1) <= 17.5.and.lat_dg(j,1) >= 37.5.and.lat_dg(j,1)
     &   <= 62.0)) .AND. fland(i,j) > 0.)then
         if(top_dev_ij(i,j) < topo_lim .AND. avg_model(i,j,n__SAT)
     &   > sat_lim .AND. avg_model(i,j,n__SW) > SW_lim .AND.
     &   avg_model(i,j,n__gwet) > gw_llim .AND. 
     &   avg_model(i,j,n__gwet) < gw_ulim) then
          if(CH4_src(i,j,kwet) <= 0.)then
           zm=0.d0; zmcount=0.d0
           select case (nn_or_zon)
           case(1) ! use zonal average
             do ii=1,im
              if(CH4_src(ii,j,kwet) > 0.d0)then
                zm=zm+CH4_src(ii,j,kwet)
                zmcount=zmcount+1.d0
              endif
             enddo
           case(0) ! use nearest neighbors
             call pack_data( grid, CH4_src(:,:,kwet), src_glob(:,:) )
             if(am_i_root( )) then             
               ix=0
               do while(zmcount == 0.)
                ix=ix+1
                if(ix>im)call stop_model('ix>im int wetl dist',255)
                do ii=i-ix,i+ix
                 do jj=j-ix,j+ix
                  if(ii > 0.and.ii <= im.and.jj > 0.and.jj <= jm.and.
     &            CH4_src(ii,jj,kwet) > 0.)then
                    zm=zm+CH4_src(ii,jj,kwet)
                    zmcount=zmcount+1.d0
                  endif
                 enddo
                enddo
               enddo   
             endif ! root
           case default
             call stop_model('problem with nn_or_zon',255)
           end select
           if(zmcount <= 0.)then
             write(out_line,*)'zmcount for wetl src <= 0 @ IJ=',I,J
             call write_parallel(trim(out_line),unit=666,crit=.true.) 
             CH4_src(i,j,kwet)=0.d0
           else
             CH4_src(i,j,kwet)=zm/zmcount
           end if
          end if          ! end of no prescribed wetlands there
         else             ! no wetlands should be here
          CH4_src(i,j,kwet)=0.d0
         end if           ! end wetlands criteria
        end if            ! end U.S./E.U./land criteria
       end do             ! i loop
       end do             ! j loop
      end if              ! Consider interactive welands distrbution?
      
C Limit wetlands source to be positive:
      CH4_src(:,J_0:J_1,kwet)=max(CH4_src(:,J_0:J_1,kwet),0.d0) 
C No emissions over glacier latitudes if desired:
      if(ice_age /= 0) then
        do j=MAX(J_0,ice_age),MIN(J_1,JM)
          CH4_src(:,j,kwet)=0.d0
        enddo
      endif
#endif
      return
      end subroutine read_CH4_sources
#endif

c
c
      SUBROUTINE get_lightning_NOx
!@sum  get_lightning_NOx to define the 3D source of NOx from lightning
!@auth Colin Price / Greg Faluvegi
!@ver  1.0 (based on CB436Tds3M23 & DB396Tds3M23)
c
C**** GLOBAL parameters and variables:
      USE FLUXES, only        : tr3Dsource
      USE TRACER_COM, only    : n_NOx,nOther
      USE LIGHTNING, only     : HGT,JS,JN,SRCLIGHT,RNOx_lgt
      USE CONSTANT, only      : bygrav
      USE MODEL_COM, only     : fland,IM,LS1,LM
      USE DYNAMICS, only      : LTROPO, PHI
      USE GEOM, only          : BYDXYP
      USE DOMAIN_DECOMP, only : GRID, GET, write_parallel
c
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
c
!@var LATINDX Pickering latitude index
!@var LANDINDX Pickering surface type index
!@var IH Pickering altitude index
!@var HEIGHT Pickering HGT variable specific to I,J,L, point
!@var LEVTROP local variable to hold the tropopause level
!@var ALTTROP altitude of tropopause
!@var ALTTOP altitude at the top of ?
!@var I,J,L loop indicies in x,y,z directions
!@param pmin2psec to convert from min-1 to sec-1
      REAL*8, PARAMETER :: pmin2psec = 1.D0/60.D0
      INTEGER :: LATINDX,LANDINDX,IH,LEVTROP,I,J,L
      REAL*8, DIMENSION(16) :: HEIGHT
      REAL*8 :: ALTTROP,ALTTOP
      
      INTEGER :: J_1, J_0, J_0H, J_1H

      CALL GET(grid, J_STRT=J_0, J_STOP=J_1) 
c
      DO J=J_0,J_1
      DO I=1,IM
C****    Lightning source function altitude dependence:
C****    Determine if latitude is tropical:
         IF(J > JS.AND.J <= JN) THEN
           LATINDX=1
         ELSE
           LATINDX=2
         END IF
C****    Determine if over land or ocean:
         IF(fland(I,J) >= 0.5) THEN
           LANDINDX=1
         ELSE
           LANDINDX=2
         END IF
C****    Choose appropriate height file:
         DO IH=1,16
          HEIGHT(IH)=HGT(LATINDX,LANDINDX,IH)*0.01d0
         ENDDO
C****    Store local tropopause
         LEVTROP=LTROPO(I,J)
         ALTTROP=PHI(I,J,LEVTROP)*BYGRAV*1.d-3
         IF(ALTTROP == 0.) GOTO 1510
C****    Zero source accumulator
         SRCLIGHT = 0.d0 ! 1 to LM levels
C        Determine which GCM level the height belongs to
         DO 1505 IH=1,16
           L=1
 1502      continue
           ALTTOP=
     &     (PHI(I,J,L)+(PHI(I,J,L+1)-PHI(I,J,L))*0.5d0)*BYGRAV*1.d-3
c          RNOx_lgt is gN/min added per grid box (lightning NOx):
           IF(IH <= ALTTOP)THEN
             SRCLIGHT(L)=SRCLIGHT(L)+RNOx_lgt(I,J)*HEIGHT(IH)
             goto 1505
           elseif(IH-1 < ALTTOP)then
             SRCLIGHT(L)=SRCLIGHT(L)+RNOx_lgt(I,J)*HEIGHT(IH)*
     &       (ALTTOP-IH+1.d0)
             SRCLIGHT(L+1)=SRCLIGHT(L+1)+RNOx_lgt(I,J)*HEIGHT(IH)*
     &       (IH-ALTTOP)
             goto 1505
           else
             L=L+1
           ENDIF
           goto 1502
 1505    continue
 1510    continue    
C        save tracer 3D source. convert from gN/min to kgN/s :
         DO L=1,LEVTROP
           tr3Dsource(I,J,L,nOther,n_NOx) = 
     &     SRCLIGHT(L)*pmin2psec*1.d-3
         END DO
         DO L=LEVTROP+1,LM
            tr3Dsource(I,J,LEVTROP,nOther,n_NOx) = tr3Dsource(I,J
     $       ,LEVTROP,nOther,n_NOx) + SRCLIGHT(L)*pmin2psec*1.d-3
         END DO
      END DO ! I
      END DO ! J
         
C
      END SUBROUTINE get_lightning_NOx
c
c
      SUBROUTINE get_aircraft_NOx
!@sum  get_aircraft_NOx to define the 3D source of NOx from aircraft
!@auth Drew Shindell? / Greg Faluvegi / Jean Learner
!@ver  1.0 (based on DB396Tds3M23)
c
C**** GLOBAL parameters and variables:
C
      USE MODEL_COM, only: itime,jday,JDperY,im,jm,lm,ptop,psf,sig
      USE DOMAIN_DECOMP, only: GRID, GET, write_parallel
      USE CONSTANT, only: sday,hrday
      USE FILEMANAGER, only: openunit,closeunit
      USE FLUXES, only: tr3Dsource
      USE GEOM, only: dxyp
      USE TRACER_COM, only: itime_tr0,trname,n_NOx,nAircraft
      use TRACER_SOURCES, only: Laircr
C
      IMPLICIT NONE
c
!@var nanns,nmons: number of annual and monthly input files
!@var l,ll dummy loop variable
!@var pres local pressure variable
      integer, parameter :: nanns=0,nmons=1
      integer, dimension(nmons) :: mon_units, imon
      integer :: jdlast=0
      integer l,i,j,iu,k,ll,nt
      character*80 :: title
      character*12, dimension(nmons) :: mon_files=(/'NOx_AIRCRAFT'/)
      character(len=300) :: out_line
      logical :: LINJECT
      logical :: ifirst=.true.
      logical, dimension(nmons) :: mon_bins=(/.true./) ! binary file?
      real*8, allocatable, dimension(:,:,:,:) :: tlca, tlcb ! for monthly sources
!     real*8 tlca(im,jm,Laircr,nmons),tlcb(im,jm,Laircr,nmons)      
      real*8 frac, bySperHr      
      REAL*8, DIMENSION(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,Laircr,1)
     &                                     :: src
      REAL*8, DIMENSION(LM)                :: pres
      REAL*4, PARAMETER, DIMENSION(Laircr) :: PAIRL =
     & (/1013.25,898.74,794.95,701.08,616.40,540.19,471.81,
     &    410.60,355.99,307.42,264.36,226.32,193.30,165.10,
     &    141.01,120.44,102.87,87.866,75.048/)
      save ifirst,jdlast,tlca,tlcb,mon_units,imon
      INTEGER :: J_1, J_0, J_0H, J_1H

      CALL GET(grid, J_STRT=J_0, J_STOP=J_1) 
c
C**** Local parameters and variables and arguments:
c
C**** Aircraft NOx source input is monthly, on 19 levels. Therefore:
C     19x12=228 records. Read it in here and interpolated each day.

      if (itime < itime_tr0(n_NOx)) return

C****
C**** Monthly sources are interpolated to the current day
C**** The titles of the input files say this is in KG/m2/s, so no
C**** conversion is necessary:
C****
      if (ifirst) then
        call GET(grid, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)
        allocate(tlca(im,j_0H:j_1H,Laircr,nmons),
     &           tlcb(im,j_0H:j_1H,Laircr,nmons))
        ifirst=.false.
      endif  
      k = 1
      call openunit(mon_files(k),mon_units(k),mon_bins(k))
      call read_monthly_3Dsources(Laircr,mon_units(k),jdlast,
     *    tlca(:,:,:,k),tlcb(:,:,:,k),src(:,:,:,k),frac,imon(k))
      call closeunit(mon_units(k))
      jdlast = jday
      write(out_line,*)
     &'NOx from aircraft interpolated to current day',frac
      call write_parallel(trim(out_line)) 
C====
C====   Place aircraft sources onto model levels:
C====
      tr3Dsource(:,J_0:J_1,:,nAircraft,n_NOx) = 0.d0
      PRES(:)=SIG(:)*(PSF-PTOP)+PTOP
      DO J=J_0,J_1
       DO I=1,IM
        tr3Dsource(i,j,1,nAircraft,n_NOx) = SRC(I,J,1,1)*dxyp(j)
        DO LL=2,Laircr
          LINJECT=.TRUE.
          DO L=1,LM
           IF(PAIRL(LL) > PRES(L).AND.LINJECT) THEN
             tr3Dsource(i,j,l,nAircraft,n_NOx) =
     &       tr3Dsource(i,j,l,nAircraft,n_NOx) + SRC(I,J,LL,1)*dxyp(j)
             LINJECT=.FALSE.
           ENDIF
          ENDDO ! L
        ENDDO   ! LL
       END DO   ! I
      END DO    ! J
C
      return
      END SUBROUTINE get_aircraft_NOx
C
C
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
!@ver  1.0 
c
C**** GLOBAL parameters and variables:
C
      USE MODEL_COM, only: jday,im,jm,lm,ptop,psf,sig
      USE DOMAIN_DECOMP, only: GRID, GET, write_parallel
      USE FILEMANAGER, only: openunit,closeunit
      use TRACER_SOURCES, only: Lsulf
C
      IMPLICIT NONE

C**** Local parameters and variables and arguments:
c
!@var nmons: number of monthly input files
!@param Psulf pressure levels of the input file
      real*8, parameter, dimension(Lsulf) :: Psulf = (/
     & 0.9720D+03,0.9445D+03,0.9065D+03,
     & 0.8515D+03,0.7645D+03,0.6400D+03,0.4975D+03,0.3695D+03,
     & 0.2795D+03,0.2185D+03,0.1710D+03,0.1335D+03,0.1016D+03,
     & 0.7120D+02,0.4390D+02,0.2470D+02,0.1390D+02,0.7315D+01,
     & 0.3045D+01,0.9605D+00,0.3030D+00,0.8810D-01,0.1663D-01/)
      integer, parameter :: ncalls=3
      integer, dimension(ncalls):: mon_units, imon
      integer i,j,iu,k,l,nc
      integer, dimension(ncalls)  :: jdlast=(/0,0,0/)  
      character*80 :: title
      CHARACTER(LEN=*),  INTENT(IN) :: fn
      character(len=300) :: out_line
      logical, dimension(ncalls) :: mon_bins=(/.true.,.true.,.true./)
CCnoCClogical, dimension(ncalls) :: ifirst=(/.true.,.true.,.true./)
      logical :: ifirst=.true.     
      real*8, allocatable, dimension(:,:,:,:) :: tlca, tlcb
      REAL*8, DIMENSION(LM)    :: pres,srcLout
      REAL*8, DIMENSION(Lsulf) :: srcLin
      REAL*8, 
     &DIMENSION(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,Lsulf,ncalls)
     &                                                         :: src
      real*8, dimension(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM) ::
     &      field
      real*8 :: frac
      save jdlast,tlca,tlcb,mon_units,imon,ifirst
      INTEGER :: J_1, J_0, J_0H, J_1H

      CALL GET(grid, J_STRT=J_0, J_STOP=J_1)

      select case(trim(fn))
      case('DMS_FIELD') ; nc=1
      case('SO2_FIELD') ; nc=2
      case('SULFATE_SA'); nc=3
      case default
        call stop_model('please address filename in read_aero',255)
      end select
c
C**** Input file is monthly, on LM levels.
C     Read it in here and interpolated each day.
C     Interpolation in the vertical.
C
      if (ifirst) then
        call GET(grid, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)
        allocate(tlca(im,j_0H:j_1H,Lsulf,ncalls),
     &           tlcb(im,j_0H:j_1H,Lsulf,ncalls))
        ifirst     = .false.
      endif       
      call openunit(trim(fn),mon_units(nc),mon_bins(nc))
      call read_monthly_3Dsources(Lsulf,mon_units(nc),jdlast(nc),
     *   tlca(:,:,:,nc),tlcb(:,:,:,nc),src(:,:,:,nc),frac,imon(nc))
      call closeunit(mon_units(nc))
      jdlast(nc) = jday
      write(out_line,*)
     &trim(fn)//' interpolated to current day',frac
      call write_parallel(trim(out_line)) 
C====
C====   Place field onto model levels
C====              
      PRES(:)=SIG(:)*(PSF-PTOP)+PTOP
      DO J=J_0,J_1; DO I=1,IM
        srcLin(1:Lsulf)=src(I,J,1:Lsulf,nc)
        call LOGPINT(Lsulf,Psulf,srcLin,LM,PRES,srcLout,.true.)
        field(I,J,1:LM)=srcLout(1:LM)
      END DO   ; END DO    
C
      return
      END SUBROUTINE read_aero
C
C
      SUBROUTINE read_monthly_3Dsources(Ldim,iu,jdlast,tlca,tlcb,data1,
     & frac,imon)
!@sum Read in monthly sources and interpolate to current day
!@+   Calling routine must have the lines:
!@+      real*8 tlca(im,jm,Ldim,nm),tlcb(im,jm,Ldim,nm)
!@+      integer imon(nm)   ! nm=number of files that will be read
!@+      data jdlast /0/
!@+      save jdlast,tlca,tlcb,imon
!@+   Input: iu, the fileUnit#; jdlast
!@+   Output: interpolated data array + two monthly data arrays
!@auth Jean Lerner and others / Greg Faluvegi
      USE MODEL_COM, only: jday,im,jm,idofm=>JDmidOfM
      USE FILEMANAGER, only : NAMEUNIT
      USE DOMAIN_DECOMP, only : GRID,GET,READT_PARALLEL,REWIND_PARALLEL
     & ,write_parallel
      implicit none
!@var Ldim how many vertical levels in the read-in file?
!@var L dummy vertical loop variable
      integer Ldim,L,imon,iu,jdlast
      character(len=300) :: out_line
      real*8 :: frac
      real*8, DIMENSION(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO) ::
     &     A2D,B2D,dummy
      real*8, DIMENSION(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,Ldim) ::
     &     tlca,tlcb,data1
     
      integer :: J_0, J_1

      CALL GET(grid, J_STRT=J_0, J_STOP=J_1)     
C
      if (jdlast == 0) then   ! NEED TO READ IN FIRST MONTH OF DATA
        imon=1                ! imon=January
        if (jday <= 16)  then ! JDAY in Jan 1-15, first month is Dec
          do L=1,Ldim*11
            CALL READT_PARALLEL(grid,iu,NAMEUNIT(iu),0,dummy,1)
          end do
          DO L=1,Ldim
            CALL READT_PARALLEL(grid,iu,NAMEUNIT(iu),0,A2D,1)
            tlca(:,J_0:J_1,L)=A2D(:,J_0:J_1)
          END DO
          CALL REWIND_PARALLEL( iu )
        else              ! JDAY is in Jan 16 to Dec 16, get first month
  120     imon=imon+1
          if (jday > idofm(imon) .AND. imon <= 12) go to 120
          do L=1,Ldim*(imon-2)
            CALL READT_PARALLEL(grid,iu,NAMEUNIT(iu),0,dummy,1)
          end do
          DO L=1,Ldim
            CALL READT_PARALLEL(grid,iu,NAMEUNIT(iu),0,A2D,1)
            tlca(:,J_0:J_1,L)=A2D(:,J_0:J_1)
          END DO
          if (imon == 13)  CALL REWIND_PARALLEL( iu )
        end if
      else                         ! Do we need to read in second month?
        if (jday /= jdlast+1) then ! Check that data is read in daily
          if (jday /= 1 .OR. jdlast /= 365) then
            write(out_line,*)'Bad day values in read_monthly_3Dsources'
     &      //': JDAY,JDLAST=',JDAY,JDLAST
            call write_parallel(trim(out_line),crit=.true.) 
            call stop_model('Bad values in read_monthly_3Dsources',255)
          end if
          imon=imon-12             ! New year
          go to 130
        end if
        if (jday <= idofm(imon)) go to 130
        imon=imon+1                ! read in new month of data
        if (imon == 13) then
          CALL REWIND_PARALLEL( iu  )
        else  
          do L=1,Ldim*(imon-1)
            CALL READT_PARALLEL(grid,iu,NAMEUNIT(iu),0,dummy,1)
          end do
        endif
      end if
      DO L=1,Ldim
        CALL READT_PARALLEL(grid,iu,NAMEUNIT(iu),0,B2D,1)
        tlcb(:,J_0:J_1,L)=B2D(:,J_0:J_1)
      END DO
  130 continue
c**** Interpolate two months of data to current day
      frac = float(idofm(imon)-jday)/(idofm(imon)-idofm(imon-1))
      data1(:,J_0:J_1,:) =
     & tlca(:,J_0:J_1,:)*frac + tlcb(:,J_0:J_1,:)*(1.-frac)
      return
      end subroutine read_monthly_3Dsources


      subroutine get_CH4_IC(icall)
!@sum get_CH4_IC to generate initial conditions for methane.
!@auth Greg Faluvegi/Drew Shindell
!@ver  1.0 (based on DB396Tds3M23)
c
C**** GLOBAL parameters and variables:
C
      USE MODEL_COM, only  : im,jm,lm,ls1,JEQ,DTsrc
      USE DOMAIN_DECOMP, only : GRID,GET, write_parallel
      USE GEOM, only       : dxyp
      USE DYNAMICS, only   : am
      USE CONSTANT, only: mair
      USE TRACER_COM, only : trm, TR_MM, n_CH4, nStratwrite
      USE FLUXES, only: tr3Dsource
      USE TRCHEM_Shindell_COM, only: CH4altT, CH4altX, ch4_init_sh,
     *     ch4_init_nh,fix_CH4_chemistry,pfix_CH4_S,pfix_CH4_N
c
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
c
!@var CH4INIT temp variable for ch4 initial conditions
!@var I,J,L dummy loop variables
!@param bymair 1/molecular wt. of air = 1/mair
!@var JN J around 30 N
!@var JS J arount 30 S
!@var icall =1 (during run) =0 (first time)
      INTEGER, PARAMETER:: JS = JM/3 + 1, JN = 2*JM/3
      REAL*8, PARAMETER :: bymair = 1.d0/mair
      REAL*8 CH4INIT,bydtsrc
      INTEGER I, J, L, icall
      integer :: J_0, J_1

      CALL GET(grid, J_STRT=J_0, J_STOP=J_1)   
      
      bydtsrc=1.d0/DTsrc
C     First, the troposphere:
      DO L=1,LS1-1
      DO J=J_0,J_1
C       Initial latitudinal gradient for CH4:
        IF(J < JEQ)THEN ! Southern Hemisphere
          select case(fix_CH4_chemistry)
          case default
            CH4INIT=ch4_init_sh*TR_MM(n_CH4)*bymair*1.d-6*DXYP(J)
          case(1)
            CH4INIT=pfix_CH4_S*TR_MM(n_CH4)*bymair*DXYP(J)
          end select
        ELSE             ! Northern Hemisphere
          select case(fix_CH4_chemistry)
          case default
            CH4INIT=ch4_init_nh*TR_MM(n_CH4)*bymair*1.d-6*DXYP(J)
          case(1)
            CH4INIT=pfix_CH4_N*TR_MM(n_CH4)*bymair*DXYP(J)
          end select
        ENDIF
        select case(icall)
        case(0) ! initial conditions
          trm(:,j,l,n_CH4)=am(L,:,J)*CH4INIT
        case(1) ! overwriting
          tr3Dsource(:,j,l,nStratwrite,n_CH4) = 
     &    (am(L,:,J)*CH4INIT-trm(:,j,l,n_CH4))*bydtsrc
        end select
      END DO
      END DO
c
C     Now, the stratosphere:
      do L=LS1,LM
      do j=J_0,J_1
c
c     Define stratospheric ch4 based on HALOE obs for tropics
c     and extratropics and scale by the ratio of initial troposphere
c     mixing ratios to 1.79 (observed):
        IF(J < JEQ)THEN ! Southern Hemisphere
          select case(fix_CH4_chemistry)
          case default
            CH4INIT=
     &      ch4_init_sh/1.79d0*TR_MM(n_CH4)*bymair*1.E-6*DXYP(J)
          case(1)
            CH4INIT=pfix_CH4_S/1.79d0*TR_MM(n_CH4)*bymair*DXYP(J)
          end select
        ELSE             ! Northern Hemisphere
          select case(fix_CH4_chemistry)
          case default
            CH4INIT=
     &      ch4_init_nh/1.79d0*TR_MM(n_CH4)*bymair*1.E-6*DXYP(J)
          case(1)
            CH4INIT=pfix_CH4_N/1.79d0*TR_MM(n_CH4)*bymair*DXYP(J)
          end select
        ENDIF
        IF((J <= JS).OR.(J > JN)) THEN                 ! extratropics
          CH4INIT=CH4INIT*CH4altX(L)
        ELSE IF((J > JS).AND.(J <= JN)) THEN           ! tropics
          CH4INIT=CH4INIT*CH4altT(L)
        END IF
        select case(icall)
        case(0) ! initial conditions
          trm(:,j,l,n_CH4)= am(L,:,J)*CH4INIT
        case(1) ! overwriting
          tr3Dsource(:,j,l,nStratwrite,n_CH4) =
     &    (am(L,:,J)*CH4INIT-trm(:,j,l,n_CH4))*bydtsrc
        end select
      end do   ! j
      end do   ! l
C
      RETURN
      end subroutine get_CH4_IC
C  
C     
      SUBROUTINE calc_lightning(I,J,LMAX,LFRZ)
!@sum calc_lightning to calculate lightning flash amount and cloud-
!@+   to-ground amount, based on cloud top height. WARNING: this 
!@+   routine is apparently resolution dependant.  See the comments.
!@auth Colin Price (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on CB436Tds3M23)
C
C**** GLOBAL parameters and variables:
C
      USE LIGHTNING, only : JN,JS,RNOx_lgt
      USE MODEL_COM, only : fland
      USE GEOM,      only : bydxyp
      USE CONSTANT,  only : bygrav
      USE DYNAMICS,  only : gz
      USE TRDIAG_COM, only : ijs_CtoG,ijs_flash,taijs=>taijs_loc
c
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
C
!@var LMAX highest layer of current convective cloud event
!@var lmax_temp local copy of LMAX (alterable)
!@var LFRZ freezing level
!@var HTCON height of convection?
!@var HTFRZ height of freezing?
!@var lightning flashes per minute
!@var th,th2,th3,th4 thickness of cold sector, squared, cubed, etc.
!@var CG fraction of lightning that is cloud-to-ground
!@var zlt ?
!@var I model longitude index
!@var J model latitude index
!@param tune_land multiplier of flash rate over land to match observ.
!@param tune_ocean multiplier of flash rate over land to match observ.
!@param tune_NOx   multiplier of NOx production rate from lightning
      INTEGER, INTENT(IN) :: LMAX,LFRZ,I,J
      INTEGER lmax_temp
      REAL*8 HTCON,HTFRZ,flash,th,th2,th3,th4,zlt,CG
      REAL*8,PARAMETER::tune_land=2.2d0,tune_ocean=3.9d0
#ifdef SHINDELL_STRAT_CHEM      
     & ,tune_NOx=0.536d0
#else
     & ,tune_NOx=0.670d0
#endif
C
c The folowing simple algorithm calculates the lightning
c frequency in each gridbox using the moist convective cloud
c scheme.  The algorithm is based on the Price and Rind (1992)
c JGR paper, and differentiates between oceanic and continental
c convective clouds.  Only the non-entraining plume is used for
c the lightning.  However, I [Greg] have removed the IC variable
C from this routine because Gavin says: "IC distinguishes between
C entraining and non-entraining, but since non-entraining always
C go higher, those are the ones represented by LMCMAX."
C
c The lightning calculation uses the maximum height of the cloud,
c or the maximum depth of the mass flux penetration, LMAX, at
c every timestep.  In each timestep there may be a number of
c different values of LMAX, associated with clouds originating
c from different levels.  However, there are always less values
c of LMAX then the number of vertical layers. gz is the
c geopotential height, used from DYNAMICS module.

      lmax_temp=lmax
      if(j < JS.or.j > JN)lmax_temp=lmax+1 !30 S to 30 N
      HTCON=gz(i,j,lmax_temp)*bygrav*1.d-3
      HTFRZ=gz(i,j,LFRZ)*bygrav*1.d-3

c IF the gridbox is over land or coastline use the continental
c parameterization.  The factor 2.17 is derived from Fig. 1
c of the Price and Rind (1994) Mon. Wea. Rev. paper, and is
c related to the horizontal resolution of the GCM (in this case
c 4x5 degrees).  We use a threshold of 1% for continental
c gridboxes. The units are flashes per minute.

      If (fland(i,j) <= 0.01) then
        flash=tune_ocean*2.17d0*6.2d-4*(htcon**1.73d0)  ! ocean
      else
        flash= tune_land*2.17d0*3.44d-5*(htcon**4.92d0) ! continent
      end if

c The formulation by Price and Rind (1993) in Geophysical Res.
c Letters is used to calculate the fraction of total lightning
c that is cloud-to-ground lightning.  This uses the thickness of
c cold sector of the cloud (Hmax-Hzero) as the determining
c parameter.  The algorithm is only valid for thicknesses greater
c than 5.5km and less than 14km.

      th=(HTCON-HTFRZ)
      th=min(max(th,5.5d0),14.d0)
      th2=th*th
      th3=th2*th
      th4=th3*th
      zlt = 0.021d0*th4 - 0.648d0*th3 + 7.493d0*th2 
     &    - 36.544d0*th + 63.088d0
      CG=flash/(1.+zlt)

c Given the number of cloud-to-ground flashes, we can now calculate
c the NOx production rate based on the Price et al. (1997) JGR paper.
c The units are grams of Nitrogen per minute:

      RNOx_lgt(i,j)=15611.d0*(CG + 0.1d0*(flash-CG))*tune_NOx

C If flash is indeed in flashes/min, accumulate it in flashes:
       TAIJS(I,J,ijs_flash)=TAIJS(I,J,ijs_flash) + flash*60.d0
       TAIJS(I,J,ijs_CtoG) =TAIJS(I,J,ijs_CtoG)  +    CG*60.d0

      END SUBROUTINE calc_lightning
c
c
      SUBROUTINE get_sza(I,J,tempsza)
!@sum get_sza calculates the solar angle.  The intention is
!@+   that this routine will only be used when the COSZ1 from the 
!@+   radiation code is < 0, i.e. SZA > 90 deg.
!@auth Greg Faluvegi, based on the Harvard CTM routine SCALERAD
!@ver 1.0

c
C**** GLOBAL parameters and variables:  
C
      USE MODEL_COM, only: IM,JM,nday,Itime,DTsrc 
      USE CONSTANT, only: PI, radian
      USE TRCHEM_Shindell_COM, only: byradian
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments
C
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
      REAL*8, INTENT(OUT) :: tempsza
      INTEGER, INTENT(IN) :: I,J
      INTEGER DX,DY

      DX   = NINT(360./REAL(IM))
      DY   = NINT(180./REAL(JM))       
      vlat = -90. + REAL((J-1)*DY)
      if (J == 1)  vlat= -90. + 0.5*REAL(DY)
      if (j == JM) vlat=  90. - 0.5*REAL(DY)
      VLON = 180. - REAL((I-1)*DX)
C     This added 0.5 is to make the instantaneous zenith angle
C     more representative throughout the 1 hour time step:
!old  TIMEC = ((JDAY*24.0) + JHOUR  + 0.5)*3600.
      TIMEC = (mod(itime,365*nday) + nday + 0.5)*DTsrc  
!     if(DTsrc /= 3600.)call stop_model
!    &('DTsrc has changed. Please check on get_sza subroutine',255)
      P1 = 15.*(TIMEC/3600. - VLON/15. - 12.)
      FACT = (TIMEC/86400. - 81.1875)*ANG1
      P2 = 23.5*SIN(FACT*radian)
      P3 = VLAT
      temp = (SIN(P3*radian)*SIN(P2*radian)) +
     &(COS(P1*radian)*COS(P2*radian)*COS(P3*radian))
      tempsza = acos(temp)*byradian

      RETURN
      END SUBROUTINE get_sza
#endif
C
C
      SUBROUTINE LOGPINT(LIN,PIN,AIN,LOUT,POUT,AOUT,min_zero)
!@sum LOGPINT does vertical interpolation of column variable,
!@+   linearly in ln(P).
!@auth Greg Faluvegi
!@ver 1.0
C
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments
C
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
C
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
C
C
#ifdef TRACERS_ON
      SUBROUTINE special_layers_init
!@sum special_layers_init determined special altitude levels that
!@+ are important to Drew Shindell's tracer code.  This is done to
!@+ so that hardcoded levels are not needed when switching
!@+ vertical resolutions.
!@auth Greg Faluvegi
!@ver 1.0
C
C**** Global variables:
c
      USE MODEL_COM, only : LM,SIG,PSF,PTOP
      USE TRCHEM_Shindell_COM, only: L75P,L75M,F75P,F75M, ! FACT1
     &                           L569P,L569M,F569P,F569M  ! CH4
C
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments
C
!@var natural log of nominal pressure for verticle interpolations
!@var log75 natural log of 75 hPa
!@var log569 natural log of 569 hPa
      REAL*8 log75,log569
      REAL*8, DIMENSION(LM) :: LOGP
      INTEGER L
c      
      LOGP(:)=LOG(SIG(:)*(PSF-PTOP)+PTOP)
      log75=LOG(75.d0)
      log569=LOG(569.d0)
c 
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
c      
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
!@ver 1.0
C
C**** Global variables:
c
      USE MODEL_COM, only: DTsrc,Itime,ItimeI
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
      
      if(Itime == ItimeI)first_mod(:,:,m)=1 ! Ititialize
            
      if(m > nra_ch4.or.m < 1)call stop_model('nra_ch4 problem',255)
      if(iH(I,J,m) < 0.or.iH(I,J,m) > maxHR_ch4) then
        write(6,*) 'IJM iH maxHR_ch4=',I,J,m,iH(I,J,m),maxHR_ch4
        call stop_model('iH or maxHR_ch4 problem',255)
      endif
      nmax=NINT(24.d0*nicall*3600.d0/DTsrc)
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

#ifdef GFED_3D_BIOMASS
      SUBROUTINE get_GFED_biomass_burning(nt)
!@sum  get_GFED_biomass_burning to read the
!@+    3D source of various tracers from biomass burning
!@auth Greg Faluvegi / Jean Learner
!@ver  1.0 (based on DB396Tds3M23 aircraft reading)
c
C**** GLOBAL parameters and variables:
C
      USE MODEL_COM, only: itime,jday,JDperY,im,jm,lm
      USE DOMAIN_DECOMP, only: GRID, GET, write_parallel
      USE DYNAMICS, only: GZ
      USE CONSTANT, only: sday,hrday,byGRAV
      USE FILEMANAGER, only: openunit,closeunit
      USE GEOM, only: bydxyp
      USE TRACER_COM, only: itime_tr0,trname,nBiomass,ntm
      use TRACER_SOURCES, only: LbbGFED,GFED_BB
C
      IMPLICIT NONE
c
!@var pres local pressure variable
      integer :: mon_units,l,i,j,iu,k,ll,luse
      integer, dimension(ntm) :: jdlast, imon
      integer, intent(in) :: nt ! the tracer index
      character*80 :: title
      character*40 :: mon_files
      character(len=300) :: out_line
      logical :: ifirst=.true.
      logical :: mon_bins=.true. ! binary file?
      real*8, allocatable, dimension(:,:,:,:) :: tlca, tlcb ! for monthly sources
      real*8 frac, PIfact, bySperHr, a, b, c, ha, hb, hc
      REAL*8,DIMENSION(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LbbGFED,ntm)
     &                                     ::src
      REAL*4, PARAMETER, DIMENSION(LbbGFED) :: zmid =
     & (/ 50.  ,300., 750.,1500.,2500., 4500./)
      save ifirst,jdlast,tlca,tlcb,mon_units,imon
      INTEGER :: J_1, J_0, J_0H, J_1H, I_0, I_1

      CALL GET(grid, J_STRT=J_0, J_STOP=J_1 ,
     &               I_STRT=I_0, I_STOP=I_1 )

      if (itime < itime_tr0(nt)) return

      mon_files=trim(trname(nt))//'_IIASA_BBURN'
      bySperHr=1.d0/3600.d0

C****
C**** Monthly sources are interpolated to the current day
C**** Conversion is from KG/4x5grid/HR to KG/m2/s:
C****
      if (ifirst) then
        call GET(grid, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)
        allocate(tlca(im,j_0H:j_1H,LbbGFED,ntm),
     &           tlcb(im,j_0H:j_1H,LbbGFED,ntm))
        ifirst=.false.
      endif
      call openunit(mon_files,mon_units,mon_bins)
      call read_monthly_3Dsources(LbbGFED,mon_units,jdlast(nt),
     &tlca(:,:,:,nt),tlcb(:,:,:,nt),src(:,:,:,nt),frac,imon(nt))
      call closeunit(mon_units)
      jdlast(nt) = jday
      write(out_line,*)
     &trname(nt)//' from biomass burn interp to current day',frac
      call write_parallel(trim(out_line))

      do L=1,LbbGFED; do j=j_0,j_1; do i=1,IM
        GFED_BB(i,j,l,nt)=src(i,j,l,nt)*bySperHr
      enddo         ; enddo       ; enddo

      return
      END SUBROUTINE get_GFED_biomass_burning


      SUBROUTINE dist_GFED_biomass_burning(nt)
!@sum  dist_GFED_biomass_burning to locate on model levels the
!@+    3D source of various tracers from biomass burning
!@auth Greg Faluvegi / Jean Learner
!@ver  1.0 (based on DB396Tds3M23 aircraft reading)
c
C**** GLOBAL parameters and variables:
C
      USE MODEL_COM, only: im,jm,lm
      USE DOMAIN_DECOMP, only: GRID, GET, write_parallel
      USE DYNAMICS, only: GZ
      USE CONSTANT, only: byGRAV
      USE TRACER_COM, only: itime_tr0,trname,nBiomass,ntm,n_SO4
     & ,n_SO2
      use TRACER_SOURCES, only: LbbGFED,GFED_BB
      USE FLUXES, only : tr3Dsource
#ifdef TRACERS_SPECIAL_Shindell
      USE TRCHEM_Shindell_COM, only: PI_run,PIratio_bburn
#endif
C
      IMPLICIT NONE
c
!@var pres local pressure variable
      integer :: l,i,j,k,ll,luse
      integer, intent(in) :: nt ! the tracer index
      real*8 :: a, b, c, ha, hb, hc, fact_SO4
      REAL*4, PARAMETER, DIMENSION(LbbGFED) :: zmid =
     & (/ 50.  ,300., 750.,1500.,2500., 4500./)

      INTEGER :: J_1, J_0, J_0H, J_1H, I_0, I_1

      CALL GET(grid, J_STRT=J_0, J_STOP=J_1 ,
     &               I_STRT=I_0, I_STOP=I_1 )

      fact_SO4=0.0375d0

      do J=J_0,J_1 ; do I=I_0,I_1
        do LL=1,LbbGFED
          Luse=0
          loop_L: do L=1,LM
            if(L==LM) call stop_model('biomass burning crazy-high',255)
            if(l > 1)then
              ha=gz(i,j,l-1)*bygrav
              a=abs(ha-zmid(LL))
            else
              a=0.d0
            endif
            hb=gz(i,j,l)*bygrav
            b=abs(hb-zmid(LL))
            hc=gz(i,j,l+1)*bygrav
            c=abs(hc-zmid(LL))
            if((L==1.and.b<=c).or.(L /=1.and.(b<a.and.c>=b)))then
              Luse=L
              exit loop_L
            endif
          enddo loop_L
          if(nt==n_SO4)then
            tr3Dsource(i,j,Luse,nBiomass,nt) =
     &      tr3Dsource(i,j,Luse,nBiomass,nt) +
     &      GFED_BB(i,j,LL,n_SO2)*fact_SO4
          else
            tr3Dsource(i,j,Luse,nBiomass,nt) =
     &      tr3Dsource(i,j,Luse,nBiomass,nt)+GFED_BB(i,j,LL,nt)
          endif
        enddo

      enddo; enddo

      return
      END SUBROUTINE dist_GFED_biomass_burning
#endif
