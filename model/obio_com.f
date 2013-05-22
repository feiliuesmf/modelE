#include "rundeck_opts.h"
      MODULE obio_com
!@sum  obio_com contains the parameters, arrays and definitions
!@+    necessary for the OceanBiology routines
!@auth NR

      USE obio_dim

#ifdef OBIO_ON_GARYocean
      USE OCEANRES, only : kdm=>lmo
#else
      USE hycom_dim_glob
      USE hycom_scalars, only: baclin
#endif

      implicit none

c --- dobio       activate Watson Gregg's ocean biology code
      logical dobio
      data dobio/.true./
c

      real, ALLOCATABLE, DIMENSION(:,:)    :: tzoo2d
      real, ALLOCATABLE, DIMENSION(:,:,:)  :: tfac3d,wshc3d
      real, ALLOCATABLE, DIMENSION(:,:,:)  :: Fescav3d
      real, ALLOCATABLE, DIMENSION(:,:,:,:):: rmuplsr3d,rikd3d
      real, ALLOCATABLE, DIMENSION(:,:,:,:):: acdom3d
      real, ALLOCATABLE, DIMENSION(:,:,:)  :: gcmax         !cocco max growth rate
      real, ALLOCATABLE, DIMENSION(:,:)    :: pCO2          !partial pressure of CO2
      real, ALLOCATABLE, DIMENSION(:,:)    :: pCO2_glob
      real, ALLOCATABLE, DIMENSION(:,:)    :: cexpij        !detritus term (Pg,C/yr)
      real, ALLOCATABLE, DIMENSION(:,:)    :: cexp_glob
      real, ALLOCATABLE, DIMENSION(:,:)    :: caexpij       !CaCO3 export term (Pg,C/yr)
      real, ALLOCATABLE, DIMENSION(:,:)    :: caexp_glob
      real, ALLOCATABLE, DIMENSION(:,:)    :: pp2tot_day    !net pp total per day
      real, ALLOCATABLE, DIMENSION(:,:)    :: pp2tot_day_glob !net pp total per day
      real, ALLOCATABLE, DIMENSION(:,:)    :: tot_chlo      !tot chlorophyl at surf. layer
      real, ALLOCATABLE, DIMENSION(:,:)    :: tot_chlo_glob !tot chlorophyl at surf. layer
      real, ALLOCATABLE, DIMENSION(:,:,:,:):: rhs_obio      !rhs matrix
      real, ALLOCATABLE, DIMENSION(:,:,:)  :: chng_by       !integr tendency for total C

#ifdef OBIO_RUNOFF

#ifdef NITR_RUNOFF
!      real, ALLOCATABLE, DIMENSION(:,:)    :: rnitrmflo_loc      ! riverine nitrate mass flow rate (kg/s)
      real, ALLOCATABLE, DIMENSION(:,:)    :: rnitrconc_loc      ! riverine nitrate concentration (kg/kg)
#endif
#ifdef DIC_RUNOFF
      real, ALLOCATABLE, DIMENSION(:,:)    :: rdicconc_loc       ! riverine dic concentration (kg/kg)
#endif
#ifdef DOC_RUNOFF
      real, ALLOCATABLE, DIMENSION(:,:)    :: rdocconc_loc       ! riverine doc concentration (kg/kg)
#endif 
#ifdef SILI_RUNOFF
      real, ALLOCATABLE, DIMENSION(:,:)    :: rsiliconc_loc      ! riverine silica concentration (kg/kg)
#endif
#ifdef IRON_RUNOFF
      real, ALLOCATABLE, DIMENSION(:,:)    :: rironconc_loc      ! riverine iron concentration (kg/kg)
#endif
#ifdef POC_RUNOFF
      real, ALLOCATABLE, DIMENSION(:,:)    :: rpocconc_loc       ! riverine poc concentration (kg/kg)
#endif
#ifdef ALK_RUNOFF
      real, ALLOCATABLE, DIMENSION(:,:)    :: ralkconc_loc       ! riverine alkalinity concentration (mol/kg)
#endif

#endif

#ifndef OBIO_ON_GARYocean   /* NOT for Russell ocean */
      real, ALLOCATABLE, DIMENSION(:,:,:,:) :: tracav, tracav_loc
      real, ALLOCATABLE, DIMENSION(:,:,:)   :: plevav, plevav_loc
      real, ALLOCATABLE, DIMENSION(:,:) :: pCO2av, pCO2av_loc
      real, ALLOCATABLE, DIMENSION(:,:) :: pp2tot_dayav
      real, ALLOCATABLE, DIMENSION(:,:) :: pp2tot_dayav_loc
      real, ALLOCATABLE, DIMENSION(:,:) :: cexpav, cexpav_loc
      real, ALLOCATABLE, DIMENSION(:,:) :: caexpav, caexpav_loc
      real, ALLOCATABLE, DIMENSION(:,:) :: ao_co2fluxav,ao_co2fluxav_loc
      real, ALLOCATABLE, DIMENSION(:,:) :: ao_co2flux_loc  !ao CO2 on the ocean grid ***NOT for GASEXCH runs****
      real, ALLOCATABLE, DIMENSION(:,:) :: ao_co2flux_glob
#endif
#ifdef OBIO_ON_GARYocean
      real, ALLOCATABLE, DIMENSION(:,:,:,:):: tracer_loc    !only for gary ocean
      real, ALLOCATABLE, DIMENSION(:,:,:,:):: tracer        !only for gary ocean

      integer nstep0


      !test point
!!    integer, parameter :: itest=16, jtest=45    !equatorial Pacific                  2deg ocean
!!    integer, parameter :: itest=32, jtest=20    !southern ocean; Pacific          
      integer, parameter :: itest=1,  jtest=90     !equator Pacific
#else
!     integer, parameter :: itest=(220,320) equator Atlant; (245,275) 0.6S;274.5E Nino3
!     integer, parameter :: itest=316, jtest=258    !257.5E;-50.7S
      integer, parameter :: itest=243, jtest=1      !equator,dateline
      real :: diag_counter
#endif


      integer, parameter :: EUZ_DEFINED=1


      real, parameter :: rlamz=1.0            !Ivlev constant
      real, parameter :: greff=0.25           !grazing efficiency     
      real, parameter :: drate=0.05/24.0      !phytoplankton death rate/hr
      real, parameter :: dratez1=0.1/24.0     !zooplankton death rate/hr
      real, parameter :: dratez2=0.5/24.0     !zooplankton death rate/hr
      real, parameter :: regen=0.25           !regeneration fraction

      real ::  obio_deltath,obio_deltat       !time steps in hours because all rates are in hrs
      real ::  co2mon(26,12)        !26 years 1979-2004, 12 months

      integer npst,npnd   !starting and ending array index for PAR
      data npst,npnd /3,17/

      real WtoQ(nlt)           !Watts/m2 to quanta/m2/s conversion

! reduced rank arrays for obio_model calculations
      integer ihra_ij

      real cexp, caexp
      real temp1d(kdm),dp1d(kdm),obio_P(kdm,ntyp+n_inert)
     .                 ,det(kdm,ndet),car(kdm,ncar),avgq1d(kdm)
     .                 ,gcmax1d(kdm),saln1d(kdm),p1d(kdm+1)
     .                 ,alk1d(kdm),flimit(kdm,nchl,5)

      real atmFe_ij,covice_ij
      integer inwst,inwnd,jnwst,jnwnd     !starting and ending indices 
                                          !for daylight 
                                          !in i and j directions
      real acdom(kdm,nlt)                 !absorption coefficient of CDOM
      real P_tend(kdm,ntyp+n_inert)       !bio tendency (dP/dt)

#ifdef TRACERS_Alkalinity
      real A_tend(kdm)
#endif
#ifdef OBIO_RUNOFF

#ifdef NITR_RUNOFF
      real rnitrconc_ij
!     .    , rnitrmflo_ij
#endif
#ifdef DIC_RUNOFF
      real rdicconc_ij
#endif
#ifdef DOC_RUNOFF
      real rdocconc_ij
#endif
#ifdef SILI_RUNOFF
      real rsiliconc_ij
#endif
#ifdef IRON_RUNOFF
      real rironconc_ij
#endif
#ifdef POC_RUNOFF
      real rpocconc_ij
#endif
#ifdef ALK_RUNOFF
      real ralkconc_ij
#endif
#endif

      real rmuplsr(kdm,nchl)                  !growth+resp 
      real D_tend(kdm,ndet)                   !detrtial tendency
      real obio_ws(kdm+1,nchl)                !phyto sinking rate
      real tfac(kdm)                          !phyto T-dependence
      real pnoice(kdm)                        !pct ice-free
      real wsdet(kdm+1,ndet)                  !detrital sinking rate
      real rikd(kdm,nchl)                     !photoadaption state
      real tzoo                               !herbivore T-dependence
      real Fescav(kdm)                        !iron scavenging rate

C if NCHL_DEFINED > 3
      real wshc(kdm)                          !cocco sinking rate
C endif

      real :: C_tend(kdm,ncar)                !carbon tendency
      real :: pCO2_ij                         !partial pressure of CO2
      real :: gro(kdm,nchl)                   !realized growth rate
      integer :: day_of_month, hour_of_day

      real :: rhs(kdm,ntrac,17)         !secord arg-refers to tracer definition 
                                        !we are not using n_inert (always ntrac-1)
                                        !second argument refers to process that
                                        !contributes to tendency 

      real :: pp2_1d(kdm,nchl)          !net primary production

#ifndef TRACERS_GASEXCH_ocean_CO2
#ifdef TRACERS_OceanBiology
      real*8 :: ao_co2flux
#endif
#endif
      integer kzc
      real*8 :: carb_old,iron_old    !prev timesetep total carbon inventory

#ifdef restoreIRON
!per year change dI/I, + for sink, - for source
      !!!real*8 :: Iron_BC = 0.002
      real*8 :: Iron_BC = -0.005
#endif

#ifdef OBIO_ON_GARYocean
      real*8, allocatable :: focean_glob(:,:)
      integer, allocatable :: lmom_glob(:,:)
#endif

      END MODULE obio_com

!------------------------------------------------------------------------------
      subroutine gather_pCO2
      USE obio_com
#ifdef OBIO_ON_GARYocean
      USE OCEANR_DIM, only : ogrid
#else
      USE HYCOM_DIM, only : ogrid
#endif
      USE DOMAIN_DECOMP_1D, ONLY: PACK_DATA
 
      call pack_data( ogrid, pCO2, pCO2_glob )
#ifndef OBIO_ON_GARYocean
      call pack_data( ogrid, pp2tot_day, pp2tot_day_glob )
      call pack_data( ogrid, cexpij, cexp_glob )
      call pack_data( ogrid, caexpij, caexp_glob )
#endif

      end subroutine gather_pCO2

      subroutine gather_chl
      USE obio_com
#ifdef OBIO_ON_GARYocean
      USE OCEANR_DIM, only : ogrid
#else
      USE HYCOM_DIM, only : ogrid
#endif
      USE DOMAIN_DECOMP_1D, ONLY: PACK_DATA

      call pack_data( ogrid, tot_chlo, tot_chlo_glob )

      end subroutine gather_chl

!------------------------------------------------------------------------------
      subroutine alloc_obio_com
      USE obio_com
      USE obio_dim

#ifdef OBIO_ON_GARYocean
      USE OCEANR_DIM, only : ogrid
      USE DOMAIN_DECOMP_1D, only : getDomainBounds
      USE OCEANRES, only :idm=>imo,jdm=>jmo,kdm=>lmo
#else
      USE hycom_dim_glob 
      USE hycom_dim, only : i_0h,i_1h,j_0h,j_1h
#endif

      implicit none

#ifdef OBIO_ON_GARYocean
c**** Extract domain decomposition info
      INTEGER :: j_0h,j_1h,i_0h,i_1h

      I_0H = ogrid%I_STRT_HALO
      I_1H = ogrid%I_STOP_HALO
      J_0H = ogrid%J_STRT_HALO
      J_1H = ogrid%J_STOP_HALO


      ALLOCATE(tracer_loc(i_0h:i_1h,j_0h:j_1h,kdm,ntrac))
      ALLOCATE(tracer(idm,jdm,kdm,ntrac))
#endif

      ALLOCATE(tzoo2d(i_0h:i_1h,j_0h:j_1h))
      ALLOCATE(wshc3d(i_0h:i_1h,j_0h:j_1h,kdm))
      ALLOCATE(Fescav3d(i_0h:i_1h,j_0h:j_1h,kdm))
      ALLOCATE(rmuplsr3d(i_0h:i_1h,j_0h:j_1h,kdm,nchl),
     &            rikd3d(i_0h:i_1h,j_0h:j_1h,kdm,nchl))
      ALLOCATE(acdom3d(i_0h:i_1h,j_0h:j_1h,kdm,nlt))
      ALLOCATE(tfac3d(i_0h:i_1h,j_0h:j_1h,kdm))
      ALLOCATE(gcmax(i_0h:i_1h,j_0h:j_1h,kdm))
      ALLOCATE(pCO2(i_0h:i_1h,j_0h:j_1h))           
      ALLOCATE(pCO2_glob(idm,jdm))           
      ALLOCATE(pp2tot_day(i_0h:i_1h,j_0h:j_1h))
      ALLOCATE(pp2tot_day_glob(idm,jdm))
      ALLOCATE(tot_chlo(i_0h:i_1h,j_0h:j_1h))
      ALLOCATE(tot_chlo_glob(idm,jdm))
      ALLOCATE(rhs_obio(i_0h:i_1h,j_0h:j_1h,ntrac,17))
      ALLOCATE(chng_by(i_0h:i_1h,j_0h:j_1h,14))

#ifdef OBIO_RUNOFF
#ifdef NITR_RUNOFF
!      ALLOCATE(rnitrmflo_loc(i_0h:i_1h,j_0h:j_1h))
      ALLOCATE(rnitrconc_loc(i_0h:i_1h,j_0h:j_1h))
#endif
#ifdef DIC_RUNOFF
      ALLOCATE(rdicconc_loc(i_0h:i_1h,j_0h:j_1h))
#endif
#ifdef DOC_RUNOFF
      ALLOCATE(rdocconc_loc(i_0h:i_1h,j_0h:j_1h))
#endif
#ifdef SILI_RUNOFF
      ALLOCATE(rsiliconc_loc(i_0h:i_1h,j_0h:j_1h))
#endif
#ifdef IRON_RUNOFF
      ALLOCATE(rironconc_loc(i_0h:i_1h,j_0h:j_1h))
#endif
#ifdef POC_RUNOFF
      ALLOCATE(rpocconc_loc(i_0h:i_1h,j_0h:j_1h))
#endif
#ifdef ALK_RUNOFF
      ALLOCATE(ralkconc_loc(i_0h:i_1h,j_0h:j_1h))
#endif
#endif

#ifndef OBIO_ON_GARYocean   /* NOT for Russell ocean */
      ALLOCATE(ao_co2flux_loc(i_0h:i_1h,j_0h:j_1h))
      ALLOCATE(ao_co2flux_glob(idm,jdm))
      ALLOCATE(tracav(idm,jdm,kdm,ntrac))
      ALLOCATE(tracav_loc(i_0h:i_1h,j_0h:j_1h,kdm,ntrac))
      ALLOCATE(plevav(idm,jdm,kdm))
      ALLOCATE(plevav_loc(i_0h:i_1h,j_0h:j_1h,kdm))
      ALLOCATE(pCO2av(idm,jdm))
      ALLOCATE(pp2tot_dayav(idm,jdm))
      ALLOCATE(cexpij(i_0h:i_1h,j_0h:j_1h))
      ALLOCATE(cexp_glob(idm,jdm))           
      ALLOCATE(caexpij(i_0h:i_1h,j_0h:j_1h))
      ALLOCATE(caexp_glob(idm,jdm))           
      ALLOCATE(cexpav(idm,jdm))
      ALLOCATE(caexpav(idm,jdm))
      ALLOCATE(pCO2av_loc(i_0h:i_1h,j_0h:j_1h))
      ALLOCATE(pp2tot_dayav_loc(i_0h:i_1h,j_0h:j_1h))
      ALLOCATE(cexpav_loc(i_0h:i_1h,j_0h:j_1h))
      ALLOCATE(caexpav_loc(i_0h:i_1h,j_0h:j_1h))
      ALLOCATE(ao_co2fluxav(idm,jdm))
      ALLOCATE(ao_co2fluxav_loc(i_0h:i_1h,j_0h:j_1h))
#endif

      end subroutine alloc_obio_com

!------------------------------------------------------------------------------

#ifdef OBIO_ON_GARYocean

      subroutine obio_exports_init
      use obio_com, only : pCO2
      use ofluxes, only : ocnatm
      implicit none
      ocnatm%pCO2(:,:) = pCO2(:,:)
      end subroutine obio_exports_init

      subroutine def_rsf_obio(fid)
!@sum  def_rsf_ocean defines ocean array structure in restart files
!@auth M. Kelley
!@ver  beta
      USE OCEANR_DIM, only : grid=>ogrid
      Use OCN_TRACER_COM, Only : ntm,trname
      USE obio_forc, only : avgq,tirrq3d,ihra
      USE obio_com, only : gcmax,nstep0
     &     ,tracer=>tracer_loc,pCO2,pp2tot_day
      use pario, only : defvar
      use domain_decomp_1d, only : getDomainBounds
      implicit none
      integer fid   !@var fid file id
      integer :: n
      call defvar(grid,fid,nstep0,'obio_nstep0')
      call defvar(grid,fid,avgq,'avgq(dist_imo,dist_jmo,lmo)')
      call defvar(grid,fid,gcmax,'gcmax(dist_imo,dist_jmo,lmo)')
      call defvar(grid,fid,tirrq3d,'tirrq3d(dist_imo,dist_jmo,lmo)')
      call defvar(grid,fid,ihra,'ihra(dist_imo,dist_jmo)')
      do n=1,ntm
        call defvar(grid,fid,tracer(:,:,:,n),
     &       'obio_'//trim(trname(n))//'(dist_imo,dist_jmo,lmo)')
      enddo
      call defvar(grid,fid,pCO2,'pCO2(dist_imo,dist_jmo)')
      call defvar(grid,fid,pp2tot_day,'pp2tot_day(dist_imo,dist_jmo)')
      return
      end subroutine def_rsf_obio

      subroutine new_io_obio(fid,iaction)
!@sum  new_io_ocean read/write ocean arrays from/to restart files
!@auth M. Kelley
!@ver  beta new_ prefix avoids name clash with the default version
      use model_com, only : ioread,iowrite
      USE OCEANR_DIM, only : grid=>ogrid
      use pario, only : write_dist_data,read_dist_data,
     &     write_data,read_data
      Use OCN_TRACER_COM, Only : ntm,trname
      use model_com, only : nstep=>itime
      USE obio_forc, only : avgq,tirrq3d,ihra
      USE obio_com, only : gcmax,nstep0
     &     ,tracer=>tracer_loc,pCO2,pp2tot_day
      use domain_decomp_1d, only : getDomainBounds
      implicit none
      integer fid   !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
      integer :: n
      select case (iaction)
      case (iowrite)            ! output to restart file
        call write_data(grid,fid,'obio_nstep0',nstep)
        call write_dist_data(grid,fid,'avgq',avgq)
        call write_dist_data(grid,fid,'gcmax',gcmax)
        call write_dist_data(grid,fid,'tirrq3d',tirrq3d)
        call write_dist_data(grid,fid,'ihra',ihra)
        do n=1,ntm
          call write_dist_data(grid,fid,'obio_'//trim(trname(n)),
     &         tracer(:,:,:,n))
        enddo
        call write_dist_data(grid,fid,'pCO2',pCO2)
        call write_dist_data(grid,fid,'pp2tot_day',pp2tot_day)
      case (ioread)            ! input from restart file
        call read_data(grid,fid,'obio_nstep0',nstep0,
     &       bcast_all=.true.)
        call read_dist_data(grid,fid,'avgq',avgq)
        call read_dist_data(grid,fid,'gcmax',gcmax)
        call read_dist_data(grid,fid,'tirrq3d',tirrq3d)
        call read_dist_data(grid,fid,'ihra',ihra)
        do n=1,ntm
          call read_dist_data(grid,fid,'obio_'//trim(trname(n)),
     &         tracer(:,:,:,n))
        enddo
        call read_dist_data(grid,fid,'pCO2',pCO2)
        call read_dist_data(grid,fid,'pp2tot_day',pp2tot_day)
      end select
      return
      end subroutine new_io_obio

#else

      subroutine obio_set_data_after_archiv
      USE obio_com, only:
     .     diag_counter
     .    ,plevav_loc,ao_co2fluxav_loc,tracav_loc, pp2tot_dayav_loc
     .    ,pCO2av_loc, cexpav_loc
#ifdef TRACERS_Alkalinity
     .    ,caexpav_loc
#endif
      implicit none
      diag_counter=0
      tracav_loc = 0
      plevav_loc = 0
      ao_co2fluxav_loc=0
      pCO2av_loc = 0
      pp2tot_dayav_loc = 0
      cexpav_loc = 0
#ifdef TRACERS_Alkalinity
      caexpav_loc = 0
#endif
      end subroutine obio_set_data_after_archiv

      subroutine obio_gather_before_archive
      USE HYCOM_DIM, only : ogrid
      USE DOMAIN_DECOMP_1D, ONLY: PACK_DATA
#if (defined TRACERS_OceanBiology) && (defined TRACERS_GASEXCH_ocean_CO2)
      USE obio_com, only: tracav,tracav_loc,
     .    plevav,plevav_loc
#endif
#ifdef TRACERS_OceanBiology
      USE obio_com, only:
     .     ao_co2fluxav_loc, ao_co2fluxav
     .    ,pCO2av_loc, pCO2av
     .    ,pp2tot_dayav_loc, pp2tot_dayav
     .    ,cexpav_loc, cexpav
#ifdef TRACERS_Alkalinity
     .    ,caexpav_loc, caexpav
#endif
#endif
      implicit none 

#ifdef TRACERS_OceanBiology
      call gather_chl

      call pack_data(ogrid, ao_co2fluxav_loc, ao_co2fluxav)
      call pack_data(ogrid, pCO2av_loc, pCO2av)
      call pack_data(ogrid, pp2tot_dayav_loc, pp2tot_dayav)
      call pack_data(ogrid, cexpav_loc, cexpav)
#ifdef TRACERS_Alkalinity
      call pack_data(ogrid, caexpav_loc, caexpav)
#endif
#endif

#if (defined TRACERS_OceanBiology) || defined (TRACERS_GASEXCH_ocean) \
      || (defined TRACERS_AGE_OCEAN) || (defined TRACERS_OCEAN_WATER_MASSES) \
      || (defined TRACERS_ZEBRA)
      call pack_data( ogrid,  tracav_loc, tracav )
      call pack_data( ogrid,  plevav_loc, plevav )
#endif
      return
      end subroutine obio_gather_before_archive

      subroutine def_rsf_obio(fid)
!@sum  def_rsf_ocean defines ocean array structure in restart files
!@auth M. Kelley
!@ver  beta
      USE HYCOM_DIM, only : grid=>ogrid
      use pario, only : defvar
      USE HYCOM_ARRAYS, only : tracer
      USE obio_forc, only : avgq,tirrq3d,ihra
      USE obio_com, only : gcmax,pCO2av=>pCO2av_loc,pp2tot_day,
     &     ao_co2fluxav=>ao_co2fluxav_loc,
     &     pp2tot_dayav=>pp2tot_dayav_loc,
     &     cexpav=>cexpav_loc,
     &     diag_counter, tracav=>tracav_loc, plevav=>plevav_loc
      use obio_dim, only : trname
      implicit none
      integer fid   !@var fid file id
      integer :: n
      character(len=14) :: str2d
      character(len=18) :: str3d
      character(len=20) :: str3d2
      str2d ='(idm,dist_jdm)'
      str3d ='(idm,dist_jdm,kdm)'
      str3d2='(idm,dist_jdm,kdmx2)'

      do n=1,size(trname)
        call defvar(grid,fid,tracer(:,:,:,n),
     &       trim(trname(n))//str3d)
      enddo

c      call defvar(grid,fid,nstep,'obio_nstep0')
      call defvar(grid,fid,diag_counter,'obio_diag_counter')
      call defvar(grid,fid,avgq,'avgq'//str3d)
      call defvar(grid,fid,gcmax,'gcmax'//str3d)
      call defvar(grid,fid,tirrq3d,'tirrq3d'//str3d)
      call defvar(grid,fid,ihra,'ihra'//str2d)
      call defvar(grid,fid,pCO2av,'pCO2av'//str2d)
      call defvar(grid,fid,pp2tot_dayav,'pp2tot_dayav'//str2d)
      call defvar(grid,fid,ao_co2fluxav,'ao_co2fluxav'//str2d)
      call defvar(grid,fid,cexpav,'cexpav'//str2d)
      call defvar(grid,fid,pp2tot_day,'pp2tot_day'//str2d)
      call defvar(grid,fid,tracav,'tracav(idm,dist_jdm,kdm,ntrcr)')
      call defvar(grid,fid,plevav,'plevav'//str3d)

      return
      end subroutine def_rsf_obio

      subroutine new_io_obio(fid,iaction)
!@sum  new_io_ocean read/write ocean arrays from/to restart files
!@auth M. Kelley
!@ver  beta new_ prefix avoids name clash with the default version
      use model_com, only : ioread,iowrite
      use pario, only : write_dist_data,read_dist_data,
     &     write_data,read_data
      USE HYCOM_DIM, only : grid=>ogrid
      USE HYCOM_ARRAYS, only : tracer
      USE obio_forc, only : avgq,tirrq3d,ihra
      USE obio_com, only : gcmax,pCO2av=>pCO2av_loc,pp2tot_day,
     &     ao_co2fluxav=>ao_co2fluxav_loc,
     &     pp2tot_dayav=>pp2tot_dayav_loc,
     &     cexpav=>cexpav_loc,
     &     diag_counter, tracav=>tracav_loc, plevav=>plevav_loc
      use obio_dim, only : trname
      implicit none
      integer fid   !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
      integer :: n
      select case (iaction)
      case (iowrite)            ! output to restart file
        do n=1,size(trname)
          call write_dist_data(grid,fid,trim(trname(n)),tracer(:,:,:,n))
        enddo
        call write_data(grid,fid,'obio_diag_counter',diag_counter)
        call write_dist_data(grid,fid,'avgq',avgq)
        call write_dist_data(grid,fid,'gcmax',gcmax)
        call write_dist_data(grid,fid,'tirrq3d',tirrq3d)
        call write_dist_data(grid,fid,'ihra',ihra)
        call write_dist_data(grid,fid,'pCO2av',pCO2av)
        call write_dist_data(grid,fid,'pp2tot_dayav',pp2tot_dayav)
        call write_dist_data(grid,fid,'ao_co2fluxav',ao_co2fluxav)
        call write_dist_data(grid,fid,'cexpav',cexpav)
        call write_dist_data(grid,fid,'pp2tot_day',pp2tot_day)
        call write_dist_data(grid,fid,'tracav',tracav)
        call write_dist_data(grid,fid,'plevav',plevav)
      case (ioread)            ! input from restart file
        do n=1,size(trname)
          call read_dist_data(grid,fid,trim(trname(n)),tracer(:,:,:,n))
        enddo
        call read_data(grid,fid,'obio_diag_counter',diag_counter)
        call read_dist_data(grid,fid,'avgq',avgq)
        call read_dist_data(grid,fid,'gcmax',gcmax)
        call read_dist_data(grid,fid,'tirrq3d',tirrq3d)
        call read_dist_data(grid,fid,'ihra',ihra)
        call read_dist_data(grid,fid,'pCO2av',pCO2av)
        call read_dist_data(grid,fid,'pp2tot_dayav',pp2tot_dayav)
        call read_dist_data(grid,fid,'ao_co2fluxav',ao_co2fluxav)
        call read_dist_data(grid,fid,'cexpav',cexpav)
        call read_dist_data(grid,fid,'pp2tot_day',pp2tot_day)
        call read_dist_data(grid,fid,'tracav',tracav)
        call read_dist_data(grid,fid,'plevav',plevav)
      end select
      return
      end subroutine new_io_obio
#endif /* which ocean */
