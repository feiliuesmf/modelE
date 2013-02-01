#include "rundeck_opts.h"
!------------------------------------------------------------------------------
      subroutine setDefaultSpec(n, tracer)
!------------------------------------------------------------------------------
      use Dictionary_mod, only: sync_param
      use RunTimeControls_mod, only: tracers_amp
      use RunTimeControls_mod, only: tracers_tomas
      use OldTracer_mod, only: trName, do_fire
      use OldTracer_mod, only: nBBsources, set_nBBsources
      use DOMAIN_DECOMP_ATM, only: am_i_root
      use TRACER_COM, only: tracers
      use TRACER_COM, only: sect_name
      use TRACER_COM, only: set_ntsurfsrc, ntsurfsrc, ntsurfsrcmax
      use TRACER_COM, only: num_sectors
      use Tracer_mod, only: Tracer_type
      use Tracer_mod, only: setProperty
      use Tracer_mod, only: findSurfaceSources
      use Tracer_mod, only: addSurfaceSource
      use TracerBundle_mod, only: getTracer
#ifdef TRACERS_SPECIAL_Shindell      
      use TRCHEM_Shindell_COM, only: use_rad_ch4
#endif
      implicit none

      integer, intent(in) :: n
      type (Tracer_type), pointer :: tracer

      logical :: checkSourceName
      integer :: val

      call setProperty(tracer, 'ntSurfSrc', 0)

!     The following section will check for rundeck file of
!     the form: trname_01, trname_02... and thereby define
!     the ntsurfsrc(n). If those files exist it reads an
!     80 char header to get information including the
!     source name (ssame-->{sname,lname,etc.}. ntsurfsrc(n)
!     get set to zero if those files aren't found:
!     (I can enclose this in an ifdef if it causes problems
!     for people). num_srf_sources routine also assigns
!     sources to sectors, if desired:
!     general case:

      if (tracers_amp .or. tracers_tomas) then
         checkSourceName = .false.
      else
         checkSourceName = .true.
      end if

      call findSurfaceSources(tracer, checkSourceName, 
     &     sect_name(1:num_sectors))

#ifdef DYNAMIC_BIOMASS_BURNING
!     allow some tracers to have biomass burning based on fire model:
        select case (trname(n))
          case('NOx','CO','Alkenes','Paraffin','BCB','OCB','NH3','SO2'
     &         'vbsAm2', 'vbsAm1', 'vbsAz',  'vbsAp1', 'vbsAp2',
     &         'vbsAp3', 'vbsAp4', 'vbsAp5', 'vbsAp6',
#ifdef TRACERS_SPECIAL_Shindell
     &         ,'CH4'           ! in here to avoid potential Lerner tracers conflict
#endif
#ifdef TRACERS_TOMAS
     &         ,'AECOB_01','AOCOB_01' !BCB and OCB hygroscopities? Need to put emission into OB and IL.

#endif
     &         )
          call set_do_fire(n, .true.)
        end select
#endif /* DYNAMIC_BIOMASS_BURNING */

!     allow some tracers to have biomass burning sources that mix over
!     PBL layers (these become 3D sources no longer within ntsurfsrc(n)):
        select case (trname(n))
          case ('Alkenes', 'CO', 'NOx', 'Paraffin',
#ifdef TRACERS_SPECIAL_Shindell
     &         'CH4',           ! in here to avoid potential Lerner tracers conflict
#endif
     &         'AECOB_01','AOCOB_01', 
     &         'NH3', 'SO2', 'BCB', 'OCB', ! do not include sulfate here
     &         'vbsAm2', 'vbsAm1', 'vbsAz',  'vbsAp1', 'vbsAp2',
     &         'vbsAp3', 'vbsAp4', 'vbsAp5', 'vbsAp6',
     &         'M_BC1_BC', 'M_OCC_OC', 'M_BOC_BC', 'M_BOC_OC')
          val = nBBsources(n)
          call sync_param(trim(trname(n))//"_nBBsources",val)
          call set_nBBsources(n, val)
          if(nBBsources(n)>0)then
            if(do_fire(n))then
              if(am_i_root())write(6,*)
     &             'nBBsource>0 for ',trim(trname(n)),' do_fire=t'
              call stop_model('nBBsource do_fire conflict',13)
            else
              call set_ntsurfsrc(n, ntsurfsrc(n)-nBBsources(n))
            end if
          end if
        end select
        if(do_fire(n) .and.  (ntsurfsrc(n)+1 > ntsurfsrcmax))then
          write(6,*)trname(n),'ntsurfsrc+1 > max of ',ntsurfsrcmax
          call stop_model('do_fire+ntsurfsrc too large',13)
        end if
        if(ntsurfsrc(n)+nBBsources(n) > ntsurfsrcmax)then
          write(6,*)trname(n),'ntsurfsrc+nBBsources > max of ',
     &         ntsurfsrcmax
          call stop_model('ntsurfsrc+nBBsources too large',13)
        end if

!     other special cases:
#ifndef TRACERS_AEROSOLS_SOA
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP) ||\
        (defined TRACERS_TOMAS)
        select case (trname(n))
        case ('OCII', 'M_OCC_OC', 'SOAgas') ! this handles OCT_src (terpene source)
          tracer => getTracer(tracers, trname(n))
          call addSurfaceSource(tracer, "Terpene_source")
        end select
#endif
#endif  /* TRACERS_AEROSOLS_SOA */
#ifdef TRACERS_SPECIAL_Shindell
        if (trname(n)=='CH4' .and. use_rad_ch4/=0) then
          call set_ntsurfsrc(n,0)
        end if
#endif

      end subroutine setDefaultSpec

!------------------------------------------------------------------------------
      subroutine initTracerMetadata()
!------------------------------------------------------------------------------
#ifdef TRACERS_SPECIAL_Shindell
      use ShindellTracersMetadata_mod
#endif   
   
#ifdef TRACERS_TOMAS
      use TomasTracersMetadata_mod
#endif    

#ifdef TRACERS_AMP
      use AmpTracersMetadata_mod
#endif   
   
#ifdef TRACERS_AEROSOLS_Koch
      use KochTracersMetadata_mod
#endif   
   
#ifdef TRACERS_NITRATE
      use sharedTracersMetadata_mod, only:
     &  NH3_setSpec, NH4_setSpec
#endif
 
#ifdef TRACER_SPECIAL_Lerner
      use LernerTracersMetadata_mod
      USE TRACERS_MPchem_COM, only: n_MPtable,tcscale
!     @dbparam dsol describes portion of solar cycle being modeled for linoz
!     @+      +1.0 = solar max, 0.0 = neutral, -1.0 = solar min
      USE LINOZ_CHEM_COM, only: dsol
#endif

#ifdef TRACERS_RADON
      use sharedTracersMetadata_mod, only: Rn222_setSpec
#endif
   
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP) ||\
    (defined TRACERS_TOMAS) 
      use TRACER_COM, only: aer_int_yr
      USE TRACER_COM, only: offline_dms_ss, offline_ss
#endif

#ifdef TRACERS_ON
      USE TRDIAG_COM
#endif

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
      (defined TRACERS_QUARZHEM) || (defined TRACERS_TOMAS)
      use tracers_dust,only : imDust,prefDustSources,fracClayPDFscheme
     &     ,fracSiltPDFscheme
#endif

#ifdef TRACERS_AEROSOLS_SOA
      USE TRACERS_SOA, only: n_soa_i,n_soa_e,soa_init
#endif  /* TRACERS_AEROSOLS_SOA */

#ifdef TRACERS_AEROSOLS_VBS
      USE AEROSOL_SOURCES, only: VBSemifact
      USE TRACERS_VBS, only: vbs_tr,vbs_init
#endif  /* TRACERS_AEROSOLS_VBS */

#if (defined TRACERS_COSMO)
      USE COSMO_SOURCES, only: be7_src_param
#endif

#ifdef TRACERS_GASEXCH_ocean_CO2
      USE obio_forc, only : atmCO2
#endif

#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_TOMAS)
      USE AEROSOL_SOURCES, only: tune_ss1, tune_ss2, BBinc
#endif

#if (defined TRACERS_OCEAN) &&  !defined(TRACERS_OCEAN_INDEP)
! atmosphere copies atmosphere-declared tracer info to ocean
! so that the ocean can "inherit" it without referencing atm. code
      use ocn_tracer_com, only : 
     &     n_Water_ocn      => n_Water,
     &     ntrocn_ocn       => ntrocn,
     &     to_per_mil_ocn   => to_per_mil,
     &     t_qlimit_ocn     => t_qlimit,
     &     conc_from_fw_ocn => conc_from_fw,
     &     trdecay_ocn      => trdecay,
     &     trw0_ocn         => trw0
#endif
      USE DOMAIN_DECOMP_ATM, only: AM_I_ROOT
      USE CONSTANT, only: mair,mwat,pi,gasc
      USE RESOLUTION, only : jm,lm
      use TracerBundle_mod, only: getNumTracers
      USE TRACER_COM, only: ntm, ef_fact3d, no_emis_over_ice

      USE TRACER_COM, only: tracers
      use TRACER_COM, only: set_ntsurfsrc
      use TRACER_COM, only: 
     &     n_Air, n_water,
     &     n_H2O18,  n_HDO,   n_HTO,   
     &     n_Ox,       n_NOx,
     &     n_Pb210,n_Be7,   n_Be10,
     &     n_CFCn,   n_CO2n,  n_Age,
     &     n_SO4_d1,  n_SO4_d2,
     &     n_SO4_d3,n_N_d1,  n_N_d2,  n_N_d3,
     &     n_NH3,   n_NH4,   n_NO3p,
     &     n_vbsGm2, n_vbsGm1, n_vbsGz,  n_vbsGp1, n_vbsGp2,
     &     n_vbsGp3, n_vbsGp4, n_vbsGp5, n_vbsGp6,
     &     n_vbsAm2, n_vbsAm1, n_vbsAz,  n_vbsAp1, n_vbsAp2,
     &     n_vbsAp3, n_vbsAp4, n_vbsAp5, n_vbsAp6,
     &     n_OCocean,
     &     n_clay,   n_silt1, n_silt2, n_silt3, n_silt4

      USE Dictionary_mod
      use OldTracer_mod, only: trName, nGAS, nWater, nPart
      use OldTracer_mod, only: tr_wd_type, do_fire, HSTAR
      use OldTracer_mod, only: tr_RKD, tr_mm
      use OldTracer_mod, only: set_vol2mass, set_mass2vol
      use OldTracer_mod, only: initializeOldTracers
      use OldTracer_mod, only: oldAddTracer
      use OldTracer_mod, only: set_tr_mm, set_ntm_power
      use OldTracer_mod, only: set_t_qlimit
      use OldTracer_mod, only: set_needtrs
      use OldTracer_mod, only: set_trdecay
      use OldTracer_mod, only: set_HSTAR
      use OldTracer_mod, only: set_F0
      use OldTracer_mod, only: set_dodrydep
      use OldTracer_mod, only: dodrydep
      use OldTracer_mod, only: F0
      use OldTracer_mod, only: HSTAR
      use OldTracer_mod, only: set_trpdens
      use OldTracer_mod, only: set_trradius
      use OldTracer_mod, only: set_tr_wd_TYPE
      use OldTracer_mod, only: set_tr_RKD
      use OldTracer_mod, only: set_tr_DHD
      use OldTracer_mod, only: set_fq_aer
      use OldTracer_mod, only: set_rc_washt
      use OldTracer_mod, only: set_isDust
      use OldTracer_mod, only: set_tr_H2ObyCH4
      use OldTracer_mod, only: set_dowetdep
      use OldTracer_mod, only: set_trw0
      use OldTracer_mod, only: set_ntrocn
      use OldTracer_mod, only: set_conc_from_fw
      use OldTracer_mod, only: set_trglac
      use OldTracer_mod, only: set_ntisurfsrc
      use OldTracer_mod, only: set_trli0
      use OldTracer_mod, only: set_trsi0
      use OldTracer_mod, only: set_to_volume_MixRat, to_volume_MixRat
      use RunTimeControls_mod, only: tracers_special_shindell
      use RunTimeControls_mod, only: tracers_terp
      use RunTimeControls_mod, only: tracers_aerosols_soa
      use RunTimeControls_mod, only: shindell_strat_extra
      use RunTimeControls_mod, only: accmip_like_diags
      use RunTimeControls_mod, only: tracers_drydep
      use RunTimeControls_mod, only: tracers_tomas
      use RunTimeControls_mod, only: tracers_water
      use RunTimeControls_mod, only: tracers_special_o18
      use RunTimeControls_mod, only: tracers_gasexch_ocean_cfc
      use RunTimeControls_mod, only: tracers_gasexch_ocean_co2
      use RunTimeControls_mod, only: tracers_gasexch_land_co2
      use RunTimeControls_mod, only: tracers_special_lerner
      use RunTimeControls_mod, only: tracers_aerosols_koch
      use RunTimeControls_mod, only: sulf_only_aerosols
      use RunTimeControls_mod, only: tracers_aerosols_vbs
      use RunTimeControls_mod, only: tracers_aerosols_ocean
      use RunTimeControls_mod, only: tracers_nitrate
      use RunTimeControls_mod, only: cpp_tracers_dust => tracers_dust
      use RunTimeControls_mod, only: tracers_dust_silt4
      use RunTimeControls_mod, only: tracers_hetchem
      use RunTimeControls_mod, only: tracers_cosmo
      use RunTimeControls_mod, only: tracers_radon
      use RunTimeControls_mod, only: tracers_minerals
      use RunTimeControls_mod, only: tracers_quarzhem
      use RunTimeControls_mod, only: tracers_on
      use RunTimeControls_mod, only: tracers_ocean
      use RunTimeControls_mod, only: htap_like_diags
      use RunTimeControls_mod, only: tracers_air
      use RunTimeControls_mod, only: tracers_amp
      use TracerConstants_mod, only: H2O18
      use Tracer_mod, only: Tracer_type

      implicit none
      type (Tracer_type), pointer :: tracer
      external setDefaultSpec
      integer :: i, n

#ifdef TRACERS_SPECIAL_O18
      real*8 fracls
#endif

! call routine to read/set up regions and sectors for emissions:
      call setup_emis_sectors_regions()
      call initializeOldTracers(tracers, setDefaultSpec)

#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP) ||\
      (defined TRACERS_TOMAS)
!**** DMS, seasalt from offline fields
      call sync_param("OFFLINE_DMS_SS",OFFLINE_DMS_SS)
!**** seasalt from offline fields
      call sync_param("OFFLINE_SS",OFFLINE_SS)
#endif

! ***  BEGIN TRACER METADATA INITIALIZATION

#ifdef TRACERS_SPECIAL_Shindell
        if (tracers_special_shindell) then
          call SHINDELL_InitMetadata(tracer)
        end if
#endif

#ifdef TRACERS_SPECIAL_LERNER
        if (tracers_special_lerner) then
          call Lerner_InitMetadata(tracer, 1)
        end if
#endif

      if ((.not. tracers_amp) .and. tracers_water) then
        call  Water_setSpec('Water')
      end if
     
#ifdef TRACERS_SPECIAL_O18
        if (tracers_special_o18) then
          call H2O18_setSpec('H2O18')
          call HDO_setSpec('HDO')
        end if
#endif

      if (tracers_gasexch_ocean_cfc) then
        call  CFCn_setSpec('CFCn')
      end if

      if (tracers_gasexch_ocean_co2 .or. tracers_gasexch_land_co2) then
        call  CO2n_setSpec('CO2n')
      end if
      
#ifdef TRACERS_SPECIAL_LERNER
      if (tracers_special_lerner) then
        call Lerner_InitMetadata(tracer, 2)
        if (tracers_special_shindell) 
     &    call stop_model('contradictory tracer specs')
      end if
#endif

#ifdef TRACERS_AEROSOLS_KOCH
      if (tracers_aerosols_koch) then
        call KOCH_InitMetadata(tracer)
      end if
#endif

      if (tracers_aerosols_ocean) then
        call  OCocean_setSpec('OCocean') !Insoluble oceanic organic mass
      end if

      if (cpp_tracers_dust) then
        call  clay_setSpec('Clay')
        call  Silt1_setSpec('Silt1')
        call  Silt2_setSpec('Silt2')
        call  Silt3_setSpec('Silt3')
        if (tracers_dust_Silt4) call  Silt4_setSpec('Silt4')
      end if

#ifdef TRACERS_NITRATE
      if (tracers_nitrate) then
        call  NH3_setSpec('NH3')
        call  NH4_setSpec('NH4')
        call  NO3p_setSpec('NO3p')
      end if
#endif

      if (tracers_hetchem) then
        call  SO4_d1_setSpec('SO4_d1')
        call  SO4_d2_setSpec('SO4_d2')
        call  SO4_d3_setSpec('SO4_d3')
        if (tracers_nitrate) then
          call  N_d1_setSpec('N_d1')
          call  N_d2_setSpec('N_d2')
          call  N_d3_setSpec('N_d3')
        end if
      end if

      if (tracers_cosmo) then
#ifdef TRACERS_RADON
      if (tracers_radon) then
        call  Pb210_setSpec('Pb210')
      end if
#endif
        call  Be7_setSpec('Be7')
        call  Be10_setSpec('Be10')
#ifdef TRACERS_RADON
      if (tracers_special_lerner) then
        call  Rn222_setSpec('Rn222') ! duplicate with Lerner
      end if
#endif
        if (tracers_special_lerner) 
     &       call stop_model('contradictory tracer specs')
      end if

#ifdef TRACERS_MINERALS
       if (tracers_minerals) then
         call Minerals_InitMetadata(tracer)
       end if
#endif

#ifdef TRACERS_QUARZHEM
      if (tracers_quarzhem) then
        call Quarzhem_InitMetadata(tracer)
      end if
#endif

      if (tracers_air .or. htap_like_diags) then
        call  air_setSpec('air')
      end if

#ifdef TRACERS_AMP
      if (tracers_amp) then
        call AMP_InitMetadata(tracer)
      end if
#endif

#ifdef TRACERS_TOMAS
      if (tracers_tomas) then
        call TOMAS_InitMetadata(tracer)
      end if
#endif

! ***  END TRACER METADATA INITIALIZATION

      ! Generic tracer work
      ! All tracers must have been declared before reaching this point!!!
      ntm = getNumTracers(tracers)

      call set_param("NTM",NTM,'o')
      call set_param("TRNAME",trName(),ntm,'o')
      call printTracerNames(trName())

      ! Generic tracer work
      do i = 1, ntm
        if (tracers_water) then
!**** Tracers that are soluble or are scavenged or are water => wet dep
          if (tr_wd_type(i).eq.nWater.or.tr_wd_type(i) .EQ. nPART .or.
     &      tr_RKD(i).gt.0) then
            call set_dowetdep(i, .true.)
          end if
        end if
        if (tracers_drydep) then
!**** If tracers are particles or have non-zero HSTAR or F0 do dry dep:
!**** Any tracers that dry deposits needs the surface concentration:
          if(HSTAR(i).GT.0..OR.F0(i).GT.0..OR.tr_wd_type(i).eq.nPART) 
     &      then
            call set_dodrydep(i, .true.)
            call set_needtrs(i, .true.)
            if (tracers_water) then
              if (tr_wd_type(i).eq.nWATER) call stop_model
     &         ('A water tracer should not undergo dry deposition.',255)
            end if
          end if
        end if

        if (tracers_on) then
!**** Define the conversion from mass to volume units here
          call set_mass2vol(i, mair/tr_mm(i))
          call set_vol2mass(i, tr_mm(i)/mair)
          call set_to_conc(i, 0)
          if (tracers_special_shindell) then
!**** Aerosol tracer output should be mass mixing ratio
            select case (tr_wd_TYPE(i))
            case (nGAS)
              call set_to_volume_MixRat(i, 1) !gas output to volume mixing ratio
            case (nPART, nWATER)
              call set_to_volume_MixRat(i, 0) ! aerosol/water output to mass mixing ratio
            case default
              call set_to_volume_MixRat(i, 0) !default output to mass mixing ratio
            end select
          end if
          if (tracers_gasexch_ocean_co2 .or. tracers_gasexch_land_co2) 
     &      then
            call set_to_volume_MixRat(i, 1) !gas output to volume mixing ratio
          end if
        end if
        
      end do

      contains

       subroutine printTracerNames(tracerNames)
       use domain_decomp_atm, only: am_i_root

       character(len=*) :: tracerNames(:)
       integer :: i
      
       if (am_i_root()) then
          do i = 1, size(tracerNames)
             write(6,*) 'TRACER',i,trim(tracerNames(i))
          end do
       end if
      
       end subroutine printTracerNames

      subroutine Water_setSpec(name)
      character(len=*), intent(in) :: name

      n = oldAddTracer(name)
      n_Water = n
      if (tracers_water .or. tracers_ocean) then
        call set_trw0(n, 1.d+0)
        call set_ntrocn(n, 0)
      end if
      if (tracers_water) then
        call set_ntm_power(n, -4)
        call set_tr_mm(n, mwat)
        call set_needtrs(n,  .true.)
        call set_tr_wd_type(n, nWater)
        call set_trli0(n, 1.d+0)
        call set_trsi0(n, 1.d+0)
        call set_tr_H2ObyCH4(n, 1.d+0)
      end if
      if (tracers_ocean) call set_trglac(n, 1.d+0)
#ifdef TRACERS_SPECIAL_O18
      if (tracers_special_o18) iso_index(n) = 1 ! indexing for isotopic fractionation calcs
#endif
      end subroutine Water_setSpec

#ifdef TRACERS_SPECIAL_O18
      subroutine H2O18_setSpec(name)
      character(len=*), intent(in) :: name

      n = oldAddTracer(name)
      n_H2O18 = n
      call set_ntm_power(n, -7)
      call set_tr_mm(n, H2O18%molMass)
      call set_needtrs(n,  .true.)
      call set_tr_wd_type(n, nwater)
      iso_index(n) = 2          ! indexing for isotopic fractionation calcs
      call set_trw0(n, 2.228d-3   ) ! SMOW mass ratio of water molecules
      call set_trli0(n, 0.980d0*trw0(n)  ) ! d=-20
      call set_trsi0(n, fracls(n)*trw0(n))
      call set_tr_H2ObyCH4(n, trw0(n)*1.023d0 ) ! d=+23 (ie. no frac from O2)
      call set_ntrocn(n, -3)
      if (tracers_ocean)  call set_trglac(n, trw0(n)*0.98d0 ) ! d=-20
      end if

      end subroutine H2O18_setSpec

      subroutine HDO_setSpec(name)
      character(len=*), intent(in) :: name

      n = oldAddTracer(name)
      n_HDO = n
      call set_ntm_power(n, -8)
      call set_tr_mm(n, 19d0)
      call set_needtrs(n,  .true.)
      call set_tr_wd_type(n, nwater)
      iso_index(n) = 3          ! indexing for isotopic fractionation calcs
      call set_trw0(n, 3.29d-4    ) ! SMOW mass ratio of water molecules
      call set_trli0(n, 0.830d0*trw0(n)  ) ! d=-170
      call set_trsi0(n, fracls(n)*trw0(n))
      call set_tr_H2ObyCH4(n, trw0(n)*0.93d0  ) ! d=-70
      call set_ntrocn(n, -4)
      if (tracers_ocean) call set_trglac(n, trw0(n)*0.84d0   ) ! d=-160
      end if

      end subroutine HDO_setSpec

      subroutine HTO_setSpec(name)
      character(len=*), intent(in) :: name

      n = oldAddTracer(name)
      n_HTO = n
      call set_ntm_power(n, -18)
      call set_tr_mm(n, 20d0)
      call set_needtrs(n,  .true.)
      call set_tr_wd_type(n, nwater)
      iso_index(n) = 4          ! indexing for isotopic fractionation calcs
      call set_trw0(n, 0d0)     !2.22d-18   ) ! SMOW mass ratio of water molecules
      call set_trli0(n, 0d0)
      call set_trsi0(n, 0d0)
      call set_tr_H2ObyCH4(n, 0d0)
      call set_trdecay(n,  1.77d-9) ! =5.59d-2 /yr
      call set_ntrocn(n, -18)
      if (tracers_ocean) call set_trglac(n, 0d0)

      end subroutine HTO_setSpec

      subroutine H2O17_setSpec(name)
      character(len=*), intent(in) :: name

      n = oldAddTracer(name)
      n_H2O17 = n
      call set_ntm_power(n, -7)
      call set_tr_mm(n, 19d0)
      call set_needtrs(n,  .true.)
      call set_tr_wd_type(n, nwater)
      iso_index(n) = 5          ! indexing for isotopic fractionation calcs
      call set_trw0(n, 4.020d-5   ) ! SMOW mass ratio of water molecules
      call set_trli0(n, 0.98937d0*trw0(n)  ) ! d=-10.63 D17O=0
      call set_trsi0(n, fracls(n)*trw0(n))
      call set_tr_H2ObyCH4(n, trw0(n)*1.011596d0 ) ! d=+11.596 (some frac from O2)
      call set_ntrocn(n, -3)
      if (tracers_ocean) call set_trglac(n, trw0(n)*0.98937d0   ) ! d=-10.63 D17O=

      end subroutine H2O17_setSpec

#endif  /* TRACERS_SPECIAL_O18 */

      subroutine CFCn_setSpec(name)
      character(len=*), intent(in) :: name

      n = oldAddTracer(name)
      n_CFCn = n
      call set_ntm_power(n, -12)
!     !call set_tr_mm(n, 136.0d0    !NCAR value
      call set_tr_mm(n, 137.37d0) !note units are in gr
      call set_ntsurfsrc(n,  1)
      call set_needtrs(n, .true.)

      end subroutine CFCn_setSpec

      subroutine CO2n_setSpec(name)
      character(len=*), intent(in) :: name

      n = oldAddTracer(name)
      n_CO2n = n
      call set_ntm_power(n, -6)
      call set_tr_mm(n, 44.d0)  !grams
      call set_t_qlimit(n,  .false.)
      call set_ntsurfsrc(n,  1)
      call set_needtrs(n, .true.)

      end subroutine CO2n_setSpec

      subroutine OCocean_setSpec(name)
      character(len=*), intent(in) :: name

      n = oldAddTracer(name)
      n_OCocean = n
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 15.6d0)
      call set_trpdens(n, 1.5d3) !kg/m3
      call set_trradius(n, 4.4d-7) !m, same as seasalt1
      call set_fq_aer(n, 1.0d0) ! same as seasalt
      call set_tr_wd_type(n, nPART)

      end subroutine OCocean_setSpec

      subroutine Clay_setSpec(name)
      character(len=*), intent(in) :: name

      n = oldAddTracer(name)
      n_clay=n
      call set_ntm_power(n, -9)
      call set_trpdens(n, 2.5d3)
      if (tracers_drydep) call set_trradius(n, 0.46D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)

      end subroutine Clay_setSpec

      subroutine Silt1_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_Silt1=n
      call set_ntm_power(n, -9)
      call set_trpdens(n, 2.65d3)
      if (tracers_drydep) call set_trradius(n, 1.47D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
      end subroutine Silt1_setSpec

      subroutine Silt2_setSpec(name)
      character(len=*), intent(in) :: name

      n = oldAddTracer(name)
      n_Silt2=n
      call set_ntm_power(n, -9)
      call set_trpdens(n, 2.65d3)
      if (tracers_drydep) call set_trradius(n, 2.94D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)

      end subroutine Silt2_setSpec

      subroutine Silt3_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_Silt3 = n
      call set_ntm_power(n, -9)
      call set_trpdens(n, 2.65d3)
      if (tracers_drydep) call set_trradius(n, 5.88D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
      end subroutine Silt3_setSpec

      subroutine Silt4_setSpec(name)
      character(len=*), intent(in) :: name

      n = oldAddTracer(name)
      n_Silt4=n
      call set_ntm_power(n, -9)
      call set_trpdens(n, 2.65d3)
      if (tracers_drydep) call set_trradius(n, 11.77D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)

      end subroutine Silt4_setSpec

      subroutine NO3p_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_NO3p = n
      call set_ntsurfsrc(n,  0)
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 62.d0)
      call set_trpdens(n, 1.7d3)
      call set_trradius(n, 3.d-7)
      call set_fq_aer(n, 1.0d0   ) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, npart)
      end subroutine NO3p_setSpec

      subroutine SO4_d1_setSpec(name)
      character(len=*), intent(in) :: name

      n = oldAddTracer(name)
      n_SO4_d1 = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 96.d0)  !!!! Sulfat
      call set_trpdens(n, 2.5d3) !kg/m3 this is clay density
      call set_trradius(n, 0.75D-06 ) !m
      call set_fq_aer(n, 1.d0   ) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, npart)

      end subroutine SO4_d1_setSpec

      subroutine SO4_d2_setSpec(name)
      character(len=*), intent(in) :: name

      n = oldAddTracer(name)
      n_SO4_d2 = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 96.d0)
      call set_trpdens(n, 2.65d3) !kg/m3 this is Silt1 value
      call set_trradius(n, 2.2D-06 ) !m
      call set_fq_aer(n, 1.d0   ) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, npart)

      end subroutine SO4_d2_setSpec

      subroutine SO4_d3_setSpec(name)
      character(len=*), intent(in) :: name

      n = oldAddTracer(name)
      n_SO4_d3 = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 96.d0)
      call set_trpdens(n, 2.65d3) !this is Silt2 value
      call set_trradius(n, 4.4D-06 ) !m this is Silt2 value
      call set_fq_aer(n, 1.d0   ) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, npart)

      end subroutine SO4_d3_setSpec

      subroutine N_d1_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_N_d1 = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 62.d+0) ! NO3
      call set_trpdens(n, 2.5d3) !kg/m3 this is clay density
      call set_trradius(n, 0.75D-06 ) !m
      call set_fq_aer(n, 1.d0  ) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, npart)
      end subroutine N_d1_setSpec

      subroutine N_d2_setSpec(name)
      character(len=*), intent(in) :: name

      n = oldAddTracer(name)
      n_N_d2 = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 62.d+0)
      call set_trpdens(n, 2.65d3) !kg/m3 this is Silt1 value
      call set_trradius(n, 2.2D-06 ) !m
      call set_fq_aer(n, 1.d0  ) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, npart)

      end subroutine N_d2_setSpec

      subroutine N_d3_setSpec(name)
      character(len=*), intent(in) :: name

      n = oldAddTracer(name)
      n_N_d3 = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 62.d0)
      call set_trpdens(n, 2.65d3) !this is Silt2 value
      call set_trradius(n, 4.4D-06 ) !m this is Silt2 value
      call set_fq_aer(n, 1.d0  ) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, npart)

      end subroutine N_d3_setSpec

      subroutine Pb210_setSpec(name)
      character(len=*), intent(in) :: name

      n = oldAddTracer(name)
      n_Pb210 = n
      call set_ntm_power(n, -23)
      call set_tr_mm(n, 210.d0)
      call set_trdecay(n,  9.85d-10)
      call set_trpdens(n, 1.7d3) !kg/m3 this is SO4 value
      call set_trradius(n, 3.d-7  ) !again S04 value
      call set_fq_aer(n, 1.d0   ) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, npart) ! same as SO4

      end subroutine Pb210_setSpec

      subroutine Be7_setSpec(name)
      character(len=*), intent(in) :: name

      n = oldAddTracer(name)
      n_Be7 = n
      call set_ntm_power(n, -23) ! power of ten for tracer
      call set_tr_mm(n, 7.d0)
      call set_trdecay(n,  1.51d-7)
      call set_trpdens(n, 1.7d3) !kg/m3 this is SO4 value
      call set_trradius(n, 1.d-7  ) !appropriate for stratosphere
      call set_fq_aer(n, 1.d0   ) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, npart) ! same as SO4

      end subroutine Be7_setSpec

      subroutine Be10_setSpec(name)
      character(len=*), intent(in) :: name

      n = oldAddTracer(name)
      n_Be10 = n
      call set_ntm_power(n, -23)
      call set_tr_mm(n, 10.d0)
      call set_trpdens(n, 1.7d3) !kg/m3 this is SO4 value
      call set_trradius(n, 1.d-7  ) !appropriate for stratosphere
      call set_fq_aer(n, 1.d0   ) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, npart) ! same as SO4

      end subroutine Be10_setSpec


      subroutine Air_setSpec(name)
      character(len=*), intent(in) :: name

      n = oldAddTracer(name)
      n_Air = n
      call set_ntm_power(n, -2)
      call set_tr_mm(n, mair)

      end subroutine Air_setSpec

      end subroutine initTracerMetadata


!------------------------------------------------------------------------------
      subroutine laterInitTracerMetadata()
!------------------------------------------------------------------------------
      USE MODEL_COM, only: itime,master_yr
      use OldTracer_mod, only: itime_tr0
      use OldTracer_mod, only: set_itime_tr0
      USE TRACER_COM, only: NTM, tracers, syncProperty
      use TRACER_COM, only: coupled_chem
      use Dictionary_mod, only: sync_param,is_set_param,get_param
#ifdef TRACERS_WATER
      use TRDIAG_com, only: to_per_mil
#endif
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP) ||\
      (defined TRACERS_TOMAS)
      USE AEROSOL_SOURCES, only: tune_ss1, tune_ss2, BBinc
#endif
      use TRDIAG_COM, only: diag_rad
      use TRACER_COM, only: ntm ! should be available by this procedure call
#ifdef TRACERS_WATER
#ifdef TRDIAG_WETDEPO
      USE CLOUDS, ONLY : diag_wetdep
#endif
#endif /* TRACERS_WATER */
#ifdef TRACERS_SPECIAL_Shindell
      use tracer_sources, only: aircraft_Tyr1,aircraft_Tyr2
      USE TRCHEM_Shindell_COM,only:LCOalt,PCOalt,
     &     CH4altINT,CH4altINX,LCH4alt,PCH4alt,
     &     CH4altX,CH4altT,ch4_init_sh,ch4_init_nh,scale_ch4_IC_file,
     &     OxICIN,OxIC,OxICINL,OxICL,
     &     fix_CH4_chemistry,which_trop,PI_run,PIratio_N,PIratio_CO_T,
     &     PIratio_CO_S,PIratio_other,allowSomeChemReinit,
     &     CH4ICIN,CH4ICX,CH4ICINL,CH4ICL,rad_FL,use_rad_ch4,
     &     COICIN,COIC,COICINL,COICL,Lmax_rad_O3,Lmax_rad_CH4
     &     ,BrOxaltIN,ClOxaltIN,ClONO2altIN,HClaltIN,BrOxalt,
     &     ClOxalt,ClONO2alt,HClalt,N2OICIN,N2OICX,N2OICINL,N2OICL,
     &     CFCICIN,CFCIC,CFCICINL,CFCICL,PIratio_N2O,PIratio_CFC,
     &     use_rad_n2o,use_rad_cfc,cfc_rad95,PltOx,Tpsc_offset_N,
     &     Tpsc_offset_S
#ifdef INTERACTIVE_WETLANDS_CH4
      USE TRACER_SOURCES, only:int_wet_dist,topo_lim,sat_lim,gw_ulim,
     & gw_llim,sw_lim,exclude_us_eu,nn_or_zon,ice_age,nday_ch4,max_days,
     & ns_wet,nra_ch4
#endif
#ifdef BIOGENIC_EMISSIONS
      use biogenic_emis, only: base_isopreneX
#endif
#endif /* TRACERS_SPECIAL_Shindell */
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP) ||\
    (defined TRACERS_TOMAS) 
      use TRACER_COM, only: aer_int_yr
      USE TRACER_COM, only: offline_dms_ss, offline_ss
#endif
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
      (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP)  ||\
      (defined TRACERS_TOMAS)
      use tracers_dust,only : imDust,prefDustSources,fracClayPDFscheme
     &     ,fracSiltPDFscheme
#endif
#ifdef TRACERS_AMP
      USE AMP_AEROSOL, only: AMP_DIAG_FC, AMP_RAD_KEY
#endif

      USE TRACER_COM, only: ef_fact3d, no_emis_over_ice
      use Model_com, only: itime
      implicit none
      integer :: n

C**** 
C**** Set some documentary parameters in the database
C**** 
      do n = 1, NTM
        call set_itime_tr0(n, itime)
      end do
      call syncProperty(tracers, "itime_tr0", set_itime_tr0,itime_tr0())

      call sync_param( "COUPLED_CHEM", COUPLED_CHEM )

#ifdef TRACERS_ON
#ifdef TRACERS_SPECIAL_Lerner
!TLC - LERNER tracer will need some changes
      LERNER tracers not supported with current changes
! Lerner defaults
      n_MPtable = 0
      tcscale = 0.
#endif
#endif /* TRACERS_ON */

C**** Synchronise tracer related parameters from rundeck

 
#ifdef TRACERS_WATER
C**** Decide on water tracer conc. units from rundeck if it exists
      call sync_param("to_per_mil",to_per_mil,ntm)
#endif
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP) ||\
      (defined TRACERS_TOMAS)
      call sync_param("tune_ss1",tune_ss1)
      call sync_param("tune_ss2",tune_ss2)
      call sync_param("BBinc",BBinc)
C**** determine year of emissions
      if (is_set_param("aer_int_yr")) then
        call get_param("aer_int_yr",aer_int_yr)
      else
        aer_int_yr=master_yr
      endif
#endif
#ifdef TRACERS_AEROSOLS_VBS
      call sync_param("VBSemifact",VBSemifact,vbs_tr%nbins)
#endif
#ifdef TRACERS_SPECIAL_O18
C**** set super saturation parameter for isotopes if needed
      call sync_param("supsatfac",supsatfac)
#endif
#ifdef TRACERS_ON
      CALL sync_param("diag_rad",diag_rad)
#if (defined TRACERS_WATER) && (defined TRDIAG_WETDEPO)
      CALL sync_param("diag_wetdep",diag_wetdep)
#endif
!     not params call sync_param("trans_emis_overr_day",trans_emis_overr_day)
!     not params call sync_param("trans_emis_overr_yr", trans_emis_overr_yr )
#endif /* TRACERS_ON */
#ifdef TRACERS_SPECIAL_Shindell
      call sync_param("allowSomeChemReinit",allowSomeChemReinit)
      call sync_param("which_trop",which_trop)
      if (is_set_param("PI_run")) then
        call get_param("PI_run",PI_run)
      else
        if (master_yr == 1850) then
          PI_run=1
        else
          PI_run=0
        endif
      endif
      call sync_param("PIratio_N",PIratio_N)
      call sync_param("PIratio_CO_T",PIratio_CO_T)
      call sync_param("PIratio_CO_S",PIratio_CO_S)
      call sync_param("PIratio_other",PIratio_other)
      call sync_param("rad_FL",rad_fl)
      call sync_param("use_rad_ch4",use_rad_ch4)
      call sync_param("Lmax_rad_O3",Lmax_rad_O3)
      call sync_param("Lmax_rad_CH4",Lmax_rad_CH4)
      if (is_set_param("aircraft_Tyr1")) then
        call get_param("aircraft_Tyr1",aircraft_Tyr1)
      else
        if (master_yr == 0) then
          call stop_model("Please provide aircraft_Tyr1 via the "//
     .                    "rundeck", 255)
        else
          aircraft_Tyr1=master_yr
        endif
      endif
      if (is_set_param("aircraft_Tyr2")) then
        call get_param("aircraft_Tyr2",aircraft_Tyr2)
      else
        if (master_yr == 0) then
          call stop_model("Please provide aircraft_Tyr2 via the "//
     .                    "rundeck", 255)
        else
          aircraft_Tyr2=master_yr
        endif
      endif
      call sync_param("use_rad_n2o",use_rad_n2o)
      call sync_param("use_rad_cfc",use_rad_cfc)
      call sync_param("PIratio_N2O",PIratio_N2O)
      call sync_param("PIratio_CFC",PIratio_CFC)
      call sync_param("PltOx",PltOx)
      call sync_param("Tpsc_offset_N",Tpsc_offset_N)
      call sync_param("Tpsc_offset_S",Tpsc_offset_S)
#ifdef BIOGENIC_EMISSIONS
      call sync_param("base_isopreneX",base_isopreneX)
#endif
#ifdef INTERACTIVE_WETLANDS_CH4
      call sync_param("ice_age",ice_age)
      call sync_param("ns_wet",ns_wet)
      call sync_param("int_wet_dist",int_wet_dist)
      call sync_param("topo_lim",topo_lim)
      call sync_param("sat_lim",sat_lim)
      call sync_param("gw_ulim",gw_ulim)
      call sync_param("gw_llim",gw_llim)
      call sync_param("sw_lim",sw_lim)
      call sync_param("exclude_us_eu",exclude_us_eu)
      call sync_param("nn_or_zon",nn_or_zon)
      do n=1,nra_ch4
        if(nday_ch4(n) > max_days .or. nday_ch4(n) < 1)
     &       call stop_model('nday_ch4 out of range',255)
      end do
#endif

#endif /* TRACERS_SPECIAL_Shindell */

#if (defined TRACERS_AMP)
C**** Decide on how many times Radiation is called for aerosols once or nmode, default one call
      call sync_param("AMP_DIAG_FC",AMP_DIAG_FC)
C**** Decide Radiative Mixing Rules - Volume - Core Shell - Maxwell Garnett, default Volume
      call sync_param("AMP_RAD_KEY",AMP_RAD_KEY)
#endif
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
      (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP)  ||\
      (defined TRACERS_TOMAS)
C**** decide on AEROCOM or interactive emissions
      CALL sync_param('imDUST',imDUST)
      call sync_param('prefDustSources', prefDustSources)
      call sync_param('fracClayPDFscheme', fracClayPDFscheme)
      call sync_param('fracSiltPDFscheme', fracSiltPDFscheme)
#endif
#ifdef TRACERS_QUARZHEM
      call sync_param( 'frHemaInQuarAggr', frHemaInQuarAggr )
      call sync_param( 'pureByTotalHematite', pureByTotalHematite )
#endif

#if (defined TRACERS_COSMO)
C**** get rundeck parameter for cosmogenic source factor
      call sync_param("be7_src_param", be7_src_param)
#endif
      call sync_param("no_emis_over_ice",no_emis_over_ice)

!     initialize 3D source factors:
      ef_fact3d(:,:)=1.d0

      end subroutine laterInitTracerMetadata

