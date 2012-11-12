#include "rundeck_opts.h"
      subroutine setDefaultSpec(n, tracer)
      use Dictionary_mod, only: sync_param
      use RunTimeControls_mod, only: tracers_amp
      use RunTimeControls_mod, only: tracers_tomas
      use OldTracer_mod, only: trName, do_fire
      use OldTracer_mod, only: nBBsources, set_nBBsources
      use DOMAIN_DECOMP_ATM, only: am_i_root
      use TRACER_COM, only: tracers
      use TRACER_COM, only: set_ntsurfsrc
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

      subroutine initTracerMetadata()
      USE DOMAIN_DECOMP_ATM, only: AM_I_ROOT
      USE CONSTANT, only: mair,mwat,pi
#ifdef TRACERS_AEROSOLS_SOA
     &     ,gasc
#endif  /* TRACERS_AEROSOLS_SOA */
      USE RESOLUTION, only : jm,lm
      USE MODEL_COM, only: dtsrc, master_yr
      USE ATM_COM, only: pmidl00
      USE GEOM, only: axyp,byaxyp
      USE ATM_COM, only: am     ! Air mass of each box (kg/m^2)
      use OldTracer_mod, only: trName, nGAS, nWater, nPart
      use OldTracer_mod, only: tr_wd_type, do_fire, HSTAR
      use OldTracer_mod, only: tr_RKD, nBBsources, tr_mm
      use OldTracer_mod, only: set_vol2mass, set_mass2vol
      use TracerBundle_mod, only: getNumTracers
      USE TRACER_COM, only: ntm, ef_fact3d, no_emis_over_ice
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP) ||\
    (defined TRACERS_TOMAS) 
      use TRACER_COM, only: aer_int_yr
      USE TRACER_COM, only: offline_dms_ss, offline_ss
#endif
      USE TRACER_COM, only: tracers
      USE TRACER_COM, only: num_sectors
      use TRACER_COM, only: set_ntsurfsrc, ntsurfsrc
      use TRACER_COM, only:
#ifdef TRACERS_TOMAS
      use TRACER_COM, only: n_ASO4, n_ANACL, n_AECOB, n_AECIL
      use TRACER_COM, only: n_AOCOB, n_AOCIL, n_ADUST, n_ANUM
      use TRACER_COM, only: n_SOAgas, xk, n_AH2O, n_ASO4, nbins
#endif
      use TRACER_COM, only: 
     *                 n_SF6_c,                                        
     *     n_Air,    n_SF6,   n_Rn222, n_CO2,      n_N2O,
     *     n_CFC11,  n_14CO2, n_CH4,   n_O3,       n_water,
     *     n_H2O18,  n_HDO,   n_HTO,   n_Ox,       n_NOx,
     *     n_N2O5,   n_HNO3,  n_H2O2,  n_CH3OOH,   n_HCHO,
     *     n_HO2NO2, n_CO,    n_PAN,   n_H2O17,
     *     n_Isoprene, n_AlkylNit, n_Alkenes, n_Paraffin,
     *     n_stratOx, n_Terpenes,n_codirect,
     *     n_isopp1g,n_isopp1a,n_isopp2g,n_isopp2a,
     *     n_apinp1g,n_apinp1a,n_apinp2g,n_apinp2a,
     *     n_DMS,    n_MSA,   n_SO2,   n_SO4,    n_H2O2_s,
     *     n_ClOx,   n_BrOx,  n_HCl,   n_HOCl,   n_ClONO2,
     *     n_HBr,    n_HOBr,  n_BrONO2,n_CFC,    n_GLT,
     *     n_Pb210,n_Be7,   n_Be10,
     .     n_CFCn,   n_CO2n,  n_Age,
     *     n_seasalt1,  n_seasalt2, n_SO4_d1,  n_SO4_d2,
     *     n_SO4_d3,n_N_d1,  n_N_d2,  n_N_d3,
     &     n_NH3,   n_NH4,   n_NO3p,
     *     n_BCII,  n_BCIA,  n_BCB,
     *     n_OCII,  n_OCIA,  n_OCB,
     *     n_vbsGm2, n_vbsGm1, n_vbsGz,  n_vbsGp1, n_vbsGp2,
     *     n_vbsGp3, n_vbsGp4, n_vbsGp5, n_vbsGp6,
     *     n_vbsAm2, n_vbsAm1, n_vbsAz,  n_vbsAp1, n_vbsAp2,
     *     n_vbsAp3, n_vbsAp4, n_vbsAp5, n_vbsAp6,
     *     n_OCocean,
     &     n_clay,   n_silt1, n_silt2, n_silt3, n_silt4,
     &     n_clayilli,n_claykaol,n_claysmec,n_claycalc,
     &     n_clayquar,n_sil1quar,n_sil1feld,n_sil1calc,
     &     n_sil1hema,n_sil1gyps,n_sil2quar,n_sil2feld,
     &     n_sil2calc,n_sil2hema,n_sil2gyps,n_sil3quar,
     &     n_sil3feld,n_sil3calc,n_sil3hema,n_sil3gyps,
     &     n_sil1quhe,n_sil2quhe,n_sil3quhe,
     *     n_M_NO3,   n_M_NH4,   n_M_H2O,   n_M_AKK_SU,
     *     n_N_AKK_1, n_M_ACC_SU,n_N_ACC_1, n_M_DD1_SU,
     *     n_M_DD1_DU,n_N_DD1_1, n_M_DS1_SU,n_M_DS1_DU,
     *     n_N_DS1_1 ,n_M_DD2_SU,n_M_DD2_DU,n_N_DD2_1 ,
     *     n_M_DS2_SU,n_M_DS2_DU,n_N_DS2_1 ,n_M_SSA_SU,
     *     n_M_SSA_SS,n_M_SSC_SS,
     *     n_M_OCC_SU,n_M_OCC_OC,n_N_OCC_1 ,
     *     n_M_BC1_SU,n_M_BC1_BC,n_N_BC1_1 ,n_M_BC2_SU,
     *     n_M_BC2_BC,n_N_BC2_1 ,n_M_BC3_SU,n_M_BC3_BC,
     *     n_N_BC3_1 ,n_M_DBC_SU,n_M_DBC_BC,n_M_DBC_DU,
     *     n_N_DBC_1 ,n_M_BOC_SU,n_M_BOC_BC,n_M_BOC_OC,
     *     n_N_BOC_1, n_M_BCS_SU,n_M_BCS_BC,n_N_BCS_1 ,
     *     n_M_MXX_SU,n_M_MXX_BC,n_M_MXX_OC,n_M_MXX_DU,
     *     n_M_MXX_SS,n_N_MXX_1 ,n_M_OCS_SU,n_M_OCS_OC,
     *     n_N_OCS_1,n_M_SSS_SS,n_M_SSS_SU,
     *     n_H2SO4, n_N_SSA_1, n_N_SSC_1
#ifdef TRACERS_ON
      USE TRDIAG_COM
#endif
      USE Dictionary_mod
#ifdef TRACERS_SPECIAL_Lerner
      USE TRACERS_MPchem_COM, only: n_MPtable,tcscale
!     @dbparam dsol describes portion of solar cycle being modeled for linoz
!     @+      +1.0 = solar max, 0.0 = neutral, -1.0 = solar min
      USE LINOZ_CHEM_COM, only: dsol
#endif
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
      (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP)  ||\
      (defined TRACERS_TOMAS)
      use tracers_dust,only : imDust,prefDustSources,fracClayPDFscheme
     &     ,fracSiltPDFscheme
#endif
#ifdef TRACERS_QUARZHEM
     &     ,DensityHematite, DensityQuartz, FreeFe, frHemaInQuarAggr,
     &     pureByTotalHematite
#endif
#ifdef TRACERS_SPECIAL_Shindell
      USE TRCHEM_Shindell_COM,only:LCOalt,PCOalt,
     &     CH4altINT,CH4altINX,LCH4alt,PCH4alt,
     &     CH4altX,CH4altT,ch4_init_sh,ch4_init_nh,scale_ch4_IC_file,
     &     OxICIN,OxIC,OxICINL,OxICL,
     &     fix_CH4_chemistry,PI_run,
     &     allowSomeChemReinit,
     &     CH4ICIN,CH4ICX,CH4ICINL,CH4ICL,rad_FL,use_rad_ch4,
     &     COICIN,COIC,COICINL,COICL,Lmax_rad_O3
     &     ,BrOxaltIN,ClOxaltIN,ClONO2altIN,HClaltIN,BrOxalt,
     &     ClOxalt,ClONO2alt,HClalt,N2OICIN,N2OICX,N2OICINL,N2OICL,
     &     CFCICIN,CFCIC,CFCICINL,CFCICL,
     &     use_rad_cfc,cfc_rad95,PltOx,Tpsc_offset_N,
     &     Tpsc_offset_S
#ifdef INTERACTIVE_WETLANDS_CH4
      USE TRACER_SOURCES, only:int_wet_dist,topo_lim,sat_lim,gw_ulim,
     &     gw_llim,sw_lim,exclude_us_eu,nn_or_zon,ice_age,nday_ch4,max_days,
     &     ns_wet,nra_ch4
#endif
#endif /* TRACERS_SPECIAL_Shindell */
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
#if (defined TRACERS_AMP)
      USE AERO_PARAM, only: DG_DD1, DG_DD2, DG_AKK,
     &     DG_DS1, DG_DS2, DG_SSA, DG_SSC, DG_ACC,
     &     DG_SSS, DG_OCC, DG_BC1, DG_BC2, DG_BC3,
     &     DG_DBC, DG_BOC, DG_BCS, DG_OCS, DG_MXX,

     &     SOLU_DD1, SOLU_DD2, SOLU_AKK, SOLU_ACC,
     &     SOLU_DS1, SOLU_DS2, SOLU_SSA, SOLU_SSC,
     &     SOLU_SSS, SOLU_OCC, SOLU_BC1, SOLU_BC2, SOLU_BC3,
     &     SOLU_DBC, SOLU_BOC, SOLU_BCS, SOLU_OCS, SOLU_MXX

      USE AERO_ACTV, only: DENS_SULF, DENS_DUST,
     &     DENS_SEAS, DENS_BCAR, DENS_OCAR
      USE AERO_CONFIG, only: nbins
      USE AMP_AEROSOL, only: AMP_DIAG_FC, AMP_RAD_KEY
      USE AERO_COAG, only : SETUP_KIJ
      USE AERO_SETUP
      USE AERO_NPF, only: SETUP_NPFMASS
      USE AERO_DIAM, only: SETUP_DIAM
#endif
#ifdef TRACERS_GASEXCH_ocean_CO2
      USE obio_forc, only : atmCO2
#endif
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP) ||\
      (defined TRACERS_TOMAS)
      USE AEROSOL_SOURCES, only: tune_ss1, tune_ss2, om2oc, BBinc
#endif
#ifdef TRACERS_TOMAS
      use TOMAS_AEROSOL, only : binact10,binact02,
     &     fraction10,fraction02
#endif
      USE FILEMANAGER, only: openunit,closeunit,nameunit
      use OldTracer_mod, only: initializeOldTracers
      use OldTracer_mod, only: oldAddTracer
      use OldTracer_mod, only: set_tr_mm, set_ntm_power
      use OldTracer_mod, only: set_t_qlimit
      use OldTracer_mod, only: set_needtrs
      use OldTracer_mod, only: set_trdecay
      use OldTracer_mod, only: set_mass2vol
      use OldTracer_mod, only: set_vol2mass
      use OldTracer_mod, only: set_HSTAR
      use OldTracer_mod, only: set_F0
      use OldTracer_mod, only: set_dodrydep

      use OldTracer_mod, only: dodrydep
      use OldTracer_mod, only: F0
      use OldTracer_mod, only: HSTAR

      use OldTracer_mod, only: set_do_fire
      use OldTracer_mod, only: set_nBBsources
      use OldTracer_mod, only: set_emisPerFireByVegType
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

#if (defined TRACERS_OCEAN) &&  !defined(TRACERS_OCEAN_INDEP)
!     atmosphere copies atmosphere-declared tracer info to ocean
!     so that the ocean can "inherit" it without referencing atm. code
      use ocn_tracer_com, only : 
     &     n_Water_ocn      => n_Water,
     &     ntrocn_ocn       => ntrocn,
     &     to_per_mil_ocn   => to_per_mil,
     &     t_qlimit_ocn     => t_qlimit,
     &     conc_from_fw_ocn => conc_from_fw,
     &     trdecay_ocn      => trdecay,
     &     trw0_ocn         => trw0
#endif
      USE FLUXES, only : atmocn
      use RunTimeControls_mod, only: tracers_special_shindell
      use RunTimeControls_mod, only: tracers_terp
      use RunTimeControls_mod, only: tracers_aerosols_soa
      use RunTimeControls_mod, only: shindell_strat_extra
      use RunTimeControls_mod, only: accmip_like_diags
      use RunTimeControls_mod, only: tracers_amp
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
      use RunTimeControls_mod, only: htap_like_diags
      use RunTimeControls_mod, only: tracers_air
      use RunTimeControls_mod, only: tracers_amp
      use RunTimeControls_mod, only: tracers_amp_m1
      use RunTimeControls_mod, only: tracers_amp_m2
      use RunTimeControls_mod, only: tracers_amp_m3
      use RunTimeControls_mod, only: tracers_amp_m4
      use RunTimeControls_mod, only: tracers_amp_m5
      use RunTimeControls_mod, only: tracers_amp_m6
      use RunTimeControls_mod, only: tracers_amp_m7
      use RunTimeControls_mod, only: tracers_amp_m8
      use TracerConstants_mod, only: H2O18
      implicit none
      external setDefaultSpec
      integer :: l,k,n,kr,m,ns
#ifdef TRACERS_SPECIAL_O18
      real*8 fracls
#endif
#ifdef TRACERS_TOMAS
      real*8 :: TOMAS_dens,TOMAS_radius
#endif
#if (defined TRACERS_WATER) || (defined TRACERS_DRYDEP)
!@param convert_HSTAR converts from mole/Joule to mole/(L*atm)
      real*8, parameter :: convert_HSTAR = 1.01325d2
#endif
#ifdef TRACERS_SPECIAL_Shindell
!     @var iu_data unit number
!     @var title header read in from file
      integer iu_data,i,j,nq
      character*80 title
      character(len=300) :: out_line
      real*8, dimension(6) :: temp_ghg
      integer :: temp_year
#endif /* TRACERS_SPECIAL_Shindell */

#if defined(TRACERS_GASEXCH_ocean) && defined(TRACERS_GASEXCH_ocean_CFC)
      integer i, iu_data
#endif
#if (                           !defined(TRACERS_GASEXCH_ocean_CO2)) && defined(TRACERS_GASEXCH_land_CO2)
      real*8 :: atmCO2 = 280.d0
#endif

!     temp storage for new tracer interfaces
      integer, allocatable :: values(:)
      integer :: val
      integer :: nn

!     call routine to read/set up regions and sectors for emissions:
      call setup_emis_sectors_regions()
      call initializeOldTracers(tracers, setDefaultSpec)

#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP) ||\
      (defined TRACERS_TOMAS)
C**** DMS, seasalt from offline fields
      call sync_param("OFFLINE_DMS_SS",OFFLINE_DMS_SS)
C**** seasalt from offline fields
      call sync_param("OFFLINE_SS",OFFLINE_SS)
#endif

      if (tracers_special_shindell) then
#ifdef TRACERS_SPECIAL_Shindell
        call  Ox_setSpec('Ox')
        call  NOx_setSpec('NOx')
        call  ClOx_setSpec('ClOx')
        call  BrOx_setSpec('BrOx')
        call  N2O5_setSpec('N2O5')
        call  HNO3_setSpec('HNO3')
        call  H2O2_setSpec('H2O2')
        call  CH3OOH_setSpec('CH3OOH')

        call  HCHO_setSpec('HCHO')
        call  HO2NO2_setSpec('HO2NO2')
        call  CO_setSpec('CO')
        call  CH4_setSpec('CH4')
        call  PAN_setSpec('PAN')
        call  Isoprene_setSpec('Isoprene')
        call  AlkylNit_setSpec('AlkylNit')
        call  Alkenes_setSpec('Alkenes')
        call  Paraffin_setSpec('Paraffin')

        if (tracers_terp) then
          call  Terpenes_setSpec('Terpenes')
        end if

#ifdef TRACERS_AEROSOLS_SOA
        if (tracers_aerosols_soa) then
          call  isopp1g_setSpec('isopp1g')
          call  isopp1a_setSpec('isopp1a')
          call  isopp2g_setSpec('isopp2g')
          call  isopp2a_setSpec('isopp2a')
          if (tracers_terp) then
            call  apinp1g_setSpec('apinp1g')
            call  apinp1a_setSpec('apinp1a')
            call  apinp2g_setSpec('apinp2g')
            call  apinp2a_setSpec('apinp2a')
          end if
        end if
#endif
        call  HCl_setSpec('HCl')
        call  HOCl_setSpec('HOCl')
        call  ClONO2_setSpec('ClONO2')
        call  HBr_setSpec('HBr')
        call  HOBr_setSpec('HOBr')
        call  BrONO2_setSpec('BrONO2')
        call  N2O_setSpec('N2O')
        call  CFC_setSpec('CFC')

        if (shindell_strat_extra) then
          if (accmip_like_diags) then
            call  codirect_setSpec('codirect')
            call  stratOx_setSpec('stratOx')
            call  GLT_setSpec('GLT') ! generic linear tracer
          end if
        end if

#endif
      end if

      if (tracers_special_lerner) then
        call  N2O_setSpec('N2O') ! warning duplicate with Shindell
        call  CH4_setSpec('CH4') ! warning duplicate with Shindell
      end if

      if ((.not. tracers_amp) .and. tracers_water) then
        call  Water_setSpec('Water')
      end if
      

#ifdef TRACERS_SPECIAL_O18
      call  H2O18_setSpec('H2O18')
      call  HDO_setSpec('HDO')
c$$$      call  HTO_setSpec('HTO')
c$$$      call  H2O17_setSpec('H2O17')
#endif  /* TRACERS_SPECIAL_O18 */

      if (tracers_gasexch_ocean_cfc) then
        call  CFCn_setSpec('CFCn')
      end if

      if (tracers_gasexch_ocean_co2 .or. tracers_gasexch_land_co2) then
        call  CO2n_setSpec('CO2n')
      end if
      
      if (tracers_special_lerner) then
        call  SF6_setSpec('SF6')
        call  Rn222_setSpec('Rn222')
        call  CO2_setSpec('CO2')
        call  CFC11_setSpec('CFC11')
        call  C_14O2_setSpec('14CO2')
        call  O3_setSpec('O3')
        call  SF6_c_setSpec('SF6_c')
        if (tracers_special_shindell) 
     &       call stop_model('contradictory tracer specs')
      end if

      if (tracers_aerosols_koch) then
        call  DMS_setSpec('DMS')
        call  MSA_setSpec('MSA')
        call  SO2_setSpec('SO2')
        call  SO4_setSpec('SO4')
#ifndef TRACERS_SPECIAL_Shindell
        call  H2O2_s_setSpec('H2O2_s')
#endif  /* TRACERS_SPECIAL_Shindell */
        if (.not. sulf_only_aerosols) then
          call  seasalt1_setSpec('seasalt1')
          call  seasalt2_setSpec('seasalt2')
          call  BCII_setSpec('BCII')
          call  BCIA_setSpec('BCIA')
          call  BCB_setSpec('BCB')
        end if
#ifdef TRACERS_AEROSOLS_VBS
          call  VBS_setSpec('vbsGm2', 1,'igas')
          call  VBS_setSpec('vbsGm1', 2,'igas')
          call  VBS_setSpec('vbsGz', 3,'igas')
          call  VBS_setSpec('vbsGp1', 4,'igas')
          call  VBS_setSpec('vbsGp2', 5,'igas')
          call  VBS_setSpec('vbsGp3', 6,'igas')
          call  VBS_setSpec('vbsGp4', 7,'igas')
          call  VBS_setSpec('vbsGp5', 8,'igas')
          call  VBS_setSpec('vbsGp6', 9,'igas')

          call  VBS_setSpec('vbsAm2', 1,'iaer')
          call  VBS_setSpec('vbsAm1', 2,'iaer')
          call  VBS_setSpec('vbsAz', 3,'iaer')
          call  VBS_setSpec('vbsAp1', 4,'iaer')
          call  VBS_setSpec('vbsAp2', 5,'iaer')
          call  VBS_setSpec('vbsAp3', 6,'iaer')
          call  VBS_setSpec('vbsAp4', 7,'iaer')
          call  VBS_setSpec('vbsAp5', 8,'iaer')
          call  VBS_setSpec('vbsAp6', 9,'iaer')
#else
          call  OCII_setSpec('OCII')   !Insoluble industrial organic mass
          call  OCIA_setSpec('OCIA')   !Aged industrial organic mass
          call  OCB_setSpec('OCB')     !Biomass organic mass
#endif /* TRACERS_AEROSOLS_VBS */
      end if

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

      if (tracers_nitrate) then
        call  NH3_setSpec('NH3')
        call  NH4_setSpec('NH4')
        call  NO3p_setSpec('NO3p')
      end if

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
        if (tracers_radon) call  Pb210_setSpec('Pb210')
        call  Be7_setSpec('Be7')
        call  Be10_setSpec('Be10')
        if (tracers_radon) call  Rn222_setSpec('Rn222') ! duplicate with Lerner
        if (tracers_special_lerner) 
     &       call stop_model('contradictory tracer specs')
      end if

#ifdef TRACERS_MINERALS
      if (tracers_minerals) then
        call  Clayilli_setSpec('Clayilli') ! http://webmineral.com/data/Illite.shtml
        call  Claykaol_setSpec('Claykaol') ! http://www.webmineral.com/data/Kaolinite.shtml
        call  Claysmec_setSpec('Claysmec') ! http://www.webmineral.com/data/Rectorite.shtml
        call  Claycalc_setSpec('Claycalc') ! http://www.webmineral.com/data/Calcite.shtml
        call  Clayquar_setSpec('Clayquar') ! http://www.webmineral.com/data/Quartz.shtml
        call  Sil1quar_setSpec('Sil1quar') ! http://www.webmineral.com/data/Quartz.shtml
        call  Sil1feld_setSpec('Sil1feld') ! http://www.mindat.org/min-1624.html
        call  Sil1calc_setSpec('Sil1calc') ! http://www.webmineral.com/data/Calcite.shtml
        call  Sil1hema_setSpec('Sil1hema') ! http://www.webmineral.com/data/Hematite.shtml
        call  Sil1gyps_setSpec('Sil1gyps') ! http://www.webmineral.com/data/Gypsum.shtml
        call  Sil2quar_setSpec('Sil2quar') ! http://www.webmineral.com/data/Quartz.shtml
        call  Sil2feld_setSpec('Sil2feld') ! http://www.mindat.org/min-1624.html
        call  Sil2calc_setSpec('Sil2calc') ! http://www.webmineral.com/data/Calcite.shtml
        call  Sil2hema_setSpec('Sil2hema') ! http://www.webmineral.com/data/Hematite.shtml
        call  Sil2gyps_setSpec('Sil2gyps') ! http://www.webmineral.com/data/Gypsum.shtml
        call  Sil3quar_setSpec('Sil3quar') ! http://www.webmineral.com/data/Quartz.shtml
        call  Sil3feld_setSpec('Sil3feld') ! http://www.mindat.org/min-1624.html
        call  Sil3calc_setSpec('Sil3calc') ! http://www.webmineral.com/data/Calcite.shtml
        call  Sil3hema_setSpec('Sil3hema') ! http://www.webmineral.com/data/Hematite.shtml
        call  Sil3gyps_setSpec('Sil3gyps') ! http://www.webmineral.com/data/Gypsum.shtml
      end if
#endif  /* TRACERS_MINERALS */

#ifdef TRACERS_QUARZHEM
      if (tracers_quarzhem) then
        call  Sil1quhe_setSpec('Sil1quhe')
        call  Sil2quhe_setSpec('Sil2quhe')
        call  Sil3quhe_setSpec('Sil3quhe')
      end if
#endif

      if (tracers_air .or. htap_like_diags) then
        call  air_setSpec('air')
      end if

#ifdef TRACERS_AMP
      if (tracers_amp) then

C**** Tracers for Scheme AMP: Aerosol Microphysics (Mechanism M1 - M8)
        call  M_NO3_setSpec('M_NO3')
        call  M_NH4_setSpec('M_NH4')
        call  M_H2O_setSpec('M_H2O')

        if ( tracers_amp_m1 .or. tracers_amp_m2 .or. tracers_amp_m3
     &       .or. tracers_amp_m5 .or. tracers_amp_m6 .or. 
     &       tracers_amp_m7) then
          call  M_AKK_SU_setSpec('M_AKK_SU')
          call  N_AKK_1_setSpec('N_AKK_1')
        end if

        call  M_ACC_SU_setSpec('M_ACC_SU')
        call  N_ACC_1_setSpec('N_ACC_1')
        call  M_DD1_SU_setSpec('M_DD1_SU')
        call  M_DD1_DU_setSpec('M_DD1_DU')
        call  N_DD1_1_setSpec('N_DD1_1')
        call  M_DS1_SU_setSpec('M_DS1_SU')
        call  M_DS1_DU_setSpec('M_DS1_DU')
        call  N_DS1_1_setSpec('N_DS1_1')

        if (tracers_amp_m1 .or. tracers_amp_m2 .or. tracers_amp_m3
     &       .or. tracers_amp_m4) then
          call  M_DD2_SU_setSpec('M_DD2_SU')
          call  M_DD2_DU_setSpec('M_DD2_DU')
          call  N_DD2_1_setSpec('N_DD2_1')
          call  M_DS2_SU_setSpec('M_DS2_SU')
          call  M_DS2_DU_setSpec('M_DS2_DU')
          call  N_DS2_1_setSpec('N_DS2_1')
        end if

        if (tracers_amp_m1 .or.  tracers_amp_m2 .or. tracers_amp_m3
     &       .or. tracers_amp_m5 .or.  tracers_amp_m6 .or. 
     &       tracers_amp_m7) then
          call  M_SSA_SU_setSpec('M_SSA_SU')
          call  M_SSA_SS_setSpec('M_SSA_SS')
          call  M_SSC_SS_setSpec('M_SSC_SS')
        end if

        if ( tracers_amp_m4 .or. tracers_amp_m8) then
          call  M_SSS_SU_setSpec('M_SSS_SU')
          call  M_SSS_SS_setSpec('M_SSS_SS')
        end if

        call  M_OCC_SU_setSpec('M_OCC_SU')
        call  M_OCC_OC_setSpec('M_OCC_OC')
        call  N_OCC_1_setSpec('N_OCC_1')
        call  M_BC1_SU_setSpec('M_BC1_SU')
        call  M_BC1_BC_setSpec('M_BC1_BC')
        call  N_BC1_1_setSpec('N_BC1_1')
        call  M_BC2_SU_setSpec('M_BC2_SU')
        call  M_BC2_BC_setSpec('M_BC2_BC')
        call  N_BC2_1_setSpec('N_BC2_1')

        if ( tracers_amp_m1 .or. tracers_amp_m5) then
          call  M_BC3_SU_setSpec('M_BC3_SU')
          call  M_BC3_BC_setSpec('M_BC3_BC')
          call  N_BC3_1_setSpec('N_BC3_1')
        end if

        if ( tracers_amp_m2 .or. tracers_amp_m6) then
          call  M_OCS_SU_setSpec('M_OCS_SU')
          call  M_OCS_OC_setSpec('M_OCS_OC')
          call  N_OCS_1_setSpec('N_OCS_1')
        end if
          
        if ( tracers_amp_m1 .or. tracers_amp_m2 .or. tracers_amp_m6)then
          call  M_DBC_SU_setSpec('M_DBC_SU')
          call  M_DBC_BC_setSpec('M_DBC_BC')
          call  M_DBC_DU_setSpec('M_DBC_DU')
          call  N_DBC_1_setSpec('N_DBC_1')
        end if

        if (tracers_amp_m1 .or. tracers_amp_m2 .or. tracers_amp_m3
     &       .or. tracers_amp_m6 .or. tracers_amp_m7) then
          call  M_BOC_SU_setSpec('M_BOC_SU')
          call  M_BOC_BC_setSpec('M_BOC_BC')
          call  M_BOC_OC_setSpec('M_BOC_OC')
          call  N_BOC_1_setSpec('N_BOC_1')
        end if

        if ( tracers_amp_m1 .or. tracers_amp_m2 .or. tracers_amp_m5
     &       .or. TRACERS_AMP_M6) then
          call  M_BCS_SU_setSpec('M_BCS_SU')
          call  M_BCS_BC_setSpec('M_BCS_BC')
          call  N_BCS_1_setSpec('N_BCS_1')
        end if

        call  M_MXX_SU_setSpec('M_MXX_SU')
        call  M_MXX_BC_setSpec('M_MXX_BC')
        call  M_MXX_OC_setSpec('M_MXX_OC')
        call  M_MXX_DU_setSpec('M_MXX_DU')
        call  M_MXX_SS_setSpec('M_MXX_SS')
        call  N_MXX_1_setSpec('N_MXX_1')

        call  H2SO4_setSpec('H2SO4')
        call  DMS_setSpec('DMS') ! duplicate with Koch
        call  SO2_setSpec('SO2') ! duplicate with Koch
#ifndef TRACERS_SPECIAL_Shindell
        call  H2O2_s_setSpec('H2O2_s') ! duplicate with Koch
#endif  /* TRACERS_SPECIAL_Shindell */
        call  NH3_setSpec('NH3') ! duplicate with nitrate
        if (tracers_aerosols_koch) then
          call stop_model('contradictory tracer specs')
        end if
        if (tracers_nitrate) then
          call stop_model('contradictory tracer specs')
        end if
      end if
#endif /* TRACERS_AMP */

#ifdef TRACERS_TOMAS
      CALL initbounds()
      call readbinact ("binact10_12.dat",binact10) 
      call readbinact ("binact02_12.dat",binact02) 
      call readfraction("fraction10_12.dat",fraction10) 
      call readfraction("fraction02_12.dat",fraction02) 
      call readmielut           ! aerosol radiation lookup table

!For aerosol tracers in TOMAS model, 
!fq_aer is determined based on kohler theory. 

      call  TOMAS_H2SO4_setSpec('H2SO4')
      call  DMS_setSpec('DMS') ! note duplicate with Koch
      call  SO2_setSpec('SO2')
#ifndef TRACERS_AEROSOLS_SOA
      call  TOMAS_SOAgas_setSpec('SOAgas')
#endif  /* TRACERS_AEROSOLS_SOA */
      call  H2O2_s_setSpec('H2O2_s') ! duplicate with Koch
      call  NH3_setSpec('NH3') ! duplicate with nitrate
      call  NH4_setSpec('NH4') ! duplicate with nitrate

! Koch SO4: 1700 kg/m3
! Koch BC : 1300 
! Koch OC : 1500
! Koch DUST : 2500 for  clay and 2650 for silt
! Koch SS : 2200 kg/m3

! so4 : 1780 kg/m3
! ss:  2165 kg/m3
!bc: 1800 kg/m3 or 2200 kg/m3
!oc:1400 kg/m3
!ddust : 2650 kg/m3 
      n_ASO4(:)  = TOMAS_setSpec(TOMAS_ASO4_setSpec,  'ASO4', nbins)
      n_ANACL(:) = TOMAS_setSpec(TOMAS_ANACL_setSpec, 'ANACL',nbins)
      n_AECOB(:) = TOMAS_setSpec(TOMAS_AECOB_setSpec, 'AECOB',nbins)
      n_AECIL(:) = TOMAS_setSpec(TOMAS_AECIL_setSpec, 'AECIL',nbins)
      n_AOCOB(:) = TOMAS_setSpec(TOMAS_AOCOB_setSpec, 'AOCOB',nbins)
      n_AOCIL(:) = TOMAS_setSpec(TOMAS_AOCIL_setSpec, 'AOCIL',nbins)
      n_ADUST(:) = TOMAS_setSpec(TOMAS_ADUST_setSpec, 'ADUST',nbins)
      n_ANUM(:) = TOMAS_setSpec(TOMAS_ANUM_setSpec,   'ANUM', nbins)
      n_AH2O(:)  = TOMAS_setSpec(TOMAS_AH2O_setSpec,  'AH2O', nbins)
#endif /* TRACERS_AMP */

      ! Generic tracer work
      ! All tracers must have been declared before reaching this point!!!
      ntm = getNumTracers(tracers)

      call set_param("NTM",NTM,'o')
      call set_param("TRNAME",trName(),ntm,'o')
      call printTracerNames(trName())

      ! Generic tracer work
      do n = 1, ntm
#ifdef TRACERS_WATER
C**** Tracers that are soluble or are scavenged or are water => wet dep
        if (tr_wd_type(n).eq.nWater.or.tr_wd_type(n) .EQ. nPART .or.
     *       tr_RKD(n).gt.0) then
          call set_dowetdep(n, .true.)
        end if
#endif
#ifdef TRACERS_DRYDEP
C**** If tracers are particles or have non-zero HSTAR or F0 do dry dep:
C**** Any tracers that dry deposits needs the surface concentration:
        if(HSTAR(n).GT.0..OR.F0(n).GT.0..OR.tr_wd_type(n).eq.nPART) then
          call set_dodrydep(n, .true.)
          call set_needtrs(n, .true.)
#ifdef TRACERS_WATER
          if (tr_wd_type(n).eq.nWATER) call stop_model
     &         ('A water tracer should not undergo dry deposition.',255)
#endif
        end if
#endif /* TRACERS_DRYDEP */

#ifdef TRACERS_ON
C**** Define the conversion from mass to volume units here
        call set_mass2vol(n, mair/tr_mm(n))
        call set_vol2mass(n, tr_mm(n)/mair)
        call set_to_conc(n, 0)
#ifdef TRACERS_SPECIAL_Shindell
C**** Aerosol tracer output should be mass mixing ratio
        select case (tr_wd_TYPE(n))
        case (nGAS)
          call set_to_volume_MixRat(n, 1) !gas output to volume mixing ratio
        case (nPART, nWATER)
          call set_to_volume_MixRat(n, 0) ! aerosol/water output to mass mixing ratio
        case default
          call set_to_volume_MixRat(n, 0) !default output to mass mixing ratio
        end select
#endif
#if defined(TRACERS_GASEXCH_ocean_CO2) || defined(TRACERS_GASEXCH_land_CO2)
        call set_to_volume_MixRat(n, 1) !gas output to volume mixing ratio
#endif
#endif /* TRACERS_ON */

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
    

#ifdef TRACERS_SPECIAL_Shindell
      subroutine Ox_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_Ox = n
      call set_ntm_power(n, -8)
      call set_tr_mm(n, 48.d0)
#ifdef TRACERS_DRYDEP
      call set_F0(n,  1.4d0)
      call set_HSTAR(n,  1.d-2)
#endif
      end subroutine Ox_setSpec

      subroutine NOx_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_NOx = n
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 14.01d0)
#ifdef TRACERS_DRYDEP
      call set_F0(n,  1.d-1)
      call set_HSTAR(n,  1.d-2)
#endif
#ifdef DYNAMIC_BIOMASS_BURNING
! 12 below are the 12 VDATA veg types or Ent remapped to them,
! from Olga Pechony's AR5_EPFC_factors_incl_SO2.xlsx file.
      emisPerFireByVegType(n,1:12)=(/0.0000000d+00, 1.1378230d-07,
     & 3.2166037d-07, 1.5559274d-07, 4.1611088d-07, 5.7316458d-07,
     & 2.1700112d-07, 3.0054335d-07, 0.0000000d+00, 0.0000000d+00,
     & 0.0000000d+00, 0.0000000d+00/)
#endif
#ifdef TRACERS_SPECIAL_Shindell
      call check_aircraft_sectors(n_NOx) ! special 3D source case
#endif
      end subroutine NOx_setSpec

      subroutine ClOx_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_ClOx = n
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 51.5d0)
      end subroutine ClOx_setSpec

      subroutine BrOx_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_BrOx = n
      call set_ntm_power(n, -14)
      call set_tr_mm(n, 95.9d0)
      end subroutine BrOx_setSpec

      subroutine N2O5_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_N2O5 = n
      call set_ntm_power(n, -12)
      call set_tr_mm(n, 108.02d0)
      end subroutine N2O5_setSpec

      subroutine HNO3_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_HNO3 = n
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 63.018d0)
      call set_tr_RKD(n, 2.073d3 ) ! in mole/J = 2.1d5 mole/(L atm)
#ifdef TRACERS_DRYDEP
      call set_HSTAR(n, 1.d14)
#endif
      end subroutine HNO3_setSpec

      subroutine H2O2_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_H2O2 = n
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 34.016d0)
      call set_tr_RKD(n, 9.869d2    ) ! in mole/J = 1.d5 mole/(L atm)
      call set_tr_DHD(n, -5.52288d4 ) ! in J/mole = -13.2 kcal/mole.
#ifdef TRACERS_DRYDEP
      call set_HSTAR(n, tr_RKD(n)*convert_HSTAR)
      call set_F0(n,  1.d0)
#endif
      end subroutine H2O2_setSpec

      subroutine CH3OOH_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_CH3OOH = n
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 48.042d0)
#ifdef TRACERS_DRYDEP
      call set_HSTAR(n,  3.d2)
#endif
      end subroutine CH3OOH_setSpec

      subroutine HCHO_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_HCHO = n
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 30.026d0)
      call set_tr_RKD(n, 6.218d1 ) ! mole/J = 6.3d3 mole/(L atm)
#ifdef TRACERS_DRYDEP
      call set_HSTAR(n, 6.d3)
#endif
      end subroutine HCHO_setSpec

      subroutine HO2NO2_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_HO2NO2 = n
      call set_ntm_power(n, -12)
      call set_tr_mm(n, 79.018d0)
      end subroutine HO2NO2_setSpec

      subroutine CO_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_CO = n
      call set_ntm_power(n, -8)
      call set_tr_mm(n, 28.01d0)
#ifdef DYNAMIC_BIOMASS_BURNING
! 12 below are the 12 VDATA veg types or Ent remapped to them,
! from Olga Pechony's AR5_EPFC_factors_incl_SO2.xlsx file.
      emisPerFireByVegType(n,1:12)=(/0.0000000d+00, 7.0401156d-06,
     & 1.8708386d-05, 1.0678024d-05, 2.6742857d-05, 4.0226296d-05,
     & 2.3661527d-05, 4.4639346d-05, 0.0000000d+00, 0.0000000d+00,
     & 0.0000000d+00, 0.0000000d+00/)
#endif
      end subroutine CO_setSpec
      subroutine PAN_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_PAN = n
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 121.054d0) ! assuming CH3COOONO2 = PAN)
#ifdef TRACERS_DRYDEP
      call set_HSTAR(n,  3.6d0)
#endif
      end subroutine PAN_setSpec

      subroutine Isoprene_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_Isoprene = n
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 60.05d0) ! i.e. 5 carbons
#ifdef TRACERS_DRYDEP
      call set_HSTAR(n,  1.3d-2)
#endif
      end subroutine Isoprene_setSpec

      subroutine AlkylNit_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_AlkylNit = n
      call set_ntm_power(n, -11)
      call set_tr_mm(n, mair)   !unknown molecular weight, so use air and make
                                ! note in the diagnostics write-out...
      end subroutine AlkylNit_setSpec

      subroutine Alkenes_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_Alkenes = n
      call set_ntm_power(n, -10)
      call set_tr_mm(n, 1.0d0)  ! So, careful: source files now in Kmole/m2/s or
                                ! equivalently, kg/m2/s for species with tr_mm=1
#ifdef DYNAMIC_BIOMASS_BURNING
! 12 below are the 12 VDATA veg types or Ent remapped to them,
! from Olga Pechony's AR5_EPFC_factors_incl_SO2.xlsx file.
      emisPerFireByVegType(n,1:12)=(/0.0000000d+00, 6.1516259d-09,
     & 1.1544214d-08, 6.9501711d-09, 1.7481154d-08, 2.5840087d-08,
     & 1.5709551d-08, 3.5913079d-08, 0.0000000d+00, 0.0000000d+00,
     & 0.0000000d+00, 0.0000000d+00/)
#endif
      end subroutine Alkenes_setSpec

      subroutine Paraffin_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_Paraffin = n
      call set_ntm_power(n, -10)
      call set_tr_mm(n, 1.0d0)  ! So, careful: source files now in Kmole/m2/s or
                                ! equivalently, kg/m2/s for species with tr_mm=1
#ifdef DYNAMIC_BIOMASS_BURNING
! 12 below are the 12 VDATA veg types or Ent remapped to them,
! from Olga Pechony's AR5_EPFC_factors_incl_SO2.xlsx file.
      emisPerFireByVegType(n,1:12)=(/0.0000000d+00, 1.5258348d-09,
     & 5.6236904d-09, 3.1752858d-09, 1.0662656d-08, 1.5271524d-08,
     & 8.0735774d-09, 2.6055675d-08, 0.0000000d+00, 0.0000000d+00,
     & 0.0000000d+00, 0.0000000d+00/)
#endif
      end subroutine Paraffin_setSpec

      subroutine Terpenes_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_Terpenes = n
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 120.10d0) ! i.e. 10 carbons
#ifdef TRACERS_DRYDEP
      call set_HSTAR(n,  1.3d-2)
#endif
      end subroutine Terpenes_setSpec

#ifdef TRACERS_AEROSOLS_SOA
      subroutine isopp1g_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_isopp1g = n
      n_soa_i = n_isopp1g       !the first from the soa species
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 15.6d0)
      call set_tr_RKD(n, 1.d4 / convert_HSTAR ) !Henry; from mole/(L atm) to mole/J
      call set_tr_DHD(n, -12.d0 * gasc        ) !Henry temp dependence (J/mole), Chung and Seinfeld, 2002
      call set_tr_wd_type(n, ngas)
#ifdef TRACERS_DRYDEP
      call set_HSTAR(n,  tr_RKD(n) * convert_HSTAR)
#endif
      end subroutine isopp1g_setSpec

      subroutine isopp1a_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_isopp1a = n
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 15.6d0)
      call set_trpdens(n, 1.5d3) !kg/m3
      call set_trradius(n, 3.d-7) !m
      call set_fq_aer(n, 0.8d0) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, nPART)
      end subroutine isopp1a_setSpec

      subroutine isopp2g_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_isopp2g = n
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 15.6d0)
      call set_tr_RKD(n, 1.d4 / convert_HSTAR ) !Henry; from mole/(L atm) to mole/J
      call set_tr_DHD(n, -12.d0 * gasc        ) !Henry temp dependence (J/mole), Chung and Seinfeld, 2002
      call set_tr_wd_type(n, ngas)
#ifdef TRACERS_DRYDEP
      call set_HSTAR(n,  tr_RKD(n) * convert_HSTAR)
#endif
      end subroutine isopp2g_setSpec

      subroutine isopp2a_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_isopp2a = n
#ifndef TRACERS_TERP
      n_soa_e = n_isopp2a       !the last from the soa species
#endif  /* TRACERS_TERP */
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 15.6d0)
      call set_trpdens(n, 1.5d3) !kg/m3
      call set_trradius(n, 3.d-7) !m
      call set_fq_aer(n, 0.8d0) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, nPART)
      end subroutine isopp2a_setSpec

      subroutine apinp1g_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_apinp1g = n
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 15.6d0)
      call set_tr_RKD(n, 1.d4 / convert_HSTAR ) !Henry; from mole/(L atm) to mole/J
      call set_tr_DHD(n, -12.d0 * gasc        ) !Henry temp dependence (J/mole), Chung and Seinfeld, 2002
      call set_tr_wd_type(n, ngas)
#ifdef TRACERS_DRYDEP
      call set_HSTAR(n,  tr_RKD(n) * convert_HSTAR)
#endif
      end subroutine apinp1g_setSpec

      subroutine apinp1a_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_apinp1a = n
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 15.6d0)
      call set_trpdens(n, 1.5d3) !kg/m3
      call set_trradius(n, 3.d-7) !m
      call set_fq_aer(n, 0.8d0) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, nPART)
      end subroutine apinp1a_setSpec

      subroutine apinp2g_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_apinp2g = n
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 15.6d0)
      call set_tr_RKD(n, 1.d4 / convert_HSTAR ) !Henry; from mole/(L atm) to mole/J
      call set_tr_DHD(n, -12.d0 * gasc        ) !Henry temp dependence (J/mole), Chung and Seinfeld, 2002
      call set_tr_wd_type(n, ngas)
#ifdef TRACERS_DRYDEP
      call set_HSTAR(n,  tr_RKD(n) * convert_HSTAR)
#endif
      end subroutine apinp2g_setSpec

      subroutine apinp2a_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_apinp2a = n
      n_soa_e = n_apinp2a       !the last from the soa species
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 15.6d0)
      call set_trpdens(n, 1.5d3) !kg/m3
      call set_trradius(n, 3.d-7) !m
      call set_fq_aer(n, 0.8d0) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, nPART)
      end subroutine apinp2a_setSpec
#endif  /* TRACERS_AEROSOLS_SOA */

      subroutine HCl_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_HCl = n
      call set_ntm_power(n, -10)
      call set_tr_mm(n, 36.5d0)
      end subroutine HCl_setSpec

      subroutine HOCl_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_HOCl = n
      call set_ntm_power(n, -12)
      call set_tr_mm(n, 52.5d0)
      end subroutine HOCl_setSpec

      subroutine ClONO2_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_ClONO2 = n
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 97.5d0)
      end subroutine ClONO2_setSpec

      subroutine HBr_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_HBr = n
      call set_ntm_power(n, -14)
      call set_tr_mm(n, 80.9d0)
      end subroutine HBr_setSpec
      
      subroutine HOBr_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_HOBr = n
      call set_ntm_power(n, -14)
      call set_tr_mm(n, 96.9d0)
      end subroutine HOBr_setSpec

      subroutine BrONO2_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_BrONO2 = n
      call set_ntm_power(n, -14)
      call set_tr_mm(n, 141.9d0)
      end subroutine BrONO2_setSpec

      subroutine CFC_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_CFC = n
      call set_ntm_power(n, -12)
      call set_tr_mm(n, 137.4d0) !CFC11
      end subroutine CFC_setSpec

      subroutine codirect_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_codirect = n
      call set_ntm_power(n, -8)
      call set_tr_mm(n, 28.01d0)
      call set_trdecay(n,  2.31482d-7) ! 1/(50 days)
! not a radiactive decay, but functionally identical
      end subroutine codirect_setSpec

      subroutine stratOx_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_stratOx = n
! assumes initial Ox conditions read in for Ox tracer
      call set_ntm_power(n, -8)
      call set_tr_mm(n, 48.d0)
#ifdef TRACERS_DRYDEP
      call set_F0(n,  1.0d0)
      call set_HSTAR(n,  1.d-2)
#endif
      end subroutine stratOx_setSpec

      subroutine GLT_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_GLT = n
      call set_ntm_power(n, -11)
      call set_tr_mm(n, mair)
      end subroutine GLT_setSpec

#endif /* TRACERS_SPECIAL_Shindell */

      subroutine CH4_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_CH4 = n
      call set_tr_mm(n, 16.d0)
#ifdef TRACERS_SPECIAL_Lerner
      call set_ntsurfsrc(n,  14)
      call set_ntm_power(n, -9)
      n_MPtable(n) = 3
      tcscale(n_MPtable(n)) = 1.
#endif
#ifdef TRACERS_SPECIAL_Shindell
      call set_ntm_power(n, -8)
#ifdef DYNAMIC_BIOMASS_BURNING
!     12 below are the 12 VDATA veg types or Ent remapped to them,
!     from Olga Pechony's AR5_EPFC_factors_incl_SO2.xlsx file.
      emisPerFireByVegType(n,1:12)=(/0.0000000d+00,1.0864168d-06,
     & 6.3624935d-07, 4.7021388d-07, 1.0293500d-06, 1.7132404d-06,
     & 1.4364367d-06, 3.0849296d-06, 0.0000000d+00, 0.0000000d+00,
     & 0.0000000d+00, 0.0000000d+00/)
#endif /* DYNAMIC_BIOMASS_BURNING */
#endif /* TRACERS_SPECIAL_Shindell */
      end subroutine CH4_setSpec

      subroutine N2O_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_N2O = n
      call set_ntm_power(n, -9)
      call set_tr_mm(n, 44.d0)
#ifdef TRACERS_SPECIAL_Lerner
      call set_ntsurfsrc(n,  1)
      n_MPtable(n) = 1
      tcscale(n_MPtable(n)) = 1.
#endif
      end subroutine N2O_setSpec

      subroutine Water_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_Water = n
#if (defined TRACERS_WATER) || (defined TRACERS_OCEAN)
      call set_trw0(n, 1.d+0)
      call set_ntrocn(n, 0)
#endif
#ifdef TRACERS_WATER
      call set_ntm_power(n, -4)
      call set_tr_mm(n, mwat)
      call set_needtrs(n,  .true.)
      call set_tr_wd_type(n, nWater)
      call set_trli0(n, 1.d+0)
      call set_trsi0(n, 1.d+0)
      call set_tr_H2ObyCH4(n, 1.d+0)
#endif
#ifdef TRACERS_OCEAN
      call set_trglac(n, 1.d+0)
#endif
#ifdef TRACERS_SPECIAL_O18
      iso_index(n) = 1          ! indexing for isotopic fractionation calcs
#endif
      end subroutine Water_setSpec

#ifdef TRACERS_SPECIAL_O18
      subroutine H2O18_setSpec(name)
      use TracerConstants_mod, only: H2O18
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
#ifdef TRACERS_OCEAN
      call set_trglac(n, trw0(n)*0.98d0   ) ! d=-20
#endif
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
#ifdef TRACERS_OCEAN
      call set_trglac(n, trw0(n)*0.84d0   ) ! d=-160
#endif
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
#ifdef TRACERS_OCEAN
      call set_trglac(n, 0d0)
#endif
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
#ifdef TRACERS_OCEAN
      call set_trglac(n, trw0(n)*0.98937d0   ) ! d=-10.63 D17O=
#endif
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
      use TRACER_COM, only: set_ntsurfsrc
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_CO2n = n
      call set_ntm_power(n, -6)
      call set_tr_mm(n, 44.d0)  !grams
      call set_t_qlimit(n,  .false.)
      call set_ntsurfsrc(n,  1)
      call set_needtrs(n, .true.)
      end subroutine CO2n_setSpec

      subroutine SF6_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_SF6 = n
      call set_ntm_power(n, -14)
      call set_tr_mm(n, 146.01d0)
      call set_ntsurfsrc(n,  1)
      end subroutine SF6_setSpec

      subroutine CO2_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_CO2 = n
      call set_ntm_power(n, -6)
      call set_tr_mm(n, 44.d0)
      call set_t_qlimit(n,  .false.)
      call set_ntsurfsrc(n,  6)
      end subroutine CO2_setSpec

      subroutine CFC11_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_CFC11 = n
      call set_ntm_power(n, -12)
      call set_tr_mm(n, 137.4d0)
      call set_ntsurfsrc(n,  1)
#ifdef TRACERS_SPECIAL_Lerner
      n_mptable(n) = 2
      tcscale(n_MPtable(n)) = 1.
#endif
      end subroutine CFC11_setSpec

      subroutine C_14O2_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_14CO2 = n
      call set_ntm_power(n, -18)
      call set_tr_mm(n, 46.d0)
      call set_ntsurfsrc(n,  1)
      end subroutine C_14O2_setSpec

      subroutine O3_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_O3 = n
      call set_ntm_power(n, -8)
      call set_tr_mm(n, 48.d0)
      call set_ntsurfsrc(n,  1)
#ifdef TRACERS_SPECIAL_Lerner
C**** Get solar variability coefficient from namelist if it exits
      dsol = 0.
      call sync_param("dsol",dsol)
#endif
      end subroutine O3_setSpec

      subroutine SF6_c_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_SF6_c = n
      call set_ntm_power(n, -14)
      call set_tr_mm(n, 146.01d0)
      call set_ntsurfsrc(n,  1)
      end subroutine SF6_c_setSpec

      subroutine DMS_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_DMS = n
      call set_ntm_power(n, -12)
!     the ocean source of DMS is actually interactive and therefore should
!     not count for ntsurfsrc....
      call set_ntsurfsrc(n,  0) ! ocean DMS concentration
      call set_ntisurfsrc(n, 1)
      call set_tr_mm(n, 62.d+0)
      call set_needtrs(n, .true.)
      end subroutine DMS_setSpec

      subroutine MSA_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_MSA = n
      call set_ntm_power(n, -13)
      call set_tr_mm(n, 96.d+0) !(H2O2 34;SO2 64)
      call set_trpdens(n, 1.7d3) !kg/m3 this is sulfate value
      call set_trradius(n, 5.d-7 ) !m (SO4 3;BC 1;OC 3)
      call set_fq_aer(n, 1.0d0   ) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, npart)
      end subroutine MSA_setSpec

      subroutine SO2_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_SO2 = n
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 64.d+0)
      call set_tr_RKD(n, 0.0118d0 ) !mole/J or  1.2  M/atm
      call set_tr_DHD(n, -2.62d4) ! in J/mole= -6.27 kcal/mol
      call set_tr_wd_type(n, ngas)
#ifdef TRACERS_DRYDEP
c     call set_HSTAR(n, tr_RKD(n)*convert_HSTAR)
      CALL SET_HSTAR(N, 1.D5)
#endif
#ifdef DYNAMIC_BIOMASS_BURNING
! 12 below are the 12 VDATA veg types or Ent remapped to them,
! from Olga Pechony's AR5_EPFC_factors_incl_SO2.xlsx file.
      emisPerFireByVegType(n,1:12)=(/0.0000000d+00, 2.6305068d-07,
     &     1.3513656d-07, 1.0093965d-07, 1.6911058d-07, 3.2019645d-07,
     &     3.1232341d-07, 4.1607765d-07, 0.0000000d+00, 0.0000000d+00,
     &     0.0000000d+00, 0.0000000d+00/)
#endif
      end subroutine SO2_setSpec

      subroutine SO4_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_SO4 = n
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 96.d+0)
      call set_trpdens(n, 1.7d3) !kg/m3 this is sulfate value
      call set_trradius(n, 3.d-7 ) !m
      call set_fq_aer(n, 1.d0   ) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, npart)
      end subroutine SO4_setSpec

      subroutine H2O2_s_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_H2O2_s = n
      call set_ntm_power(n, -10)
      call set_tr_mm(n, 34.016d0)
      call set_tr_RKD(n, 986.9d0)
      call set_tr_DHD(n, -5.52288d4 ) ! in J/mole = -13.2 kcal/mole.
      call set_tr_wd_type(n, ngas)
#ifdef TRACERS_DRYDEP
      call set_HSTAR(n, tr_RKD(n)*convert_HSTAR)
      call set_F0(n,  1.d0)
#endif
      end subroutine H2O2_s_setSpec

      subroutine seasalt1_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_seasalt1 = n
      call set_ntsurfsrc(n,  0) ! ocean bubbles
      call set_ntisurfsrc(n, 1)
      call set_ntm_power(n, -10)
      call set_tr_mm(n, 75.d0)  !Na x 3.256
      call set_trpdens(n, 2.2d3) !kg/m3 This is for non-hydrated
      call set_trradius(n, 4.4d-7 ) ! This is non-hydrated
      call set_fq_aer(n, 1.0d0   ) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, npart)
      end subroutine seasalt1_setSpec

      subroutine seasalt2_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_seasalt2 = n
      call set_ntsurfsrc(n,  0) ! ocean bubbles
      call set_ntisurfsrc(n, 1)
      call set_ntm_power(n, -9)
      call set_tr_mm(n, 75.d0)  !Na x 3.256
      call set_trpdens(n, 2.2d3) !kg/m3 This is for non-hydrated
      call set_trradius(n, 5.0d-6) ! This is non-hydrated
#ifdef TRACERS_AEROSOLS_Koch
      if (OFFLINE_DMS_SS.ne.1 .and. OFFLINE_SS.ne.1) then
        call set_trradius(n, 1.7d-6 ) ! This is non-hydrated
      end if
#endif
      call set_fq_aer(n, 1.0d0   ) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, npart)
      end subroutine seasalt2_setSpec

      subroutine BCII_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_BCII = n
      call set_ntm_power(n, -12)
      call set_tr_mm(n, 12.d0)
      call set_trpdens(n, 1.3d3) !kg/m3
      call set_trradius(n, 1.d-7 ) !m
      call set_fq_aer(n, 0.0d0   ) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, npart)
      end subroutine BCII_setSpec

      subroutine BCIA_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_BCIA = n
      call set_ntm_power(n, -12)
      call set_tr_mm(n, 12.d0)
      call set_trpdens(n, 1.3d3) !kg/m3
      call set_trradius(n, 1.d-7 ) !m
      call set_fq_aer(n, 1.d0   ) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, npart)
      end subroutine BCIA_setSpec

      subroutine BCB_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_BCB = n
      call set_ntm_power(n, -12)
      call set_tr_mm(n, 12.d0)
      call set_trpdens(n, 1.3d3) !kg/m3
      call set_trradius(n, 1.d-7 ) !m
      call set_fq_aer(n, 0.6d0 ) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, npart)
#ifdef DYNAMIC_BIOMASS_BURNING
! 12 below are the 12 VDATA veg types or Ent remapped to them,
! from Olga Pechony's AR5_EPFC_factors_incl_SO2.xlsx file.
      emisPerFireByVegType(n,1:12)=(/0.0000000d+00, 8.2731990d-09,
     & 1.7817767d-07, 8.5456378d-08, 2.5662467d-07, 3.3909114d-07,
     & 1.1377826d-07, 2.9145593d-07, 0.0000000d+00, 0.0000000d+00,
     & 0.0000000d+00, 0.0000000d+00/)
#endif
      end subroutine BCB_setSpec

#ifdef TRACERS_AEROSOLS_VBS
      subroutine VBS_setSpec(name, n, index, type) result(label)
      character(len=*), intent(in) :: name
      integer, intent(in) :: index
      character(len=4), intent(in) :: type

      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      label = n

      select case (type)
      case ('igas')
        vbs_tr%igas(index) = n
        call num_srf_sources(n,.false.)
      case ('iaer')
        vbs_tr%iaer(index) = n
        call num_srf_sources(n,.true.)
      end selec t case

      call set_ntm_power(n, -11)
      call set_tr_mm(n, 15.6d0)
      select case(name)
      case ('vbsGm2', 'vbsGm1', 'vbsGz',  'vbsGp1', 'vbsGp2', ! VBS gas-phase
     &      'vbsGp3', 'vbsGp4', 'vbsGp5', 'vbsGp6')
        call set_tr_wd_type(n, ngas)
        call set_tr_RKD(n, 1.d4 / convert_HSTAR ) !Henry; from mole/(L atm) to mole/J
!     call set_tr_DHD(n, -12.d0 * gasc        ) !Henry temp dependence (J/mole), Chung and Seinfeld, 2002
#ifdef TRACERS_DRYDEP
        call set_HSTAR(n,  tr_RKD(n) * convert_HSTAR)
#endif
      case ('vbsAm2', 'vbsAm1', 'vbsAz',  'vbsAp1', 'vbsAp2', ! VBS aerosol-phase
     &      'vbsAp3', 'vbsAp4', 'vbsAp5', 'vbsAp6')
#ifdef DYNAMIC_BIOMASS_BURNING
!     12 below are the 12 VDATA veg types or Ent remapped to them,
!     from Olga Pechony's AR5_EPFC_factors_incl_SO2.xlsx file.
        emisPerFireByVegType(n,1:12)=(/0.0000000d+00, 6.4230818d-07,
     &   1.6844633d-06, 8.7537586d-07, 2.5902559d-06, 4.3200689d-06,
     &   2.6824284d-06, 3.9549395d-06, 0.0000000d+00, 0.0000000d+00,
     &   0.0000000d+00, 0.0000000d+00/)
#endif
        call set_tr_wd_type(n, npart)
        call set_trpdens(n, 1.5d3) !kg/m3
        call set_trradius(n, 3.d-7 ) !m
        call set_fq_aer(n, 1.0d0   ) !fraction of aerosol that dissolves
      end select
      end subroutine VBS_setSpec
#endif /* TRACERS_AEROSOLS_VBS */

      subroutine OCII_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_OCII = n
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 15.6d0)
      call set_trpdens(n, 1.5d3) !kg/m3
      call set_trradius(n, 3.d-7 ) !m
      call set_fq_aer(n, 0.0d0   ) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, npart)
      end subroutine OCII_setSpec

      subroutine OCIA_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_OCIA = n
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 15.6d0)
      call set_trpdens(n, 1.5d3) !kg/m3
      call set_trradius(n, 3.d-7 ) !m
      call set_fq_aer(n, 1.d0   ) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, npart)
      end subroutine OCIA_setSpec

      subroutine OCB_setSpec(name)
      use OldTracer_mod, only: om2oc, set_om2oc
      character(len=*), intent(in) :: name
      real*8 :: tmp
      n = oldAddTracer(name)
      n_OCB = n
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)
      tmp = om2oc(n_OCB)
      call sync_param("OCB_om2oc",tmp)
      call set_om2oc(n_OCB, tmp)
#endif
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 15.6d0)
      call set_trpdens(n, 1.5d3) !kg/m3
      call set_trradius(n, 3.d-7 ) !m
      call set_fq_aer(n, 0.8d0   ) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, npart)
#ifdef DYNAMIC_BIOMASS_BURNING
! 12 below are the 12 VDATA veg types or Ent remapped to them,
! from Olga Pechony's AR5_EPFC_factors_incl_SO2.xlsx file.
      emisPerFireByVegType(n,1:12)=(/0.0000000d+00, 6.4230818d-07,
     &     1.6844633d-06, 8.7537586d-07, 2.5902559d-06, 4.3200689d-06,
     &     2.6824284d-06, 3.9549395d-06, 0.0000000d+00, 0.0000000d+00,
     &     0.0000000d+00, 0.0000000d+00/)
#endif
      end subroutine OCB_setSpec

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
#ifdef TRACERS_DRYDEP
      call set_trradius(n, 0.46D-06)
#endif
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
#ifdef TRACERS_DRYDEP
      call set_trradius(n, 1.47D-06)
#endif
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
#ifdef TRACERS_DRYDEP
      call set_trradius(n, 2.94D-06)
#endif
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
      end subroutine Silt2_setSpec

      subroutine Silt3_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_Silt3=n
      call set_ntm_power(n, -9)
      call set_trpdens(n, 2.65d3)
#ifdef TRACERS_DRYDEP
      call set_trradius(n, 5.88D-06)
#endif
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
#ifdef TRACERS_DRYDEP
      call set_trradius(n, 11.77D-06)
#endif
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
      end subroutine Silt4_setSpec

      subroutine NH3_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_NH3 = n
      call set_ntm_power(n, -10)
      call set_tr_mm(n, 17.d0)
      call set_tr_RKD(n, 0.7303d0   ) !tr_RKD=74 M/atm
      call set_tr_DHD(n, -2.84d4  ) !tr_DHD=-6.80 kcal/mole
      call set_tr_wd_type(n, ngas)
#ifdef TRACERS_DRYDEP
      call set_HSTAR(n, tr_RKD(n)*convert_HSTAR)
#endif
#ifdef DYNAMIC_BIOMASS_BURNING
! 12 below are the 12 VDATA veg types or Ent remapped to them,
! from Olga Pechony's AR5_EPFC_factors_corrected_NH3.xlsx file.
      emisPerFireByVegType(n,1:12)=(/0.0000000d+00, 1.2274993d-06,
     &     3.8813269d-07, 3.5230462d-07, 4.0781484d-07, 8.8901584d-07,
     &     1.1341459d-06, 1.4117913d-06, 0.0000000d+00, 0.0000000d+00,
     &     0.0000000d+00, 0.0000000d+00/)
#endif
      end subroutine NH3_setSpec

      subroutine NH4_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_NH4 = n
      call set_ntsurfsrc(n,  0)
      call set_ntm_power(n, -10)
      call set_tr_mm(n, 18.d0)
      call set_trpdens(n, 1.7d3)
      call set_trradius(n, 3.d-7)
      call set_fq_aer(n, 1.0d0   ) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, npart)
      end subroutine NH4_setSpec

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

      subroutine Rn222_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_Rn222 = n
      call set_ntm_power(n, -21)
      call set_tr_mm(n, 222.d0)
      call set_trdecay(n,  2.1d-6)
      call set_ntsurfsrc(n,  1)
      end subroutine Rn222_setSpec

#ifdef TRACERS_MINERALS
      subroutine Clayilli_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_Clayilli=n
          call set_ntm_power(n, -9)
          call set_trpdens(n, 2.795d3) ! measured; http://www.mindat.org/min-2011.html
#ifdef TRACERS_DRYDEP
      call set_trradius(n, 0.46D-06)
#endif
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)

      end subroutine Clayilli_setSpec

      subroutine Claykaol_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_Claykaol=n
      call set_ntm_power(n, -9)
      call set_trpdens(n, 2.63d3) ! calculated; http://www.mindat.org/min-2011.html
#ifdef TRACERS_DRYDEP
      call set_trradius(n, 0.46D-06)
#endif
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)

      end subroutine Claykaol_setSpec

      subroutine Claysmec_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_Claysmec=n
      call set_ntm_power(n, -9)
      call set_trpdens(n, 2.35D3) ! for Montmorillonite
#ifdef TRACERS_DRYDEP
      call set_trradius(n, 0.46D-06)
#endif
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)

      end subroutine Claysmec_setSpec

      subroutine Claycalc_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_Claycalc=n
          call set_ntm_power(n, -9)
          call set_trpdens(n, 2.71d3) ! measured; http://www.mindat.org/min-859.html
#ifdef TRACERS_DRYDEP
      call set_trradius(n, 0.46D-06)
#endif
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
      end subroutine Claycalc_setSpec

      subroutine Clayquar_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_Clayquar=n
      call set_ntm_power(n, -9)
      call set_trpdens(n, DensityQuartz)
#ifdef TRACERS_DRYDEP
      call set_trradius(n, 0.46D-06)
#endif
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
      end subroutine Clayquar_setSpec

      subroutine Sil1quar_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_Sil1quar=n
      call set_ntm_power(n, -9)
      call set_trpdens(n, DensityQuartz)
#ifdef TRACERS_DRYDEP
      call set_trradius(n, 1.47D-06)
#endif
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
      end subroutine Sil1quar_setSpec

      subroutine Sil1feld_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_Sil1feld=n
      call set_ntm_power(n, -9)
      call set_trpdens(n, 2.68d3) ! average Plagioclase
#ifdef TRACERS_DRYDEP
      call set_trradius(n, 1.47D-06)
#endif
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
      end subroutine Sil1feld_setSpec

      subroutine Sil1calc_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_Sil1calc=n
      call set_ntm_power(n, -9)
      call set_trpdens(n, 2.71d3) ! measured; http://www.mindat.org/min-859.html
#ifdef TRACERS_DRYDEP
      call set_trradius(n, 1.47D-06)
#endif
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
      end subroutine Sil1calc_setSpec

      subroutine Sil1hema_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_Sil1hema=n
      call set_ntm_power(n, -9)
      call set_trpdens(n, DensityHematite)
#ifdef TRACERS_DRYDEP
      call set_trradius(n, 1.47D-06)
#endif
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
      end subroutine Sil1hema_setSpec

      subroutine Sil1gyps_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_Sil1gyps=n
      call set_ntm_power(n, -9)
      call set_trpdens(n, 2.312d3) ! measured; http://www.mindat.org/min-1784.html
#ifdef TRACERS_DRYDEP
      call set_trradius(n, 1.47D-06)
#endif
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
      end subroutine Sil1gyps_setSpec

      subroutine Sil2quar_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_Sil2quar=n
      call set_ntm_power(n, -9)
      call set_trpdens(n, DensityQuartz)
#ifdef TRACERS_DRYDEP
      call set_trradius(n, 2.94D-06)
#endif
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
      end subroutine Sil2quar_setSpec

      subroutine Sil2feld_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_Sil2feld=n
      call set_ntm_power(n, -9)
      call set_trpdens(n, 2.68d3) ! average Plagioclase
#ifdef TRACERS_DRYDEP
      call set_trradius(n, 2.94D-06)
#endif
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
      end subroutine Sil2feld_setSpec

      subroutine Sil2calc_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_Sil2calc=n
      call set_ntm_power(n, -9)
      call set_trpdens(n, 2.71d3) ! measured; http://www.mindat.org/min-859.html
#ifdef TRACERS_DRYDEP
      call set_trradius(n, 2.94D-06)
#endif
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
      end subroutine Sil2calc_setSpec

      subroutine Sil2hema_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_Sil2hema=n
      call set_ntm_power(n, -9)
      call set_trpdens(n, DensityHematite)
#ifdef TRACERS_DRYDEP
      call set_trradius(n, 2.94D-06)
#endif
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
      end subroutine Sil2hema_setSpec

      subroutine Sil2gyps_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_Sil2gyps=n
      call set_ntm_power(n, -9)
      call set_trpdens(n, 2.312d3) ! measured; http://www.mindat.org/min-1784.html
#ifdef TRACERS_DRYDEP
      call set_trradius(n, 2.94D-06)
#endif
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
      end subroutine Sil2gyps_setSpec

      subroutine Sil3quar_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_Sil3quar=n
      call set_ntm_power(n, -9)
      call set_trpdens(n, DensityQuartz)
#ifdef TRACERS_DRYDEP
      call set_trradius(n, 5.88D-06)
#endif
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
      end subroutine Sil3quar_setSpec

      subroutine Sil3feld_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_Sil3feld=n
      call set_ntm_power(n, -9)
      call set_trpdens(n, 2.68d3) ! average Plagioclase
#ifdef TRACERS_DRYDEP
      call set_trradius(n, 5.88D-06)
#endif
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
      end subroutine Sil3feld_setSpec

      subroutine Sil3calc_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_Sil3calc=n
      call set_ntm_power(n, -9)
      call set_trpdens(n, 2.71d3) ! measured; http://www.mindat.org/min-859.html
#ifdef TRACERS_DRYDEP
      call set_trradius(n, 5.88D-06)
#endif
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
      end subroutine Sil3calc_setSpec

      subroutine Sil3hema_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_Sil3hema=n
      call set_ntm_power(n, -9)
      call set_trpdens(n, DensityHematite)
#ifdef TRACERS_DRYDEP
      call set_trradius(n, 5.88D-06)
#endif
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
      end subroutine Sil3hema_setSpec

      subroutine Sil3gyps_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_Sil3gyps=n
      call set_ntm_power(n, -9)
      call set_trpdens(n, 2.312d3) ! measured; http://www.mindat.org/min-1784.html
#ifdef TRACERS_DRYDEP
      call set_trradius(n, 5.88D-06)
#endif
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
      end subroutine Sil3gyps_setSpec

#endif  /* TRACERS_MINERALS */

#ifdef TRACERS_QUARZHEM
      subroutine Sil1quhe_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_Sil1quhe=n
      call set_ntm_power(n, -9)
      call set_trpdens(n, (1.d0 - frHemaInQuarAggr) * DensityQuartz
     &     + frHemaInQuarAggr * DensityHematite)
#ifdef TRACERS_DRYDEP
      call set_trradius(n, 1.47D-06)
#endif
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
      end subroutine Sil1quhe_setSpec

      subroutine Sil2quhe_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_Sil2quhe=n
      call set_ntm_power(n, -9)
      call set_trpdens(n, (1.d0 - frHemaInQuarAggr) * DensityQuartz
     &     + frHemaInQuarAggr * DensityHematite)
#ifdef TRACERS_DRYDEP
      call set_trradius(n, 2.94D-06)
#endif
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
      end subroutine Sil2quhe_setSpec

      subroutine Sil3quhe_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_Sil3quhe=n
      call set_ntm_power(n, -9)
      call set_trpdens(n, (1.d0 - frHemaInQuarAggr) * DensityQuartz
     &     + frHemaInQuarAggr * DensityHematite)
#ifdef TRACERS_DRYDEP
      call set_trradius(n, 5.88D-06)
#endif
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
      end subroutine Sil3quhe_setSpec
#endif  /* TRACERS_QUARZHEM */

      subroutine Air_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_Air = n
      call set_ntm_power(n, -2)
      call set_tr_mm(n, mair)
      end subroutine Air_setSpec


#ifdef TRACERS_AMP
      subroutine H2SO4_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_H2SO4 = n
          call set_ntm_power(n, -11)
          call set_ntsurfsrc(n,  0)
          call set_tr_mm(n, 98.d0)
          call set_trpdens(n, DENS_SULF)
          call set_trradius(n, DG_ACC * .5d-6)
          call set_fq_aer(n, SOLU_ACC)
          call set_tr_wd_type(n, npart)
      end subroutine H2SO4_setSpec

      subroutine M_NO3_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_M_NO3 = n
      ntmAMPi=n                 ! always the first tracer in AMP
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 62.d0)
      call set_trpdens(n, 1.7d3)
      call set_trradius(n, 3.d-7 ) !m
      call set_fq_aer(n, 1.d0)  !fraction of aerosol that dissolves
      call set_tr_wd_type(n, npart)
      end subroutine M_NO3_setSpec

      subroutine M_NH4_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_M_NH4 = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 18.d0)
      call set_trpdens(n, 1.7d3)
      call set_trradius(n, 3.d-7)
      call set_fq_aer(n, 1.d+0)
      call set_tr_wd_type(n, npart)
      end subroutine M_NH4_setSpec

      subroutine M_H2O_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_M_H2O = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, mwat)
      call set_trpdens(n, 1.d3)
      call set_trradius(n, 3.d-7)
      call set_fq_aer(n, 1.d+0)
      call set_tr_wd_type(n, npart) !nWater
      end subroutine M_H2O_setSpec

      subroutine M_AKK_SU_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_M_AKK_SU = n
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 96.d0)
      call set_trpdens(n, DENS_SULF)
      call set_trradius(n, DG_AKK * .5d-6)
      call set_fq_aer(n,  SOLU_AKK)
      call set_tr_wd_type(n, npart)
      end subroutine M_AKK_SU_setSpec

      subroutine N_AKK_1_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_N_AKK_1 = n
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 1.d+0)
      call set_trpdens(n, DENS_SULF)
      call set_trradius(n, DG_AKK * .5d-6)
      call set_fq_aer(n, SOLU_AKK)
      call set_tr_wd_type(n, npart)
      end subroutine N_AKK_1_setSpec

      subroutine M_ACC_SU_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_M_ACC_SU = n
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 96.d0)
      call set_trpdens(n, DENS_SULF)
      call set_trradius(n, DG_ACC * .5d-6)
      call set_fq_aer(n, SOLU_ACC)
      call set_tr_wd_type(n, npart)
      end subroutine M_ACC_SU_setSpec

      subroutine N_ACC_1_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_N_ACC_1 = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 1.d+0)
      call set_trpdens(n, DENS_SULF)
      call set_trradius(n, DG_ACC * .5d-6)
      call set_fq_aer(n, SOLU_ACC)
      call set_tr_wd_type(n, npart)
      end subroutine N_ACC_1_setSpec

      subroutine M_DD1_SU_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_M_DD1_SU = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 96.d0)
      call set_trpdens(n, DENS_SULF)
      call set_trradius(n, DG_DD1 * .5d-6)
      call set_fq_aer(n, SOLU_DD1)
      call set_tr_wd_type(n, npart)
      end subroutine M_DD1_SU_setSpec

      subroutine M_DD1_DU_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_M_DD1_DU = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 1.d+0)
      call set_trpdens(n, DENS_DUST)
      call set_trradius(n, DG_DD1 * .5d-6)
      call set_fq_aer(n,  SOLU_DD1)
      call set_tr_wd_type(n, npart)
      end subroutine M_DD1_DU_setSpec

      subroutine N_DD1_1_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_N_DD1_1 = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 1.d+0)
      call set_trpdens(n, DENS_DUST)
      call set_trradius(n, DG_DD1 * .5d-6)
      call set_fq_aer(n,  SOLU_DD1)
      call set_tr_wd_type(n, npart)
      end subroutine N_DD1_1_setSpec

      subroutine M_DS1_SU_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_M_DS1_SU = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 96.d0)
      call set_trpdens(n, DENS_SULF)
      call set_trradius(n, DG_DS1 * .5d-6)
      call set_fq_aer(n, SOLU_DS1)
      call set_tr_wd_type(n, npart)
      end subroutine M_DS1_SU_setSpec

      subroutine M_DS1_DU_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_M_DS1_DU = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 1.d+0)
      call set_trpdens(n, DENS_DUST)
      call set_trradius(n, DG_DS1 * .5d-6)
      call set_fq_aer(n, SOLU_DS1)
      call set_tr_wd_type(n, npart)
      end subroutine M_DS1_DU_setSpec

      subroutine N_DS1_1_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_N_DS1_1= n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 1.d+0)
      call set_trpdens(n, DENS_DUST)
      call set_trradius(n, DG_DS1 * .5d-6)
      call set_fq_aer(n, SOLU_DS1)
      call set_tr_wd_type(n, npart)
      end subroutine N_DS1_1_setSpec

      subroutine M_DD2_SU_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_M_DD2_SU = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 96.d0)
      call set_trpdens(n, DENS_SULF)
      call set_trradius(n, DG_DD2 * .5d-6)
      call set_fq_aer(n, SOLU_DD2)
      call set_tr_wd_type(n, npart)
      end subroutine M_DD2_SU_setSpec

      subroutine M_DD2_DU_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_M_DD2_DU = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 1.d+0)
      call set_trpdens(n, DENS_DUST)
      call set_trradius(n, DG_DD2 * .5d-6)
      call set_fq_aer(n,  SOLU_DD2)
      call set_tr_wd_type(n, npart)
      end subroutine M_DD2_DU_setSpec

      subroutine N_DD2_1_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_N_DD2_1 = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 1.d+0)
      call set_trpdens(n, DENS_DUST)
      call set_trradius(n, DG_DD2 * .5d-6)
      call set_fq_aer(n,  SOLU_DD2)
      call set_tr_wd_type(n, npart)
      end subroutine N_DD2_1_setSpec

      subroutine M_DS2_SU_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_M_DS2_SU = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 96.d0)
      call set_trpdens(n, DENS_SULF)
      call set_trradius(n, DG_DS2 * .5d-6)
      call set_fq_aer(n, SOLU_DS2)
      call set_tr_wd_type(n, npart)
      end subroutine M_DS2_SU_setSpec

      subroutine M_DS2_DU_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_M_DS2_DU = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 1.d+0)
      call set_trpdens(n, DENS_DUST)
      call set_trradius(n, DG_DS2 * .5d-6)
      call set_fq_aer(n, SOLU_DS2)
      call set_tr_wd_type(n, npart)
      end subroutine M_DS2_DU_setSpec

      subroutine N_DS2_1_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_N_DS2_1 = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 1.d+0)
      call set_trpdens(n, DENS_DUST)
      call set_trradius(n, DG_DS2 * .5d-6)
      call set_fq_aer(n, SOLU_DS2)
      call set_tr_wd_type(n, npart)
      end subroutine N_DS2_1_setSpec

      subroutine M_SSA_SU_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_M_SSA_SU = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 96.d0)
      call set_trpdens(n, DENS_SULF)
      call set_trradius(n, DG_SSA * .5d-6)
      call set_fq_aer(n, SOLU_SSA)
      call set_tr_wd_type(n, npart)
      end subroutine M_SSA_SU_setSpec

      subroutine M_SSA_SS_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_M_SSA_SS = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 75.d0)
      call set_trpdens(n, DENS_SEAS)
      call set_trradius(n, DG_SSA * .5d-6)
      call set_fq_aer(n, SOLU_SSA)
      call set_tr_wd_type(n, npart)
      end subroutine M_SSA_SS_setSpec

      subroutine M_SSC_SS_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_M_SSC_SS = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 75.d0)
      call set_trpdens(n, DENS_SEAS)
      call set_trradius(n, DG_SSC * .5d-6)
      call set_fq_aer(n, SOLU_SSC)
      call set_tr_wd_type(n, npart)
      end subroutine M_SSC_SS_setSpec

      subroutine M_SSS_SS_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_M_SSS_SS = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 75.d0)
      call set_trpdens(n, DENS_SEAS)
      call set_trradius(n, DG_SSS * .5d-6)
      call set_fq_aer(n, SOLU_SSS)
      call set_tr_wd_type(n, npart)
      end subroutine M_SSS_SS_setSpec

      subroutine M_SSS_SU_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_M_SSS_SU = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 75.d0)
      call set_trpdens(n, DENS_SULF)
      call set_trradius(n, DG_SSS * .5d-6)
      call set_fq_aer(n, SOLU_SSS)
      call set_tr_wd_type(n, npart)
      end subroutine M_SSS_SU_setSpec

      subroutine M_OCC_SU_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_M_OCC_SU = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 96.d0)
      call set_trpdens(n, DENS_SULF)
      call set_trradius(n, DG_OCC * .5d-6)
      call set_fq_aer(n, SOLU_OCC)
      call set_tr_wd_type(n, npart)
      end subroutine M_OCC_SU_setSpec

      subroutine M_OCC_OC_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_M_OCC_OC = n
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 15.6d0)
      call set_trpdens(n, DENS_OCAR)
      call set_trradius(n, DG_OCC * .5d-6)
      call set_fq_aer(n,  SOLU_OCC)
      call set_tr_wd_type(n, npart)
      end subroutine M_OCC_OC_setSpec

      subroutine N_OCC_1_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_N_OCC_1 = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 1.d+0)
      call set_trpdens(n, DENS_OCAR)
      call set_trradius(n, DG_OCC * .5d-6)
      call set_fq_aer(n,  SOLU_OCC)
      call set_tr_wd_type(n, npart)
      end subroutine N_OCC_1_setSpec

      subroutine M_BC1_SU_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_M_BC1_SU = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 96.d0)
      call set_trpdens(n, DENS_SULF)
      call set_trradius(n, DG_BC1 * .5d-6)
      call set_fq_aer(n, SOLU_BC1)
      call set_tr_wd_type(n, npart)
      end subroutine M_BC1_SU_setSpec

      subroutine M_BC1_BC_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_M_BC1_BC = n
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 12.d0)
      call set_trpdens(n, DENS_BCAR)
      call set_trradius(n, DG_BC1 * .5d-6)
      call set_fq_aer(n,  SOLU_BC1)
      call set_tr_wd_type(n, npart)
      end subroutine M_BC1_BC_setSpec

      subroutine N_BC1_1_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_N_BC1_1 = n
          call set_ntm_power(n, -11)
          call set_ntsurfsrc(n,  0)
          call set_tr_mm(n, 1.d+0)
          call set_trpdens(n, DENS_BCAR)
          call set_trradius(n, DG_BC1 * .5d-6)
          call set_fq_aer(n,  SOLU_BC1)
          call set_tr_wd_type(n, npart)
      end subroutine N_BC1_1_setSpec

      subroutine M_BC2_SU_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_M_BC2_SU = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 96.d0)
      call set_trpdens(n, DENS_SULF)
      call set_trradius(n, DG_BC2 * .5d-6)
      call set_fq_aer(n, SOLU_BC2)
      call set_tr_wd_type(n, npart)
      end subroutine M_BC2_SU_setSpec

      subroutine M_BC2_BC_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_M_BC2_BC = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 12.d0)
      call set_trpdens(n, DENS_BCAR)
      call set_trradius(n, DG_BC2 * .5d-6)
      call set_fq_aer(n, SOLU_BC2)
      call set_tr_wd_type(n, npart)
      end subroutine M_BC2_BC_setSpec

      subroutine N_BC2_1_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_N_BC2_1 = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 1.d+0)
      call set_trpdens(n, DENS_BCAR)
      call set_trradius(n, DG_BC2 * .5d-6)
      call set_fq_aer(n, SOLU_BC2)
      call set_tr_wd_type(n, npart)
      end subroutine N_BC2_1_setSpec

      subroutine M_BC3_SU_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_M_BC3_SU = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 96.d0)
      call set_trpdens(n, DENS_SULF)
      call set_trradius(n, DG_BC3 * .5d-6)
      call set_fq_aer(n, SOLU_BC3)
      call set_tr_wd_type(n, npart)
      end subroutine M_BC3_SU_setSpec

      subroutine M_BC3_BC_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_M_BC3_BC = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 12.d0)
      call set_trpdens(n, DENS_BCAR)
      call set_trradius(n, DG_BC3 * .5d-6)
      call set_fq_aer(n, SOLU_BC3)
      call set_tr_wd_type(n, npart)
      end subroutine M_BC3_BC_setSpec

      subroutine N_BC3_1_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_N_BC3_1 = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 1.d+0)
      call set_trpdens(n, DENS_BCAR)
      call set_trradius(n, DG_BC3 * .5d-6)
      call set_fq_aer(n, SOLU_BC3)
      call set_tr_wd_type(n, npart)
      end subroutine N_BC3_1_setSpec

      subroutine M_DBC_SU_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_M_DBC_SU = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 96.d0)
      call set_trpdens(n, DENS_SULF)
      call set_trradius(n, DG_DBC * .5d-6)
      call set_fq_aer(n, SOLU_DBC)
      call set_tr_wd_type(n, npart)
      end subroutine M_DBC_SU_setSpec

      subroutine M_DBC_BC_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_M_DBC_BC = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 96.d0)
      call set_trpdens(n, DENS_BCAR)
      call set_trradius(n, DG_DBC * .5d-6)
      call set_fq_aer(n,  SOLU_DBC)
      call set_tr_wd_type(n, npart)
      end subroutine M_DBC_BC_setSpec

      subroutine M_DBC_DU_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_M_DBC_DU = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 1.d+0)
      call set_trpdens(n, DENS_DUST)
      call set_trradius(n, DG_DBC * .5d-6)
      call set_fq_aer(n,  SOLU_DBC)
      call set_tr_wd_type(n, npart)
      end subroutine M_DBC_DU_setSpec

      subroutine N_DBC_1_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_N_DBC_1 = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 1.d+0)
      call set_trpdens(n, DENS_DUST)
      call set_trradius(n, DG_DBC * .5d-6)
      call set_fq_aer(n,  SOLU_DBC)
      call set_tr_wd_type(n, npart)
      end subroutine N_DBC_1_setSpec

      subroutine M_BOC_SU_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_M_BOC_SU = n
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 96.d0)
      call set_trpdens(n, DENS_SULF)
      call set_trradius(n, DG_BOC * .5d-6)
      call set_fq_aer(n, SOLU_BOC)
      call set_tr_wd_type(n, npart)
      end subroutine M_BOC_SU_setSpec

      subroutine M_BOC_BC_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_M_BOC_BC = n
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 12.d0)
      call set_trpdens(n, DENS_BCAR)
      call set_trradius(n, DG_BOC * .5d-6)
      call set_fq_aer(n,  SOLU_BOC)
      call set_tr_wd_type(n, npart)
      end subroutine M_BOC_BC_setSpec

      subroutine M_BOC_OC_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_M_BOC_OC = n
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 15.6d0)
      call set_trpdens(n, DENS_OCAR)
      call set_trradius(n, DG_BOC * .5d-6)
      call set_fq_aer(n,  SOLU_BOC)
      call set_tr_wd_type(n, npart)
      end subroutine M_BOC_OC_setSpec

      subroutine N_BOC_1_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_N_BOC_1 = n
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 1.d+0)
      call set_trpdens(n, DENS_BCAR)
      call set_trradius(n, DG_BOC * .5d-6)
      call set_fq_aer(n,  SOLU_BOC)
      call set_tr_wd_type(n, npart)
      end subroutine N_BOC_1_setSpec

      subroutine M_BCS_SU_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_M_BCS_SU = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 96.d0)
      call set_trpdens(n, DENS_SULF)
      call set_trradius(n, DG_BCS * .5d-6)
      call set_fq_aer(n, SOLU_BCS)
      call set_tr_wd_type(n, npart)
      end subroutine M_BCS_SU_setSpec

      subroutine M_BCS_BC_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_M_BCS_BC = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 12.d0)
      call set_trpdens(n, DENS_BCAR)
      call set_trradius(n, DG_BCS * .5d-6)
      call set_fq_aer(n, SOLU_BCS)
      call set_tr_wd_type(n, npart)
      end subroutine M_BCS_BC_setSpec

      subroutine N_BCS_1_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_N_BCS_1 = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 1.d+0)
      call set_trpdens(n, DENS_BCAR)
      call set_trradius(n, DG_BCS * .5d-6)
      call set_fq_aer(n, SOLU_BCS)
      call set_tr_wd_type(n, npart)
      end subroutine N_BCS_1_setSpec

      subroutine M_MXX_SU_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_M_MXX_SU = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 96.d0)
      call set_trpdens(n, DENS_SULF)
      call set_trradius(n, DG_MXX * .5d-6)
      call set_fq_aer(n, SOLU_MXX)
      call set_tr_wd_type(n, npart)
      end subroutine M_MXX_SU_setSpec

      subroutine M_MXX_BC_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_M_MXX_BC = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 12.d0)
      call set_trpdens(n, DENS_BCAR)
      call set_trradius(n, DG_MXX * .5d-6)
      call set_fq_aer(n, SOLU_MXX)
      call set_tr_wd_type(n, npart)
      end subroutine M_MXX_BC_setSpec

      subroutine M_MXX_OC_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_M_MXX_OC = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 15.6d0)
      call set_trpdens(n, DENS_OCAR)
      call set_trradius(n, DG_MXX * .5d-6)
      call set_fq_aer(n, SOLU_MXX)
      call set_tr_wd_type(n, npart)
      end subroutine M_MXX_OC_setSpec

      subroutine M_MXX_DU_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_M_MXX_DU = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 1.d+0)
      call set_trpdens(n, DENS_DUST)
      call set_trradius(n, DG_MXX * .5d-6)
      call set_fq_aer(n, SOLU_MXX)
      call set_tr_wd_type(n, npart)
      end subroutine M_MXX_DU_setSpec

      subroutine M_MXX_SS_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_M_MXX_SS = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 75.d0)
      call set_trpdens(n, DENS_SEAS)
      call set_trradius(n, DG_MXX * .5d-6)
      call set_fq_aer(n, SOLU_MXX)
      call set_tr_wd_type(n, npart)
      end subroutine M_MXX_SS_setSpec

      subroutine N_MXX_1_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_N_MXX_1 = n
#ifdef TRACERS_AMP
      ntmAMPe=n                 ! always the last tracer in AMP
      if (ntmAMPi+ntmAMP-1 /= ntmAMPe)
     *     call stop_model( 'ntmAMPi+ntmAMP-1 /= ntmAMPe', 255 )
#endif
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 1.d+0)
      call set_trpdens(n, DENS_SEAS)
      call set_trradius(n, DG_MXX * .5d-6)
      call set_fq_aer(n, SOLU_MXX)
      call set_tr_wd_type(n, npart)
      end subroutine N_MXX_1_setSpec

      subroutine M_OCS_SU_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_M_OCS_SU = n
          call set_ntm_power(n, -11)
          call set_ntsurfsrc(n,  0)
          call set_tr_mm(n, 96.d0)
          call set_trpdens(n, DENS_SULF)
          call set_trradius(n, DG_OCS * .5d-6)
          call set_fq_aer(n, SOLU_OCS)
          call set_tr_wd_type(n, npart)
      end subroutine M_OCS_SU_setSpec

      subroutine M_OCS_OC_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_M_OCS_OC = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 15.6d0)
      call set_trpdens(n, DENS_OCAR)
      call set_trradius(n, DG_OCS * .5d-6)
      call set_fq_aer(n, SOLU_OCS)
      call set_tr_wd_type(n, npart)
      end subroutine M_OCS_OC_setSpec

      subroutine N_OCS_1_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_N_OCS_1 = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 1.d+0)
      call set_trpdens(n, DENS_OCAR)
      call set_trradius(n, DG_OCS * .5d-6)
      call set_fq_aer(n, SOLU_OCS)
      call set_tr_wd_type(n, npart)
      end subroutine N_OCS_1_setSpec

#endif /* TRACERS_AMP */


#ifdef TRACERS_TOMAS
      subroutine TOMAS_H2SO4_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_H2SO4 = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 98.d0)
      call set_trpdens(n, 1.78d0)
      call set_fq_aer(n, 1.d0)
      call set_tr_wd_type(n, nGas)
      end subroutine TOMAS_H2SO4_setSpec

      subroutine TOMAS_SOAgas_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_SOAgas = n
      call set_ntm_power(n, -11)
!     call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 120.10d0) ! i.e. 10 carbons
#ifdef TRACERS_DRYDEP
      call set_HSTAR(n,  0.D0)  !no dry dep
#endif
      end subroutine TOMAS_SOAgas_setSpec

      function TOMAS_setSpec(func, name, nbins) result (indices)
      interface
        integer function func(name, bin)
        character(len=*), intent(in) :: name
        integer, intent(in) :: bin
        end function func
      end interface
      character(len=*), intent(in) :: name
      integer, intent(in) :: nbins
      integer :: indices(nbins)

      integer :: bin
      character(len=len_trim(name) + 4) :: fullName

      do bin = 1, nbins
        if (len_trim(name) == 5) then
          write(fullName,'(a,"_",i2.2)') trim(name), bin
        else
          write(fullName,'(a,"__",i2.2)') trim(name), bin
        end if
        indices(bin) = func(fullName, bin)
      end do
      
      end function TOMAS_setSpec

      integer function TOMAS_ANUM_setSpec(name, bin) result(n_ANUM)
      character(len=*), intent(in) :: name
      integer, intent(in) :: bin
      n = oldAddTracer(name)
      n_ANUM = n  
      
      TOMAS_dens=1.5d3
      TOMAS_radius=(sqrt(xk(bin)*xk(bin+1))/TOMAS_dens/pi/4.*3.)
     *     **(1./3.)  
      
!#ifdef TOMAS_12_10NM 
      if(bin.le.5) call set_ntm_power(n, 10)
      if(bin.gt.5) call set_ntm_power(n, 8) 
!#endif
!#ifdef TOMAS_12_3NM 
!      if(bin.le.8) call set_ntm_power(n, 10)
!      if(bin.gt.8) call set_ntm_power(n, 8) 
!#endif

      call set_ntsurfsrc(n,  3) ! 1 for SO4,2 for EC, 3 for OC (4 for SS and 5 for DU)
      call set_tr_mm(n, 1.d+0)
      call set_trpdens(n, TOMAS_dens)
      call set_trradius(n, TOMAS_radius)
      call set_fq_aer(n, 1.d0)  !not used in wet deposition
      call set_tr_wd_type(n, npart)
      end function TOMAS_ANUM_setSpec

      integer function TOMAS_ASO4_setSpec(name, bin) result(n_ASO4)
      character(len=*), intent(in) :: name
      integer, intent(in) :: bin
      n = oldAddTracer(name)
      n_ASO4 = n 
      TOMAS_dens=1.78d3
      TOMAS_radius=(sqrt(xk(bin)*xk(bin+1))/TOMAS_dens/pi/4.*3.) 
     &     **(1./3.) 

      call set_ntm_power(n, -11)
      call set_tr_mm(n, 96.d+0)
      call set_trpdens(n, TOMAS_dens) !kg/m3 this is sulfate value
      call set_trradius(n, TOMAS_radius) !m
      call set_fq_aer(n, 1.d0   ) !not used in wet deposition
      call set_tr_wd_type(n, npart)
      end function TOMAS_ASO4_setSpec

      integer function TOMAS_ANACL_setSpec(name, bin) result(n_ANACL)
      character(len=*), intent(in) :: name
      integer, intent(in) :: bin
      n = oldAddTracer(name)
      n_ANACL = n         
      call set_ntsurfsrc(n,  0) ! ocean bubbles
      call set_ntisurfsrc(n, 1)
      TOMAS_dens= 2.165d3
      TOMAS_radius=(sqrt(xk(bin)*xk(bin+1))/TOMAS_dens/pi/4.*3.)
     &     **(1./3.) 
      if(bin.le.10) call set_ntm_power(n, -10)
      if(bin.gt.10) call set_ntm_power(n, -8)
      call set_tr_mm(n, 75.d+0)
      call set_trpdens(n, TOMAS_dens) !kg/m3 this is sulfate value
      call set_trradius(n, TOMAS_radius) !m
      call set_fq_aer(n, 1.d0   ) !not used in wet deposition
      call set_tr_wd_type(n, npart)
      end function TOMAS_ANACL_setSpec

      integer function TOMAS_AECOB_setSpec(name, bin) result(n_AECOB)
      character(len=*), intent(in) :: name
      integer, intent(in) :: bin
      n = oldAddTracer(name)
      n_AECOB = n          

      TOMAS_dens= 1.8d3
      TOMAS_radius=(sqrt(xk(bin)*xk(bin+1))/TOMAS_dens/pi/4.*3.)
     &     **(1./3.) 
      call set_ntm_power(n, -12)
      call set_tr_mm(n, 12.d+0)
      call set_trpdens(n, TOMAS_dens) !kg/m3 this is sulfate value
      call set_trradius(n, TOMAS_radius) !m
      call set_fq_aer(n, 1.d0   ) !not used in wet deposition
      call set_tr_wd_type(n, npart)
      end function TOMAS_AECOB_setSpec

      integer function TOMAS_AECIL_setSpec(name, bin) result(n_AECIL)
      character(len=*), intent(in) :: name
      integer, intent(in) :: bin
      n = oldAddTracer(name)
      n_AECIL = n        
      TOMAS_dens= 1.8d3
      TOMAS_radius=(sqrt(xk(bin)*xk(bin+1))/TOMAS_dens/pi/4.*3.)
     &     **(1./3.) 
      call set_ntm_power(n, -12)
      call set_tr_mm(n, 12.d+0)
      call set_trpdens(n, TOMAS_dens) !kg/m3 this is sulfate value
      call set_trradius(n, TOMAS_radius) !m
      call set_fq_aer(n, 1.d0   ) !not used in wet deposition
      call set_tr_wd_type(n, npart)
      end function TOMAS_AECIL_setSpec

      integer function TOMAS_AOCOB_setSpec(name, bin) result(n_AOCOB)
      use OldTracer_mod, only: om2oc, set_om2oc
      character(len=*), intent(in) :: name
      integer, intent(in) :: bin
      real*8 :: tmp

      n = oldAddTracer(name)
      n_AOCOB = n 
      TOMAS_dens= 1.4d3
      TOMAS_radius=(sqrt(xk(bin)*xk(bin+1))/TOMAS_dens/pi/4.*3.) 
     &     **(1./3.) 
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 200.d+0)
      call set_trpdens(n, TOMAS_dens) !kg/m3 this is sulfate value
      call set_trradius(n, TOMAS_radius) !m
      call set_fq_aer(n, 1.d0   ) !not used in wet deposition
      call set_tr_wd_type(n, npart)        

      if (bin==1) then
        tmp = om2oc(n_AOCOB)
        call sync_param("OCB_om2oc",tmp)
        call set_om2oc(n_AOCOB, tmp)
      end if
      end function TOMAS_AOCOB_setSpec

      integer function TOMAS_AOCIL_setSpec(name, bin) result(n_AOCIL)
      character(len=*), intent(in) :: name
      integer, intent(in) :: bin
      n = oldAddTracer(name)
      n_AOCIL = n  
      TOMAS_dens= 1.4d3
      TOMAS_radius=(sqrt(xk(bin)*xk(bin+1))/TOMAS_dens/pi/4.*3.)
     &     **(1./3.) 
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 200.d+0)
      call set_trpdens(n, TOMAS_dens) !kg/m3 this is sulfate value
      call set_trradius(n, TOMAS_radius) !m
      call set_fq_aer(n, 1.d0   ) !not used in wet deposition
      call set_tr_wd_type(n, npart)        
      end function TOMAS_AOCIL_setSpec

      integer function TOMAS_ADUST_setSpec(name, bin) result(n_ADUST)
      character(len=*), intent(in) :: name
      integer, intent(in) :: bin
      n = oldAddTracer(name)
      n_ADUST = n  
      call set_ntsurfsrc(n,  0)
      call set_ntisurfsrc(n, 1)
      if(bin.le.10) TOMAS_dens= 2.5d3 !clay 
      if(bin.gt.10) TOMAS_dens= 2.65d3 !Silt
      TOMAS_radius=(sqrt(xk(bin)*xk(bin+1))/TOMAS_dens/pi/4.*3.) 
     &     **(1./3.) 

      if(bin.le.9) call set_ntm_power(n, -11)
      if(bin.gt.9) call set_ntm_power(n, -9)
      call set_tr_mm(n, 1.d+0)
      call set_trpdens(n, TOMAS_dens) !kg/m3 this is sulfate value
      call set_trradius(n, TOMAS_radius) !m
      call set_fq_aer(n, 1.d0   ) !not used in wet deposition
      call set_tr_wd_type(n, npart)  
      end function TOMAS_ADUST_setSpec

      integer function TOMAS_AH2O_setSpec(name, bin) result(n_AH2O)
      character(len=*), intent(in) :: name
      integer, intent(in) :: bin
      n = oldAddTracer(name)
      n_AH2O = n         
      TOMAS_dens= 1.d3
      TOMAS_radius=(sqrt(xk(bin)*xk(bin+1))/TOMAS_dens/pi/4.*3.)
     *     **(1./3.) 
      call set_ntsurfsrc(n,  0)
      call set_ntm_power(n, -8)
      call set_tr_mm(n, 18.d+0)
      call set_trpdens(n, TOMAS_dens) !kg/m3 this is sulfate value
      call set_trradius(n, TOMAS_radius) !m
      call set_fq_aer(n, 1.d0   ) !not used in wet deposition
      call set_tr_wd_type(n, npart)  
      end function TOMAS_AH2O_setSpec

#endif /* TRACERS_TOMAS */

      end subroutine initTracerMetadata


      subroutine laterInitTracerMetadata()
      USE MODEL_COM, only: itime,master_yr
      use OldTracer_mod, only: itime_tr0
      use OldTracer_mod, only: set_itime_tr0
      USE TRACER_COM, only: NTM, tracers, syncProperty
#ifdef TRACERS_TOMAS
      use TOMAS_AEROSOL, only : binact10,binact02,
     &     fraction10,fraction02
#endif
      use TRACER_COM, only: coupled_chem
      use Dictionary_mod, only: sync_param,is_set_param,get_param
#ifdef TRACERS_WATER
      use TRDIAG_com, only: to_per_mil
#endif
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP) ||\
      (defined TRACERS_TOMAS)
      USE AEROSOL_SOURCES, only: tune_ss1, tune_ss2, om2oc, BBinc
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
     &     gw_llim,sw_lim,exclude_us_eu,nn_or_zon,ice_age,nday_ch4,max_days,
     &     ns_wet,nra_ch4
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
