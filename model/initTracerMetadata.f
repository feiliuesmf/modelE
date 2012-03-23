#include "rundeck_opts.h"
      SUBROUTINE initTracerMetadata()
!     @sum init_tracer initializes trace gas attributes and diagnostics
!     @auth J. Lerner
!     @calls sync_param, SET_TCON, RDLAND, RDDRYCF
      USE DOMAIN_DECOMP_ATM, only: AM_I_ROOT
      USE CONSTANT, only: mair,mwat,sday,pi
#ifdef TRACERS_AEROSOLS_SOA
     &     ,gasc
#endif  /* TRACERS_AEROSOLS_SOA */
      USE RESOLUTION, only : jm,lm
      USE MODEL_COM, only: dtsrc,itime
      USE ATM_COM, only: pmidl00
      USE GEOM, only: axyp,byaxyp
      USE ATM_COM, only: am     ! Air mass of each box (kg/m^2)
      USE TRACER_COM
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
#ifdef TRACERS_WATER
#ifdef TRDIAG_WETDEPO
      USE CLOUDS, ONLY : diag_wetdep
#endif
#endif /* TRACERS_WATER */
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
      use OldTracer_mod, only: addTracer
      use OldTracer_mod, only: set_tr_mm, set_ntm_power
      use OldTracer_mod, only: set_t_qlimit
      use OldTracer_mod, only: set_ntsurfsrc
      use OldTracer_mod, only: set_needtrs
      use OldTracer_mod, only: set_trdecay
      use OldTracer_mod, only: set_itime_tr0
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

#if (defined TRACERS_OCEAN) &&  !defined(TRACERS_OCEAN_INDEP)
!     atmosphere copies atmosphere-declared tracer info to ocean
!     so that the ocean can "inherit" it without referencing atm. code
      use ocn_tracer_com, only : 
     &     n_Water_ocn      => n_Water,
     &     itime_tr0_ocn    => itime_tr0,
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
      integer :: l,k,n,kr,m,ns,numTr
#ifdef TRACERS_SPECIAL_O18
      real*8 fracls
#endif
#ifdef TRACERS_TOMAS
      real*8 :: TOMAS_dens,TOMAS_radius
#endif
#if (defined TRACERS_WATER) || (defined TRACERS_DRYDEP)
!     @param convert_HSTAR converts from mole/Joule to mole/(L*atm)
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
      integer :: values(ntm)
      integer :: val

      INTEGER J_0, J_1, I_0, I_1

C**** 
C**** Set some documentary parameters in the database
C**** 
      call set_param("NTM",NTM,'o')
      call set_param("TRNAME",TRNAME,ntm,'o')

      call sync_param( "COUPLED_CHEM", COUPLED_CHEM )

#ifdef TRACERS_TOMAS
      CALL initbounds()
      call readbinact ("binact10_12.dat",binact10) 
      call readbinact ("binact02_12.dat",binact02) 
      call readfraction("fraction10_12.dat",fraction10) 
      call readfraction("fraction02_12.dat",fraction02) 
      call readmielut           ! aerosol radiation lookup table
#endif

      call initializeOldTracers()
      do n=1,ntm
        call addTracer(trname(n))
      end do


!     TODO remove this comment
C**** Set defaults for tracer attributes (all dimensioned ntm)
C**** Many defaults are now set in OldTracers_mod.F90

#ifdef TRACERS_ON
#ifdef TRACERS_SPECIAL_Lerner
      n_MPtable = 0
      tcscale = 0.
#endif
#endif /* TRACERS_ON */
#ifdef TRACERS_SPECIAL_O18
      iso_index = 1             ! act like water by default
#endif

C**** Synchronise tracer related paramters from rundeck

C**** Get itime_tr0 from rundeck if it exists
      values = iTime
      call sync_param("itime_tr0",values,ntm)
      do n = 1, ntm
        call set_itime_tr0(n, values(n))
      end do

#ifdef TRACERS_WATER
C**** Decide on water tracer conc. units from rundeck if it exists
      call sync_param("to_per_mil",to_per_mil,ntm)
#endif
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP) ||\
      (defined TRACERS_TOMAS)
      call sync_param("tune_ss1",tune_ss1)
      call sync_param("tune_ss2",tune_ss2)
      call sync_param("BBinc",BBinc)
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
      call sync_param("PI_run",PI_run)
      call sync_param("PIratio_N",PIratio_N)
      call sync_param("PIratio_CO_T",PIratio_CO_T)
      call sync_param("PIratio_CO_S",PIratio_CO_S)
      call sync_param("PIratio_other",PIratio_other)
      call sync_param("rad_FL",rad_fl)
      call sync_param("use_rad_ch4",use_rad_ch4)
      call sync_param("Lmax_rad_O3",Lmax_rad_O3)
      call sync_param("Lmax_rad_CH4",Lmax_rad_CH4)
      call sync_param("aircraft_Tyr1",aircraft_Tyr1)
      call sync_param("aircraft_Tyr2",aircraft_Tyr2)
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

#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP) ||\
      (defined TRACERS_TOMAS)
C**** DMS, seasalt from offline fields
      call sync_param("OFFLINE_DMS_SS",OFFLINE_DMS_SS)
C**** seasalt from offline fields
      call sync_param("OFFLINE_SS",OFFLINE_SS)
C**** decide if preindustrial emissions
      call sync_param("imPI",imPI)
C**** determine year of emissions
      call sync_param("aer_int_yr",aer_int_yr)
#endif
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

!     call routine to read/set up regions and sectors for emissions:
      call setup_emis_sectors_regions

!     initialize 3D source factors:
      ef_fact3d(:,:)=1.d0

!     The following section will check for rundeck file of
!     the form: trname_01, trname_02... and thereby define
!     the ntsurfsrc(n). If those files exist it reads an
!     80 char header to get information including the
!     source name (ssame-->{sname,lname,etc.}. ntsurfsrc(n)
!     get set to zero if those files aren't found:
!     (I can enclose this in an ifdef if it causes problems
!     for people). num_srf_sources routine also assigns
!     sources to sectors, if desired:
      do n=1,ntm
!     general case:
#if (defined TRACERS_AMP) || (defined TRACERS_TOMAS) 
        call num_srf_sources(n,.false.)
#else
        call num_srf_sources(n,.true.)
#endif

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
          call set_ntsurfsrc(n, ntsurfsrc(n)+1)
          ssname(n,ntsurfsrc(n))="Terpene_source"
        end select
#endif
#endif  /* TRACERS_AEROSOLS_SOA */
#ifdef TRACERS_SPECIAL_Shindell
        if (trname(n)=='CH4' .and. use_rad_ch4/=0) then
          call set_ntsurfsrc(n,0)
        end if
#endif
      enddo

C**** Define individual tracer characteristics
      numTr = 0
      if (tracers_special_shindell) then
#ifdef TRACERS_SPECIAL_Shindell
        n_Ox = Ox_setSpec(numTr)
        n_NOx = NOx_setSpec(numTr)
        n_ClOx = ClOx_setSpec(numTr)
        n_BrOx = BrOx_setSpec(numTr)
        n_N2O5 = N2O5_setSpec(numTr)
        n_HNO3 = HNO3_setSpec(numTr)
        n_H2O2 = H2O2_setSpec(numTr)
        n_CH3OOH = CH3OOH_setSpec(numTr)
        n_HCHO = HCHO_setSpec(numTr)
        n_HO2NO2 = HO2NO2_setSpec(numTr)
        n_CO = CO_setSpec(numTr)
        n_CH4 = CH4_setSpec(numTr)
        n_PAN = PAN_setSpec(numTr)
        n_Isoprene = Isoprene_setSpec(numTr)
        n_AlkylNit = AlkylNit_setSpec(numTr)
        n_Alkenes = Alkenes_setSpec(numTr)
        n_Paraffin = Paraffin_setSpec(numTr)

        if (tracers_terp) then
          n_Terpenes = Terpenes_setSpec(numTr)
        end if

#ifdef TRACERS_AEROSOLS_SOA
        if (tracers_aerosols_soa) then
          n_isopp1g = isopp1g_setSpec(numTr)
          n_isopp1a = isopp1a_setSpec(numTr)
          n_isopp2g = isopp2g_setSpec(numTr)
          n_isopp2a = isopp2a_setSpec(numTr)
          if (tracers_terp) then
            n_apinp1g = apinp1g_setSpec(numTr)
            n_apinp1a = apinp1a_setSpec(numTr)
            n_apinp2g = apinp2g_setSpec(numTr)
            n_apinp2a = apinp2a_setSpec(numTr)
          end if
        end if
#endif

        n_HCl = HCl_setSpec(numTr)
        n_HOCl = HOCl_setSpec(numTr)
        n_ClONO2 = ClONO2_setSpec(numTr)
        n_HBr = HBr_setSpec(numTr)
        n_HOBr = HOBr_setSpec(numTr)
        n_BrONO2 = BrONO2_setSpec(numTr)
        n_N2O = N2O_setSpec(numTr)
        n_CFC = CFC_setSpec(numTr)

        if (shindell_strat_extra) then
          if (accmip_like_diags) then
            n_codirect = codirect_setSpec(numTr)
            n_stratOx = stratOx_setSpec(numTr)
            n_GLT = GLT_setSpec(numTr) ! generic linear tracer
          end if
        end if

#endif /* TRACERS_SPECIAL_Shindell */
      end if  ! tracers_special_shindell

      if ((.not. tracers_amp) .and. tracers_water) then
        n_Water = Water_setSpec(numTr)
      end if

      

#ifdef TRACERS_SPECIAL_O18
      n_H2O18 = H2O18_setSpec(numTr)
      n_HDO = HDO_setSpec(numTr)
c$$$      n_HTO = HTO_setSpec(numTr)
c$$$      n_H2O17 = H2O17_setSpec(numTr)
#endif  /* TRACERS_SPECIAL_O18 */

      if (tracers_gasexch_ocean_cfc) then
        n_CFCn = CFCn_setSpec(numTr)
      end if

      if (tracers_gasexch_ocean_co2 .or. tracers_gasexch_land_co2) then
        n_CO2n = CO2n_setSpec(numTr)
      end if
      
      if (tracers_special_lerner) then
        n_SF6 = SF6_setSpec(numTr)
        n_Rn222 = Rn222_setSpec(numTr)
        n_CO2 = CO2_setSpec(numTr)
        n_N2O = N2O_setSpec(numTr) ! warning duplicate with Shindell
        n_CFC11 = CFC11_setSpec(numTr)
        n_14CO2 = C_14O2_setSpec(numTr)
        n_CH4 = CH4_setSpec(numTr) ! warning duplicate with Shindell
        n_O3 = O3_setSpec(numTr)
        n_SF6_c = SF6_c_setSpec(numTr)
        if (tracers_special_shindell) 
     &       call stop_model('contradictory tracer specs')
      end if

      if (tracers_aerosols_koch) then
        n_DMS = DMS_setSpec(numTr)
        n_MSA = MSA_setSpec(numTr)
        n_SO2 = SO2_setSpec(numTr)
        n_SO4 = SO4_setSpec(numTr)
        n_H2O2_s = H2O2_s_setSpec(numTr)
        if (.not. sulf_only_aerosols) then
          n_seasalt1 = seasalt1_setSpec(numTr)
          n_seasalt2 = seasalt2_setSpec(numTr)
          n_BCII = BCII_setSpec(numTr)
          n_BCIA = BCIA_setSpec(numTr)
          n_BCB = BCB_setSpec(numTr)
        end if
#ifdef TRACERS_AEROSOLS_VBS
          n_vbsGm2 = VBS_setSpec('vbsGm2', numTr, 1,'igas')
          n_vbsGm1 = VBS_setSpec('vbsGm1', numTr, 2,'igas')
          n_vbsGz  = VBS_setSpec('vbsGz',  numTr, 3,'igas')
          n_vbsGp1 = VBS_setSpec('vbsGp1', numTr, 4,'igas')
          n_vbsGp2 = VBS_setSpec('vbsGp2', numTr, 5,'igas')
          n_vbsGp3 = VBS_setSpec('vbsGp3', numTr, 6,'igas')
          n_vbsGp4 = VBS_setSpec('vbsGp4', numTr, 7,'igas')
          n_vbsGp5 = VBS_setSpec('vbsGp5', numTr, 8,'igas')
          n_vbsGp6 = VBS_setSpec('vbsGp6', numTr, 9,'igas')

          n_vbsAm2 = VBS_setSpec('vbsAm2', numTr, 1,'iaer')
          n_vbsAm1 = VBS_setSpec('vbsAm1', numTr, 2,'iaer')
          n_vbsAz  = VBS_setSpec('vbsAz',  numTr, 3,'iaer')
          n_vbsAp1 = VBS_setSpec('vbsAp1', numTr, 4,'iaer')
          n_vbsAp2 = VBS_setSpec('vbsAp2', numTr, 5,'iaer')
          n_vbsAp3 = VBS_setSpec('vbsAp3', numTr, 6,'iaer')
          n_vbsAp4 = VBS_setSpec('vbsAp4', numTr, 7,'iaer')
          n_vbsAp5 = VBS_setSpec('vbsAp5', numTr, 8,'iaer')
          n_vbsAp6 = VBS_setSpec('vbsAp6', numTr, 9,'iaer')
#else
          n_OCII = OCII_setSpec(numTr)   !Insoluble industrial organic mass
          n_OCIA = OCIA_setSpec(numTr)   !Aged industrial organic mass
          n_OCB = OCB_setSpec(numTr)     !Biomass organic mass
#endif /* TRACERS_AEROSOLS_VBS */
      end if        

      if (tracers_aerosols_ocean) then
        n_OCocean = OCocean_setSpec(numTr) !Insoluble oceanic organic mass
      end if

      if (cpp_tracers_dust) then
        n_clay = clay_setSpec(numTr)
        n_silt1 = silt1_setSpec(numTr)
        n_silt2 = silt2_setSpec(numTr)
        n_silt3 = silt3_setSpec(numTr)
        if (tracers_dust_silt4) n_silt4 = silt4_setSpec(numTr)
      end if

      if (tracers_nitrate) then
        n_NH3 = NH3_setSpec(numTr)
        n_NH4 = NH4_setSpec(numTr)
        n_NO3p = NO3p_setSpec(numTr)
      end if

      if (tracers_hetchem) then
        n_SO4_d1 = SO4_d1_setSpec(numTr)
        n_SO4_d2 = SO4_d2_setSpec(numTr)
        n_SO4_d3 = SO4_d3_setSpec(numTr)
        if (tracers_nitrate) then
          n_N_d1 = N_d1_setSpec(numTr)
          n_N_d2 = N_d2_setSpec(numTr)
          n_N_d3 = N_d3_setSpec(numTr)
        end if
      end if

      if (tracers_cosmo) then
        if (tracers_radon) n_Pb210 = Pb210_setSpec(numTr)
        n_Be7 = Be7_setSpec(numTr)
        n_Be10 = Be10_setSpec(numTr)
        if (tracers_radon) n_Rn222 = Rn222_setSpec(numTr) ! duplicate with Lerner
        if (tracers_special_lerner) 
     &       call stop_model('contradictory tracer specs')
      end if

#ifdef TRACERS_MINERALS
      if (tracers_minerals) then
        n_clayilli = clayilli_setSpec(numTr) ! http://webmineral.com/data/Illite.shtml
        n_claykaol = claykaol_setSpec(numTr) ! http://www.webmineral.com/data/Kaolinite.shtml
        n_claysmec = claysmec_setSpec(numTr) ! http://www.webmineral.com/data/Rectorite.shtml
        n_claycalc = claycalc_setSpec(numTr) ! http://www.webmineral.com/data/Calcite.shtml
        n_clayquar = clayquar_setSpec(numTr) ! http://www.webmineral.com/data/Quartz.shtml
        n_sil1quar = sil1quar_setSpec(numTr) ! http://www.webmineral.com/data/Quartz.shtml
        n_sil1feld = sil1feld_setSpec(numTr) ! http://www.mindat.org/min-1624.html
        n_sil1calc = sil1calc_setSpec(numTr) ! http://www.webmineral.com/data/Calcite.shtml
        n_sil1hema = sil1hema_setSpec(numTr) ! http://www.webmineral.com/data/Hematite.shtml
        n_sil1gyps = sil1gyps_setSpec(numTr) ! http://www.webmineral.com/data/Gypsum.shtml
        n_sil2quar = sil2quar_setSpec(numTr) ! http://www.webmineral.com/data/Quartz.shtml
        n_sil2feld = sil2feld_setSpec(numTr) ! http://www.mindat.org/min-1624.html
        n_sil2calc = sil2calc_setSpec(numTr) ! http://www.webmineral.com/data/Calcite.shtml
        n_sil2hema = sil2hema_setSpec(numTr) ! http://www.webmineral.com/data/Hematite.shtml
        n_sil2gyps = sil2gyps_setSpec(numTr) ! http://www.webmineral.com/data/Gypsum.shtml
        n_sil3quar = sil3quar_setSpec(numTr) ! http://www.webmineral.com/data/Quartz.shtml
        n_sil3feld = sil3feld_setSpec(numTr) ! http://www.mindat.org/min-1624.html
        n_sil3calc = sil3calc_setSpec(numTr) ! http://www.webmineral.com/data/Calcite.shtml
        n_sil3hema = sil3hema_setSpec(numTr) ! http://www.webmineral.com/data/Hematite.shtml
        n_sil3gyps = sil3gyps_setSpec(numTr) ! http://www.webmineral.com/data/Gypsum.shtml
      end if
#endif  /* TRACERS_MINERALS */

#ifdef TRACERS_QUARZHEM
      if (tracers_quarzhem) then
        n_sil1quhe = sil1quhe_setSpec(numTr)
        n_sil2quhe = sil2quhe_setSpec(numTr)
        n_sil3quhe = sil3quhe_setSpec(numTr)
      end if
#endif

      if (tracers_air .or. htap_like_diags) then
        n_air = air_setSpec(numTr)
      end if

#ifdef TRACERS_AMP
      if (tracers_amp) then

C**** Tracers for Scheme AMP: Aerosol Microphysics (Mechanism M1 - M8)
        n_M_NO3 = M_NO3_setSpec(numTr)
        n_M_NH4 = M_NH4_setSpec(numTr)
        n_M_H2O = M_H2O_setSpec(numTr)

        if ( tracers_amp_m1 .or. tracers_amp_m2 .or. tracers_amp_m3
     &       .or. tracers_amp_m5 .or. tracers_amp_m6 .or. 
     &       tracers_amp_m7) then
          n_M_AKK_SU = M_AKK_SU_setSpec(numTr)
          n_N_AKK_1 = N_AKK_1_setSpec(numTr)
        end if

        n_M_ACC_SU = M_ACC_SU_setSpec(numTr)
        n_N_ACC_1 = N_ACC_1_setSpec(numTr)
        n_M_DD1_SU = M_DD1_SU_setSpec(numTr)
        n_M_DD1_DU = M_DD1_DU_setSpec(numTr)
        n_N_DD1_1 = N_DD1_1_setSpec(numTr)
        n_M_DS1_SU = M_DS1_SU_setSpec(numTr)
        n_M_DS1_DU = M_DS1_DU_setSpec(numTr)
        n_N_DS1_1 = N_DS1_1_setSpec(numTr)

        if (tracers_amp_m1 .or. tracers_amp_m2 .or. tracers_amp_m3
     &       .or. tracers_amp_m4) then
          n_M_DD2_SU = M_DD2_SU_setSpec(numTr)
          n_M_DD2_DU = M_DD2_DU_setSpec(numTr)
          n_N_DD2_1 = N_DD2_1_setSpec(numTr)
          n_M_DS2_SU = M_DS2_SU_setSpec(numTr)
          n_M_DS2_DU = M_DS2_DU_setSpec(numTr)
          n_N_DS2_1 = N_DS2_1_setSpec(numTr)
        end if

        if (tracers_amp_m1 .or.  tracers_amp_m2 .or. tracers_amp_m3
     &       .or. tracers_amp_m5 .or.  tracers_amp_m6 .or. 
     &       tracers_amp_m7) then
          n_M_SSA_SU = M_SSA_SU_setSpec(numTr)
          n_M_SSA_SS = M_SSA_SS_setSpec(numTr)
          n_M_SSC_SS = M_SSC_SS_setSpec(numTr)
        end if

        if ( tracers_amp_m4 .or. tracers_amp_m8) then
          n_M_SSS_SU = M_SSS_SU_setSpec(numTr)
          n_M_SSS_SS = M_SSS_SS_setSpec(numTr)
        end if

        n_M_OCC_SU = M_OCC_SU_setSpec(numTr)
        n_M_OCC_OC = M_OCC_OC_setSpec(numTr)
        n_N_OCC_1 = N_OCC_1_setSpec(numTr)
        n_M_BC1_SU = M_BC1_SU_setSpec(numTr)
        n_M_BC1_BC = M_BC1_BC_setSpec(numTr)
        n_N_BC1_1 = N_BC1_1_setSpec(numTr)
        n_M_BC2_SU = M_BC2_SU_setSpec(numTr)
        n_M_BC2_BC = M_BC2_BC_setSpec(numTr)
        n_N_BC2_1 = N_BC2_1_setSpec(numTr)

        if ( tracers_amp_m1 .or. tracers_amp_m5) then
          n_M_BC3_SU = M_BC3_SU_setSpec(numTr)
          n_M_BC3_BC = M_BC3_BC_setSpec(numTr)
          n_N_BC3_1 = N_BC3_1_setSpec(numTr)
        end if

        if ( tracers_amp_m2 .or. tracers_amp_m6) then
          n_M_OCS_SU = M_OCS_SU_setSpec(numTr)
          n_M_OCS_OC = M_OCS_OC_setSpec(numTr)
          n_N_OCS_1 = N_OCS_1_setSpec(numTr)
        end if
          
        if ( tracers_amp_m1 .or. tracers_amp_m2 .or. tracers_amp_m6)then
          n_M_DBC_SU = M_DBC_SU_setSpec(numTr)
          n_M_DBC_BC = M_DBC_BC_setSpec(numTr)
          n_M_DBC_DU = M_DBC_DU_setSpec(numTr)
          n_N_DBC_1 = N_DBC_1_setSpec(numTr)
        end if

        if (tracers_amp_m1 .or. tracers_amp_m2 .or. tracers_amp_m3
     &       .or. tracers_amp_m6 .or. tracers_amp_m7) then
          n_M_BOC_SU = M_BOC_SU_setSpec(numTr)
          n_M_BOC_BC = M_BOC_BC_setSpec(numTr)
          n_M_BOC_OC = M_BOC_OC_setSpec(numTr)
          n_N_BOC_1 = N_BOC_1_setSpec(numTr)
        end if

        if ( tracers_amp_m1 .or. tracers_amp_m2 .or. tracers_amp_m5
     &       .or. TRACERS_AMP_M6) then
          n_M_BCS_SU = M_BCS_SU_setSpec(numTr)
          n_M_BCS_BC = M_BCS_BC_setSpec(numTr)
          n_N_BCS_1 = N_BCS_1_setSpec(numTr)
        end if

        n_M_MXX_SU = M_MXX_SU_setSpec(numTr)
        n_M_MXX_BC = M_MXX_BC_setSpec(numTr)
        n_M_MXX_OC = M_MXX_OC_setSpec(numTr)
        n_M_MXX_DU = M_MXX_DU_setSpec(numTr)
        n_M_MXX_SS = M_MXX_SS_setSpec(numTr)
        n_N_MXX_1 = N_MXX_1_setSpec(numTr)

        n_H2SO4 = H2SO4_setSpec(numTr)
        n_DMS = DMS_setSpec(numTr) ! duplicate with Koch
        n_SO2 = SO2_setSpec(numTr) ! duplicate with Koch
        n_H2O2_s = H2O2_s_setSpec(numTr) ! duplicate with Koch
        n_NH3 = NH3_setSpec(numTr) ! duplicate with nitrate
        if (tracers_aerosols_koch) then
          call stop_model('contradictory tracer specs')
        end if
        if (tracers_nitrate) then
          call stop_model('contradictory tracer specs')
        end if
      end if
#endif /* TRACERS_AMP */

#ifdef TRACERS_TOMAS
!For aerosol tracers in TOMAS model, 
!fq_aer is determined based on kohler theory. 

      n_H2SO4 = TOMAS_H2SO4_setSpec(numTr)
      n_DMS = DMS_setSpec(numTr) ! note duplicate with Koch
      n_SO2 = SO2_setSpec(numTr)
#ifndef TRACERS_AEROSOLS_SOA
      n_SOAgas = TOMAS_SOAgas_setSpec(numTr)
#endif  /* TRACERS_AEROSOLS_SOA */
      n_H2O2_s = H2O2_s_setSpec(numTr) ! duplicate with Koch
      n_NH3 = NH3_setSpec(numTr) ! duplicate with nitrate
      n_NH4 = NH4_setSpec(numTr) ! duplicate with nitrate

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
      n_ASO4(:)  = TOMAS_setSpec(TOMAS_ASO4_setSpec,  numTr, nbins)
      n_ANACL(:) = TOMAS_setSpec(TOMAS_ANACL_setSpec, numTr, nbins)
      n_AECOB(:) = TOMAS_setSpec(TOMAS_AECOB_setSpec, numTr, nbins)
      n_AECIL(:) = TOMAS_setSpec(TOMAS_AECIL_setSpec, numTr, nbins)
      n_AOCOB(:) = TOMAS_setSpec(TOMAS_AOCOB_setSpec, numTr, nbins)
      n_AOCIL(:) = TOMAS_setSpec(TOMAS_AOCOB_setSpec, numTr, nbins)
      n_ADUST(:) = TOMAS_setSpec(TOMAS_ADUST_setSpec, numTr, nbins)
      n_ANUM(:) = TOMAS_setSpec(TOMAS_ANUM_setSpec, numTr, nbins)
      n_AH2O(:)  = TOMAS_setSpec(TOMAS_AH2O_setSpec,  numTr, nbins)
#endif /* TRACERS_AMP */

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
        to_conc(n) = 0
#ifdef TRACERS_SPECIAL_Shindell
C**** Aerosol tracer output should be mass mixing ratio
        select case (tr_wd_TYPE(n))
        case (nGAS)
          to_volume_MixRat(n) = 1 !gas output to volume mixing ratio
        case (nPART)
          to_volume_MixRat(n) = 0 !aerosol output to mass mixing ratio
        case (nWATER)
          to_volume_MixRat(n) = 0 !water output to mass mixing ratio
        case default
          to_volume_MixRat(n) = 0 !default output to mass mixing ratio
        end select
#endif
#if defined(TRACERS_GASEXCH_ocean_CO2) || defined(TRACERS_GASEXCH_land_CO2)
        to_volume_MixRat(n) = 1 !gas output to volume mixing ratio
#endif
#endif /* TRACERS_ON */

      end do
      return

      contains

      integer function Air_setSpec(n) result(n_Air)
      integer, intent(inout) :: n
      n = n + 1
      n_Air = n
      call set_ntm_power(n, -2)
      call set_tr_mm(n, mair)
      end function Air_setSpec

      integer function CO2n_setSpec(n) result(n_CO2n)
      integer, intent(inout) :: n
      n = n + 1
      n_CO2n = n
      call set_ntm_power(n, -6)
      call set_tr_mm(n, 44.d0)  !grams
      call set_t_qlimit(n,  .false.)
      call set_ntsurfsrc(n,  1)
      call set_needtrs(n, .true.)
      end function CO2n_setSpec

      integer function CFCn_setSpec(n) result(n_CFCn)
      integer, intent(inout) :: n
      n = n + 1
      n_CFCn = n
      call set_ntm_power(n, -12)
!     !call set_tr_mm(n, 136.0d0    !NCAR value
      call set_tr_mm(n, 137.37d0) !note units are in gr
      call set_ntsurfsrc(n,  1)
      call set_needtrs(n, .true.)
      end function CFCn_setSpec

      integer function SF6_setSpec(n) result(n_SF6)
      integer, intent(inout) :: n
      n = n + 1
      n_SF6 = n
      call set_ntm_power(n, -14)
      call set_tr_mm(n, 146.01d0)
      call set_ntsurfsrc(n,  1)
      end function SF6_setSpec

      integer function Rn222_setSpec(n) result(n_Rn222)
      integer, intent(inout) :: n
      n = n + 1
      n_Rn222 = n
      call set_ntm_power(n, -21)
      call set_tr_mm(n, 222.d0)
      call set_trdecay(n,  2.1d-6)
      call set_ntsurfsrc(n,  1)
      end function Rn222_setSpec

      integer function CO2_setSpec(n) result(n_CO2)
      integer, intent(inout) :: n
      n = n + 1
      n_CO2 = n
      call set_ntm_power(n, -6)
      call set_tr_mm(n, 44.d0)
      call set_t_qlimit(n,  .false.)
      call set_ntsurfsrc(n,  6)
      end function CO2_setSpec

      integer function N2O_setSpec(n) result(n_N2O)
      integer, intent(inout) :: n
      n = n + 1
      n_N2O = n
      call set_ntm_power(n, -9)
      call set_tr_mm(n, 44.d0)
#ifdef TRACERS_SPECIAL_Lerner
      call set_ntsurfsrc(n,  1)
      n_mptable(n) = 1
      tcscale(n_MPtable(n)) = 1.
#endif
      end function N2O_setSpec

      integer function CFC11_setSpec(n) result(n_CFC11)
      integer, intent(inout) :: n
      n = n + 1
      n_CFC11 = n
      call set_ntm_power(n, -12)
      call set_tr_mm(n, 137.4d0)
      call set_ntsurfsrc(n,  1)
#ifdef TRACERS_SPECIAL_Lerner
      n_mptable(n) = 2
      tcscale(n_MPtable(n)) = 1.
#endif
      end function CFC11_setSpec

      integer function C_14O2_setSpec(n) result(n_C_14O2)
      integer, intent(inout) :: n
      n = n + 1
      n_C_14O2 = n
      call set_ntm_power(n, -18)
      call set_tr_mm(n, 46.d0)
      call set_ntsurfsrc(n,  1)
      end function C_14O2_setSpec

      integer function CH4_setSpec(n) result(n_CH4)
      integer, intent(inout) :: n
      n = n + 1
      n_CH4 = n
      call set_tr_mm(n, 16.d0)
#ifdef TRACERS_SPECIAL_Lerner
      call set_ntsurfsrc(n,  14)
      call set_ntm_power(n, -9)
      n_mptable(n) = 3
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
      end function CH4_setSpec

      integer function O3_setSpec(n) result(n_O3)
      integer, intent(inout) :: n
      n = n + 1
      n_O3 = n
      call set_ntm_power(n, -8)
      call set_tr_mm(n, 48.d0)
      call set_ntsurfsrc(n,  1)
#ifdef TRACERS_SPECIAL_Lerner
C**** Get solar variability coefficient from namelist if it exits
      dsol = 0.
      call sync_param("dsol",dsol)
#endif
      end function O3_setSpec

      integer function SF6_c_setSpec(n) result(n_SF6_c)
      integer, intent(inout) :: n
      n = n + 1
      n_SF6_c = n
      call set_ntm_power(n, -14)
      call set_tr_mm(n, 146.01d0)
      call set_ntsurfsrc(n,  1)
      end function SF6_c_setSpec

      integer function Water_setSpec(n) result(n_Water)
      integer, intent(inout) :: n
      n = n + 1
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
      end function Water_setSpec

#ifdef TRACERS_SPECIAL_O18
      integer function H2O18_setSpec(n) result(n_H2O18)
      use TracerConstants_mod, only: H2O18
      integer, intent(inout) :: n
      n = n + 1
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
      end function H2O18_setSpec

      integer function HDO_setSpec(n) result(n_HDO)
      integer, intent(inout) :: n
      n = n + 1
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
      end function HDO_setSpec

      integer function HTO_setSpec(n) result(n_HTO)
      integer, intent(inout) :: n
      n = n + 1
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
      end function HTO_setSpec

      integer function H2O17_setSpec(n) result(n_H2O17)
      integer, intent(inout) :: n
      n = n + 1
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
      end function H2O17_setSpec
#endif  /* TRACERS_SPECIAL_O18 */

      integer function Ox_setSpec(n) result(n_Ox)
      integer, intent(inout) :: n
      n = n + 1
      n_Ox = n
      call set_ntm_power(n, -8)
      call set_tr_mm(n, 48.d0)
#ifdef TRACERS_DRYDEP
      call set_F0(n,  1.4d0)
      call set_HSTAR(n,  1.d-2)
#endif
      end function Ox_setSpec

      integer function NOx_setSpec(n) result(n_NOx)
      integer, intent(inout) :: n
      n = n + 1
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
      end function NOx_setSpec

      integer function N2O5_setSpec(n) result(n_N2O5)
      integer, intent(inout) :: n
      n = n + 1
      n_N2O5 = n
      call set_ntm_power(n, -12)
      call set_tr_mm(n, 108.02d0)
      end function N2O5_setSpec

#ifdef TRACERS_SPECIAL_Shindell
      integer function ClOx_setSpec(n) result(n_ClOx)
      integer, intent(inout) :: n
      n = n + 1
      n_ClOx = n
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 51.5d0)
C     Interpolate ClOx altitude-dependence to model resolution:
      CALL LOGPINT(LCOalt,PCOalt,ClOxaltIN,LM,PMIDL00,ClOxalt,.true.)
      end function ClOx_setSpec

      integer function BrOx_setSpec(n) result(n_BrOx)
      integer, intent(inout) :: n
      n = n + 1
      n_BrOx = n
      call set_ntm_power(n, -14)
      call set_tr_mm(n, 95.9d0)
C     Interpolate BrOx altitude-dependence to model resolution:
      CALL LOGPINT(LCOalt,PCOalt,BrOxaltIN,LM,PMIDL00,BrOxalt,.true.)
      end function BrOx_setSpec

      integer function HCl_setSpec(n) result(n_HCl)
      integer, intent(inout) :: n
      n = n + 1
      n_HCl = n
      call set_ntm_power(n, -10)
      call set_tr_mm(n, 36.5d0)
C     Interpolate HCl altitude-dependence to model resolution:
      CALL LOGPINT(LCOalt,PCOalt,HClaltIN,LM,PMIDL00,HClalt,.true.)
      end function HCl_setSpec

      integer function ClONO2_setSpec(n) result(n_ClONO2)
      integer, intent(inout) :: n
      n = n + 1
      n_ClONO2 = n
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 97.5d0)
C     Interpolate ClONO2 altitude-dependence to model resolution:
      CALL
     &    LOGPINT(LCOalt,PCOalt,ClONO2altIN,LM,PMIDL00,ClONO2alt,.true.)
      end function ClONO2_setSpec

      integer function HOCl_setSpec(n) result(n_HOCl)
      integer, intent(inout) :: n
      n = n + 1
      n_HOCl = n
      call set_ntm_power(n, -12)
      call set_tr_mm(n, 52.5d0)
      end function HOCl_setSpec
#endif /* TRACERS_SPECIAL_Shindell */ 

      integer function HBr_setSpec(n) result(n_HBr)
      integer, intent(inout) :: n
      n = n + 1
      n_HBr = n
      call set_ntm_power(n, -14)
      call set_tr_mm(n, 80.9d0)
      end function HBr_setSpec
      
      integer function HOBr_setSpec(n) result(n_HOBr)
      integer, intent(inout) :: n
      n = n + 1
      n_HOBr = n
      call set_ntm_power(n, -14)
      call set_tr_mm(n, 96.9d0)
      end function HOBr_setSpec

      integer function BrONO2_setSpec(n) result(n_BrONO2)
      integer, intent(inout) :: n
      n = n + 1
      n_BrONO2 = n
      call set_ntm_power(n, -14)
      call set_tr_mm(n, 141.9d0)
      end function BrONO2_setSpec

      integer function CFC_setSpec(n) result(n_CFC)
      integer, intent(inout) :: n
      n = n + 1
      n_CFC = n
      call set_ntm_power(n, -12)
      call set_tr_mm(n, 137.4d0) !CFC11
      end function CFC_setSpec


        
      integer function HNO3_setSpec(n) result(n_HNO3)
      integer, intent(inout) :: n
      n = n + 1
      n_HNO3 = n
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 63.018d0)
      call set_tr_RKD(n, 2.073d3 ) ! in mole/J = 2.1d5 mole/(L atm)
#ifdef TRACERS_DRYDEP
      call set_HSTAR(n, 1.d14)
#endif
      end function HNO3_setSpec

      integer function H2O2_setSpec(n) result(n_H2O2)
      integer, intent(inout) :: n
      n = n + 1
      n_H2O2 = n
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 34.016d0)
      call set_tr_RKD(n, 9.869d2    ) ! in mole/J = 1.d5 mole/(L atm)
      call set_tr_DHD(n, -5.52288d4 ) ! in J/mole = -13.2 kcal/mole.
#ifdef TRACERS_DRYDEP
      call set_HSTAR(n, tr_RKD(n)*convert_HSTAR)
      call set_F0(n,  1.d0)
#endif
      end function H2O2_setSpec

      integer function GLT_setSpec(n) result(n_GLT)
      integer, intent(inout) :: n
      n = n + 1
      n_GLT = n
      call set_ntm_power(n, -11)
      call set_tr_mm(n, mair)
      end function GLT_setSpec

      integer function stratOx_setSpec(n) result(n_stratOx)
      integer, intent(inout) :: n
      n = n + 1
      n_stratOx = n
! assumes initial Ox conditions read in for Ox tracer
      call set_ntm_power(n, -8)
      call set_tr_mm(n, 48.d0)
#ifdef TRACERS_DRYDEP
      call set_F0(n,  1.0d0)
      call set_HSTAR(n,  1.d-2)
#endif
      end function stratOx_setSpec


      integer function codirect_setSpec(n) result(n_codirect)
      integer, intent(inout) :: n
      n = n + 1
      n_codirect = n
      call set_ntm_power(n, -8)
      call set_tr_mm(n, 28.01d0)
      call set_trdecay(n,  2.31482d-7) ! 1/(50 days)
! not a radiactive decay, but functionally identical
      end function codirect_setSpec

      integer function CH3OOH_setSpec(n) result(n_CH3OOH)
      integer, intent(inout) :: n
      n = n + 1
      n_CH3OOH = n
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 48.042d0)
#ifdef TRACERS_DRYDEP
      call set_HSTAR(n,  3.d2)
#endif
      end function CH3OOH_setSpec

      integer function HCHO_setSpec(n) result(n_HCHO)
      integer, intent(inout) :: n
      n = n + 1
      n_HCHO = n
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 30.026d0)
      call set_tr_RKD(n, 6.218d1 ) ! mole/J = 6.3d3 mole/(L atm)
#ifdef TRACERS_DRYDEP
      call set_HSTAR(n, 6.d3)
#endif
      end function HCHO_setSpec

      integer function HO2NO2_setSpec(n) result(n_HO2NO2)
      integer, intent(inout) :: n
      n = n + 1
      n_HO2NO2 = n
      call set_ntm_power(n, -12)
      call set_tr_mm(n, 79.018d0)
      end function HO2NO2_setSpec

      integer function CO_setSpec(n) result(n_CO)
      integer, intent(inout) :: n
      n = n + 1
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
      end function CO_setSpec

      integer function PAN_setSpec(n) result(n_PAN)
      integer, intent(inout) :: n
      n = n + 1
      n_PAN = n
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 121.054d0) ! assuming CH3COOONO2 = PAN)
#ifdef TRACERS_DRYDEP
      call set_HSTAR(n,  3.6d0)
#endif
      end function PAN_setSpec

      integer function Isoprene_setSpec(n) result(n_Isoprene)
      integer, intent(inout) :: n
      n = n + 1
      n_Isoprene = n
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 60.05d0) ! i.e. 5 carbons
#ifdef TRACERS_DRYDEP
      call set_HSTAR(n,  1.3d-2)
#endif
      end function Isoprene_setSpec

      integer function AlkylNit_setSpec(n) result(n_AlkylNit)
      integer, intent(inout) :: n
      n = n + 1
      n_AlkylNit = n
      call set_ntm_power(n, -11)
      call set_tr_mm(n, mair)   !unknown molecular weight, so use air and make
                                ! note in the diagnostics write-out...
      end function AlkylNit_setSpec

      integer function Alkenes_setSpec(n) result(n_Alkenes)
      integer, intent(inout) :: n
      n = n + 1
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
      end function Alkenes_setSpec

      integer function Paraffin_setSpec(n) result(n_Paraffin)
      integer, intent(inout) :: n
      n = n + 1
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
      end function Paraffin_setSpec

      integer function Terpenes_setSpec(n) result(n_Terpenes)
      integer, intent(inout) :: n
      n = n + 1
      n_Terpenes = n
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 120.10d0) ! i.e. 10 carbons
#ifdef TRACERS_DRYDEP
      call set_HSTAR(n,  1.3d-2)
#endif
      end function Terpenes_setSpec

#ifdef TRACERS_AEROSOLS_SOA
      integer function isopp1g_setSpec(n) result(n_isopp1g)
      integer, intent(inout) :: n
      n = n + 1
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
      end function isopp1g_setSpec

      integer function isopp1a_setSpec(n) result(n_isopp1a)
      integer, intent(inout) :: n
      n = n + 1
      n_isopp1a = n
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 15.6d0)
      call set_trpdens(n, 1.5d3) !kg/m3
      call set_trradius(n, 3.d-7) !m
      call set_fq_aer(n, 0.8d0) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, nPART)
      end function isopp1a_setSpec

      integer function isopp2g_setSpec(n) result(n_isopp2g)
      integer, intent(inout) :: n
      n = n + 1
      n_isopp2g = n
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 15.6d0)
      call set_tr_RKD(n, 1.d4 / convert_HSTAR ) !Henry; from mole/(L atm) to mole/J
      call set_tr_DHD(n, -12.d0 * gasc        ) !Henry temp dependence (J/mole), Chung and Seinfeld, 2002
      call set_tr_wd_type(n, ngas)
#ifdef TRACERS_DRYDEP
      call set_HSTAR(n,  tr_RKD(n) * convert_HSTAR)
#endif
      end function isopp2g_setSpec

      integer function isopp2a_setSpec(n) result(n_isopp2a)
      integer, intent(inout) :: n
      n = n + 1
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
      end function isopp2a_setSpec

      integer function apinp1g_setSpec(n) result(n_apinp1g)
      integer, intent(inout) :: n
      n = n + 1
      n_apinp1g = n
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 15.6d0)
      call set_tr_RKD(n, 1.d4 / convert_HSTAR ) !Henry; from mole/(L atm) to mole/J
      call set_tr_DHD(n, -12.d0 * gasc        ) !Henry temp dependence (J/mole), Chung and Seinfeld, 2002
      call set_tr_wd_type(n, ngas)
#ifdef TRACERS_DRYDEP
      call set_HSTAR(n,  tr_RKD(n) * convert_HSTAR)
#endif
      end function apinp1g_setSpec

      integer function apinp1a_setSpec(n) result(n_apinp1a)
      integer, intent(inout) :: n
      n = n + 1
      n_apinp1a = n
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 15.6d0)
      call set_trpdens(n, 1.5d3) !kg/m3
      call set_trradius(n, 3.d-7) !m
      call set_fq_aer(n, 0.8d0) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, nPART)
      end function apinp1a_setSpec

      integer function apinp2g_setSpec(n) result(n_apinp2g)
      integer, intent(inout) :: n
      n = n + 1
      n_apinp2g = n
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 15.6d0)
      call set_tr_RKD(n, 1.d4 / convert_HSTAR ) !Henry; from mole/(L atm) to mole/J
      call set_tr_DHD(n, -12.d0 * gasc        ) !Henry temp dependence (J/mole), Chung and Seinfeld, 2002
      call set_tr_wd_type(n, ngas)
#ifdef TRACERS_DRYDEP
      call set_HSTAR(n,  tr_RKD(n) * convert_HSTAR)
#endif
      end function apinp2g_setSpec

      integer function apinp2a_setSpec(n) result(n_apinp2a)
      integer, intent(inout) :: n
      n = n + 1
      n_apinp2a = n
      n_soa_e = n_apinp2a       !the last from the soa species
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 15.6d0)
      call set_trpdens(n, 1.5d3) !kg/m3
      call set_trradius(n, 3.d-7) !m
      call set_fq_aer(n, 0.8d0) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, nPART)
      end function apinp2a_setSpec
#endif  /* TRACERS_AEROSOLS_SOA */

      integer function DMS_setSpec(n) result(n_DMS)
      integer, intent(inout) :: n
      n = n + 1
      n_DMS = n
      call set_ntm_power(n, -12)
!     the ocean source of DMS is actually interactive and therefore should
!     not count for ntsurfsrc....
      call set_ntsurfsrc(n,  0) ! ocean DMS concentration
      call set_ntisurfsrc(n, 1)
      call set_tr_mm(n, 62.d+0)
      call set_needtrs(n, .true.)
      end function DMS_setSpec

      integer function MSA_setSpec(n) result(n_MSA)
      integer, intent(inout) :: n
      n = n + 1
      n_MSA = n
      call set_ntm_power(n, -13)
      call set_tr_mm(n, 96.d+0) !(H2O2 34;SO2 64)
      call set_trpdens(n, 1.7d3) !kg/m3 this is sulfate value
      call set_trradius(n, 5.d-7 ) !m (SO4 3;BC 1;OC 3)
      call set_fq_aer(n, 1.0d0   ) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, npart)
      end function MSA_setSpec

      integer function SO2_setSpec(n) result(n_SO2)
      integer, intent(inout) :: n
      n = n + 1
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
      end function SO2_setSpec

      integer function SO4_setSpec(n) result(n_SO4)
      integer, intent(inout) :: n
      n = n + 1
      n_SO4 = n
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 96.d+0)
      call set_trpdens(n, 1.7d3) !kg/m3 this is sulfate value
      call set_trradius(n, 3.d-7 ) !m
      call set_fq_aer(n, 1.d0   ) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, npart)
      end function SO4_setSpec

      integer function SO4_d1_setSpec(n) result(n_SO4_d1)
      integer, intent(inout) :: n
      n = n + 1
      n_SO4_d1 = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 96.d0)  !!!! Sulfat
      call set_trpdens(n, 2.5d3) !kg/m3 this is clay density
      call set_trradius(n, 0.75D-06 ) !m
      call set_fq_aer(n, 1.d0   ) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, npart)
      end function SO4_d1_setSpec

      integer function SO4_d2_setSpec(n) result(n_SO4_d2)
      integer, intent(inout) :: n
      n = n + 1
      n_SO4_d2 = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 96.d0)
      call set_trpdens(n, 2.65d3) !kg/m3 this is silt1 value
      call set_trradius(n, 2.2D-06 ) !m
      call set_fq_aer(n, 1.d0   ) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, npart)
      end function SO4_d2_setSpec

      integer function SO4_d3_setSpec(n) result(n_SO4_d3)
      integer, intent(inout) :: n
      n = n + 1
      n_SO4_d3 = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 96.d0)
      call set_trpdens(n, 2.65d3) !this is silt2 value
      call set_trradius(n, 4.4D-06 ) !m this is silt2 value
      call set_fq_aer(n, 1.d0   ) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, npart)
      end function SO4_d3_setSpec

      integer function N_d1_setSpec(n) result(n_N_d1)
      integer, intent(inout) :: n
      n = n + 1
      n_N_d1 = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 62.d+0) ! NO3
      call set_trpdens(n, 2.5d3) !kg/m3 this is clay density
      call set_trradius(n, 0.75D-06 ) !m
      call set_fq_aer(n, 1.d0  ) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, npart)
      end function N_d1_setSpec

      integer function N_d2_setSpec(n) result(n_N_d2)
      integer, intent(inout) :: n
      n = n + 1
      n_N_d2 = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 62.d+0)
      call set_trpdens(n, 2.65d3) !kg/m3 this is silt1 value
      call set_trradius(n, 2.2D-06 ) !m
      call set_fq_aer(n, 1.d0  ) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, npart)
      end function N_d2_setSpec

      integer function N_d3_setSpec(n) result(n_N_d3)
      integer, intent(inout) :: n
      n = n + 1
      n_N_d3 = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 62.d0)
      call set_trpdens(n, 2.65d3) !this is silt2 value
      call set_trradius(n, 4.4D-06 ) !m this is silt2 value
      call set_fq_aer(n, 1.d0  ) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, npart)
      end function N_d3_setSpec

      integer function BCII_setSpec(n) result(n_BCII)
      integer, intent(inout) :: n
      n = n + 1
      n_BCII = n
      call set_ntm_power(n, -12)
      call set_tr_mm(n, 12.d0)
      call set_trpdens(n, 1.3d3) !kg/m3
      call set_trradius(n, 1.d-7 ) !m
      call set_fq_aer(n, 0.0d0   ) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, npart)
      end function BCII_setSpec

      integer function BCIA_setSpec(n) result(n_BCIA)
      integer, intent(inout) :: n
      n = n + 1
      n_BCIA = n
      call set_ntm_power(n, -12)
      call set_tr_mm(n, 12.d0)
      call set_trpdens(n, 1.3d3) !kg/m3
      call set_trradius(n, 1.d-7 ) !m
      call set_fq_aer(n, 1.d0   ) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, npart)
      end function BCIA_setSpec

      integer function BCB_setSpec(n) result(n_BCB)
      integer, intent(inout) :: n
      n = n + 1
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
      end function BCB_setSpec

      integer function OCII_setSpec(n) result(n_OCII)
      integer, intent(inout) :: n
      n = n + 1
      n_OCII = n
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 15.6d0)
      call set_trpdens(n, 1.5d3) !kg/m3
      call set_trradius(n, 3.d-7 ) !m
      call set_fq_aer(n, 0.0d0   ) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, npart)
      end function OCII_setSpec

      integer function OCIA_setSpec(n) result(n_OCIA)
      integer, intent(inout) :: n
      n = n + 1
      n_OCIA = n
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 15.6d0)
      call set_trpdens(n, 1.5d3) !kg/m3
      call set_trradius(n, 3.d-7 ) !m
      call set_fq_aer(n, 1.d0   ) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, npart)
      end function OCIA_setSpec

      integer function OCB_setSpec(n) result(n_OCB)
      integer, intent(inout) :: n
      n = n + 1
      n_OCB = n
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)
      call sync_param("OCB_om2oc",om2oc(n_OCB))
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
      end function OCB_setSpec

      integer function Be7_setSpec(n) result(n_Be7)
      integer, intent(inout) :: n
      n = n + 1
      n_Be7 = n
      call set_ntm_power(n, -23) ! power of ten for tracer
      call set_tr_mm(n, 7.d0)
      call set_trdecay(n,  1.51d-7)
      call set_trpdens(n, 1.7d3) !kg/m3 this is SO4 value
      call set_trradius(n, 1.d-7  ) !appropriate for stratosphere
      call set_fq_aer(n, 1.d0   ) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, npart) ! same as SO4
      end function Be7_setSpec

      integer function Be10_setSpec(n) result(n_Be10)
      integer, intent(inout) :: n
      n = n + 1
      n_Be10 = n
      call set_ntm_power(n, -23)
      call set_tr_mm(n, 10.d0)
      call set_trpdens(n, 1.7d3) !kg/m3 this is SO4 value
      call set_trradius(n, 1.d-7  ) !appropriate for stratosphere
      call set_fq_aer(n, 1.d0   ) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, npart) ! same as SO4
      end function Be10_setSpec

      integer function Pb210_setSpec(n) result(n_Pb210)
      integer, intent(inout) :: n
      n = n + 1
      n_Pb210 = n
      call set_ntm_power(n, -23)
      call set_tr_mm(n, 210.d0)
      call set_trdecay(n,  9.85d-10)
      call set_trpdens(n, 1.7d3) !kg/m3 this is SO4 value
      call set_trradius(n, 3.d-7  ) !again S04 value
      call set_fq_aer(n, 1.d0   ) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, npart) ! same as SO4
      end function Pb210_setSpec

      integer function H2O2_s_setSpec(n) result(n_H2O2_s)
      integer, intent(inout) :: n
      n = n + 1
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
      end function H2O2_s_setSpec

      integer function seasalt1_setSpec(n) result(n_seasalt1)
      integer, intent(inout) :: n
      n = n + 1
      n_seasalt1 = n
      call set_ntsurfsrc(n,  0) ! ocean bubbles
      call set_ntisurfsrc(n, 1)
      call set_ntm_power(n, -10)
      call set_tr_mm(n, 75.d0)  !Na x 3.256
      call set_trpdens(n, 2.2d3) !kg/m3 This is for non-hydrated
      call set_trradius(n, 4.4d-7 ) ! This is non-hydrated
      call set_fq_aer(n, 1.0d0   ) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, npart)
      end function seasalt1_setSpec

      integer function seasalt2_setSpec(n) result(n_seasalt2)
      integer, intent(inout) :: n
      n = n + 1
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
      end function seasalt2_setSpec

      integer function OCocean_setSpec(n) result(n_OCocean)
      integer, intent(inout) :: n
      n = n + 1
      n_OCocean = n
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 15.6d0)
      call set_trpdens(n, 1.5d3) !kg/m3
      call set_trradius(n, 4.4d-7) !m, same as seasalt1
      call set_fq_aer(n, 1.0d0) ! same as seasalt
      call set_tr_wd_type(n, nPART)
      end function OCocean_setSpec

      integer function NH3_setSpec(n) result(n_NH3)
      integer, intent(inout) :: n
      n = n + 1
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
      end function NH3_setSpec

      integer function NH4_setSpec(n) result(n_NH4)
      integer, intent(inout) :: n
      n = n + 1
      n_NH4 = n
      call set_ntsurfsrc(n,  0)
      call set_ntm_power(n, -10)
      call set_tr_mm(n, 18.d0)
      call set_trpdens(n, 1.7d3)
      call set_trradius(n, 3.d-7)
      call set_fq_aer(n, 1.0d0   ) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, npart)
      end function NH4_setSpec

      integer function NO3p_setSpec(n) result(n_NO3p)
      integer, intent(inout) :: n
      n = n + 1
      n_NO3p = n
      call set_ntsurfsrc(n,  0)
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 62.d0)
      call set_trpdens(n, 1.7d3)
      call set_trradius(n, 3.d-7)
      call set_fq_aer(n, 1.0d0   ) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, npart)
      end function NO3p_setSpec

      integer function clay_setSpec(n) result(n_clay)
      integer, intent(inout) :: n
      n = n + 1
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
      end function clay_setSpec

      integer function silt1_setSpec(n) result(n_silt1)
      integer, intent(inout) :: n
      n = n + 1
      n_silt1=n
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
      end function silt1_setSpec

      integer function silt2_setSpec(n) result(n_silt2)
      integer, intent(inout) :: n
      n = n + 1
      n_silt2=n
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
      end function silt2_setSpec

      integer function silt3_setSpec(n) result(n_silt3)
      integer, intent(inout) :: n
      n = n + 1
      n_silt3=n
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
      end function silt3_setSpec

      integer function silt4_setSpec(n) result(n_silt4)
      integer, intent(inout) :: n
      n = n + 1
      n_silt4=n
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
      end function silt4_setSpec

#ifdef TRACERS_MINERALS
      integer function clayilli_setSpec(n) result(n_clayilli)
      integer, intent(inout) :: n
      n = n + 1
      n_clayilli=n
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

      end function clayilli_setSpec

      integer function claykaol_setSpec(n) result(n_claykaol)
      integer, intent(inout) :: n
      n = n + 1
      n_claykaol=n
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

      end function claykaol_setSpec

      integer function claysmec_setSpec(n) result(n_claysmec)
      integer, intent(inout) :: n
      n = n + 1
      n_claysmec=n
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

      end function claysmec_setSpec

      integer function claycalc_setSpec(n) result(n_claycalc)
      integer, intent(inout) :: n
      n = n + 1
      n_claycalc=n
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
      end function claycalc_setSpec

      integer function clayquar_setSpec(n) result(n_clayquar)
      integer, intent(inout) :: n
      n = n + 1
      n_clayquar=n
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
      end function clayquar_setSpec

      integer function sil1quar_setSpec(n) result(n_sil1quar)
      integer, intent(inout) :: n
      n = n + 1
      n_sil1quar=n
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
      end function sil1quar_setSpec

      integer function sil1feld_setSpec(n) result(n_sil1feld)
      integer, intent(inout) :: n
      n = n + 1
      n_sil1feld=n
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
      end function sil1feld_setSpec

      integer function sil1calc_setSpec(n) result(n_sil1calc)
      integer, intent(inout) :: n
      n = n + 1
      n_sil1calc=n
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
      end function sil1calc_setSpec

      integer function sil1hema_setSpec(n) result(n_sil1hema)
      integer, intent(inout) :: n
      n = n + 1
      n_sil1hema=n
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
      end function sil1hema_setSpec

      integer function sil1gyps_setSpec(n) result(n_sil1gyps)
      integer, intent(inout) :: n
      n = n + 1
      n_sil1gyps=n
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
      end function sil1gyps_setSpec

      integer function sil2quar_setSpec(n) result(n_sil2quar)
      integer, intent(inout) :: n
      n = n + 1
      n_sil2quar=n
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
      end function sil2quar_setSpec

      integer function sil2feld_setSpec(n) result(n_sil2feld)
      integer, intent(inout) :: n
      n = n + 1
      n_sil2feld=n
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
      end function sil2feld_setSpec

      integer function sil2calc_setSpec(n) result(n_sil2calc)
      integer, intent(inout) :: n
      n = n + 1
      n_sil2calc=n
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
      end function sil2calc_setSpec

      integer function sil2hema_setSpec(n) result(n_sil2hema)
      integer, intent(inout) :: n
      n = n + 1
      n_sil2hema=n
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
      end function sil2hema_setSpec

      integer function sil2gyps_setSpec(n) result(n_sil2gyps)
      integer, intent(inout) :: n
      n = n + 1
      n_sil2gyps=n
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
      end function sil2gyps_setSpec

      integer function sil3quar_setSpec(n) result(n_sil3quar)
      integer, intent(inout) :: n
      n = n + 1
      n_sil3quar=n
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
      end function sil3quar_setSpec

      integer function sil3feld_setSpec(n) result(n_sil3feld)
      integer, intent(inout) :: n
      n = n + 1
      n_sil3feld=n
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
      end function sil3feld_setSpec

      integer function sil3calc_setSpec(n) result(n_sil3calc)
      integer, intent(inout) :: n
      n = n + 1
      n_sil3calc=n
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
      end function sil3calc_setSpec

      integer function sil3hema_setSpec(n) result(n_sil3hema)
      integer, intent(inout) :: n
      n = n + 1
      n_sil3hema=n
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
      end function sil3hema_setSpec

      integer function sil3gyps_setSpec(n) result(n_sil3gyps)
      integer, intent(inout) :: n
      n = n + 1
      n_sil3gyps=n
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
      end function sil3gyps_setSpec

#endif  /* TRACERS_MINERALS */

#ifdef TRACERS_QUARZHEM
      integer function sil1quhe_setSpec(n) result(n_sil1quhe)
      integer, intent(inout) :: n
      n = n + 1
      n_sil1quhe=n
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
      end function sil1quhe_setSpec

      integer function sil2quhe_setSpec(n) result(n_sil2quhe)
      integer, intent(inout) :: n
      n = n + 1
      n_sil2quhe=n
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
      end function sil2quhe_setSpec

      integer function sil3quhe_setSpec(n) result(n_sil3quhe)
      integer, intent(inout) :: n
      n = n + 1
      n_sil3quhe=n
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
      end function sil3quhe_setSpec
#endif  /* TRACERS_QUARZHEM */

#ifdef TRACERS_AMP
      integer function H2SO4_setSpec(n) result(n_H2SO4)
      integer, intent(inout) :: n
      n = n + 1
      n_H2SO4 = n
          call set_ntm_power(n, -11)
          call set_ntsurfsrc(n,  0)
          call set_tr_mm(n, 98.d0)
          call set_trpdens(n, DENS_SULF)
          call set_trradius(n, DG_ACC * .5d-6)
          call set_fq_aer(n, SOLU_ACC)
          call set_tr_wd_type(n, npart)
      end function H2SO4_setSpec

      integer function M_NO3_setSpec(n) result(n_M_NO3)
      integer, intent(inout) :: n
      n = n + 1
      n_M_NO3 = n
      ntmAMPi=n                 ! always the first tracer in AMP
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 62.d0)
      call set_trpdens(n, 1.7d3)
      call set_trradius(n, 3.d-7 ) !m
      call set_fq_aer(n, 1.d0)  !fraction of aerosol that dissolves
      call set_tr_wd_type(n, npart)
      end function M_NO3_setSpec

      integer function M_NH4_setSpec(n) result(n_M_NH4)
      integer, intent(inout) :: n
      n = n + 1
      n_M_NH4 = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 18.d0)
      call set_trpdens(n, 1.7d3)
      call set_trradius(n, 3.d-7)
      call set_fq_aer(n, 1.d+0)
      call set_tr_wd_type(n, npart)
      end function M_NH4_setSpec

      integer function M_H2O_setSpec(n) result(n_M_H2O)
      integer, intent(inout) :: n
      n = n + 1
      n_M_H2O = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, mwat)
      call set_trpdens(n, 1.d3)
      call set_trradius(n, 3.d-7)
      call set_fq_aer(n, 1.d+0)
      call set_tr_wd_type(n, npart) !nWater
      end function M_H2O_setSpec

      integer function M_AKK_SU_setSpec(n) result(n_M_AKK_SU)
      integer, intent(inout) :: n
      n = n + 1
      n_M_AKK_SU = n
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 96.d0)
      call set_trpdens(n, DENS_SULF)
      call set_trradius(n, DG_AKK * .5d-6)
      call set_fq_aer(n,  SOLU_AKK)
      call set_tr_wd_type(n, npart)
      end function M_AKK_SU_setSpec

      integer function N_AKK_1_setSpec(n) result(n_N_AKK_1)
      integer, intent(inout) :: n
      n = n + 1
      n_N_AKK_1 = n
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 1.d+0)
      call set_trpdens(n, DENS_SULF)
      call set_trradius(n, DG_AKK * .5d-6)
      call set_fq_aer(n, SOLU_AKK)
      call set_tr_wd_type(n, npart)
      end function N_AKK_1_setSpec

      integer function M_ACC_SU_setSpec(n) result(n_M_ACC_SU)
      integer, intent(inout) :: n
      n = n + 1
      n_M_ACC_SU = n
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 96.d0)
      call set_trpdens(n, DENS_SULF)
      call set_trradius(n, DG_ACC * .5d-6)
      call set_fq_aer(n, SOLU_ACC)
      call set_tr_wd_type(n, npart)
      end function M_ACC_SU_setSpec

      integer function N_ACC_1_setSpec(n) result(n_N_ACC_1)
      integer, intent(inout) :: n
      n = n + 1
      n_N_ACC_1 = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 1.d+0)
      call set_trpdens(n, DENS_SULF)
      call set_trradius(n, DG_ACC * .5d-6)
      call set_fq_aer(n, SOLU_ACC)
      call set_tr_wd_type(n, npart)
      end function N_ACC_1_setSpec

      integer function M_DD1_SU_setSpec(n) result(n_M_DD1_SU)
      integer, intent(inout) :: n
      n = n + 1
      n_M_DD1_SU = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 96.d0)
      call set_trpdens(n, DENS_SULF)
      call set_trradius(n, DG_DD1 * .5d-6)
      call set_fq_aer(n, SOLU_DD1)
      call set_tr_wd_type(n, npart)
      end function M_DD1_SU_setSpec

      integer function M_DD1_DU_setSpec(n) result(n_M_DD1_DU)
      integer, intent(inout) :: n
      n = n + 1
      n_M_DD1_DU = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 1.d+0)
      call set_trpdens(n, DENS_DUST)
      call set_trradius(n, DG_DD1 * .5d-6)
      call set_fq_aer(n,  SOLU_DD1)
      call set_tr_wd_type(n, npart)
      end function M_DD1_DU_setSpec

      integer function N_DD1_1_setSpec(n) result(n_N_DD1_1)
      integer, intent(inout) :: n
      n = n + 1
      n_N_DD1_1 = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 1.d+0)
      call set_trpdens(n, DENS_DUST)
      call set_trradius(n, DG_DD1 * .5d-6)
      call set_fq_aer(n,  SOLU_DD1)
      call set_tr_wd_type(n, npart)
      end function N_DD1_1_setSpec

      integer function M_DS1_SU_setSpec(n) result(n_M_DS1_SU)
      integer, intent(inout) :: n
      n = n + 1
      n_M_DS1_SU = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 96.d0)
      call set_trpdens(n, DENS_SULF)
      call set_trradius(n, DG_DS1 * .5d-6)
      call set_fq_aer(n, SOLU_DS1)
      call set_tr_wd_type(n, npart)
      end function M_DS1_SU_setSpec

      integer function M_DS1_DU_setSpec(n) result(n_M_DS1_DU)
      integer, intent(inout) :: n
      n = n + 1
      n_M_DS1_DU = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 1.d+0)
      call set_trpdens(n, DENS_DUST)
      call set_trradius(n, DG_DS1 * .5d-6)
      call set_fq_aer(n, SOLU_DS1)
      call set_tr_wd_type(n, npart)
      end function M_DS1_DU_setSpec

      integer function N_DS1_1_setSpec(n) result(n_N_DS1_1)
      integer, intent(inout) :: n
      n = n + 1
      n_N_DS1_1= n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 1.d+0)
      call set_trpdens(n, DENS_DUST)
      call set_trradius(n, DG_DS1 * .5d-6)
      call set_fq_aer(n, SOLU_DS1)
      call set_tr_wd_type(n, npart)
      end function N_DS1_1_setSpec

      integer function M_DD2_SU_setSpec(n) result(n_M_DD2_SU)
      integer, intent(inout) :: n
      n = n + 1
      n_M_DD2_SU = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 96.d0)
      call set_trpdens(n, DENS_SULF)
      call set_trradius(n, DG_DD2 * .5d-6)
      call set_fq_aer(n, SOLU_DD2)
      call set_tr_wd_type(n, npart)
      end function M_DD2_SU_setSpec

      integer function M_DD2_DU_setSpec(n) result(n_M_DD2_DU)
      integer, intent(inout) :: n
      n = n + 1
      n_M_DD2_DU = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 1.d+0)
      call set_trpdens(n, DENS_DUST)
      call set_trradius(n, DG_DD2 * .5d-6)
      call set_fq_aer(n,  SOLU_DD2)
      call set_tr_wd_type(n, npart)
      end function M_DD2_DU_setSpec

      integer function N_DD2_1_setSpec(n) result(n_N_DD2_1)
      integer, intent(inout) :: n
      n = n + 1
      n_N_DD2_1 = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 1.d+0)
      call set_trpdens(n, DENS_DUST)
      call set_trradius(n, DG_DD2 * .5d-6)
      call set_fq_aer(n,  SOLU_DD2)
      call set_tr_wd_type(n, npart)
      end function N_DD2_1_setSpec

      integer function M_DS2_SU_setSpec(n) result(n_M_DS2_SU)
      integer, intent(inout) :: n
      n = n + 1
      n_M_DS2_SU = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 96.d0)
      call set_trpdens(n, DENS_SULF)
      call set_trradius(n, DG_DS2 * .5d-6)
      call set_fq_aer(n, SOLU_DS2)
      call set_tr_wd_type(n, npart)
      end function M_DS2_SU_setSpec

      integer function M_DS2_DU_setSpec(n) result(n_M_DS2_DU)
      integer, intent(inout) :: n
      n = n + 1
      n_M_DS2_DU = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 1.d+0)
      call set_trpdens(n, DENS_DUST)
      call set_trradius(n, DG_DS2 * .5d-6)
      call set_fq_aer(n, SOLU_DS2)
      call set_tr_wd_type(n, npart)
      end function M_DS2_DU_setSpec

      integer function N_DS2_1_setSpec(n) result(n_N_DS2_1)
      integer, intent(inout) :: n
      n = n + 1
      n_N_DS2_1 = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 1.d+0)
      call set_trpdens(n, DENS_DUST)
      call set_trradius(n, DG_DS2 * .5d-6)
      call set_fq_aer(n, SOLU_DS2)
      call set_tr_wd_type(n, npart)
      end function N_DS2_1_setSpec

      integer function M_SSA_SU_setSpec(n) result(n_M_SSA_SU)
      integer, intent(inout) :: n
      n = n + 1
      n_M_SSA_SU = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 96.d0)
      call set_trpdens(n, DENS_SULF)
      call set_trradius(n, DG_SSA * .5d-6)
      call set_fq_aer(n, SOLU_SSA)
      call set_tr_wd_type(n, npart)
      end function M_SSA_SU_setSpec

      integer function M_SSA_SS_setSpec(n) result(n_M_SSA_SS)
      integer, intent(inout) :: n
      n = n + 1
      n_M_SSA_SS = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 75.d0)
      call set_trpdens(n, DENS_SEAS)
      call set_trradius(n, DG_SSA * .5d-6)
      call set_fq_aer(n, SOLU_SSA)
      call set_tr_wd_type(n, npart)
      end function M_SSA_SS_setSpec

      integer function M_SSC_SS_setSpec(n) result(n_M_SSC_SS)
      integer, intent(inout) :: n
      n = n + 1
      n_M_SSC_SS = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 75.d0)
      call set_trpdens(n, DENS_SEAS)
      call set_trradius(n, DG_SSC * .5d-6)
      call set_fq_aer(n, SOLU_SSC)
      call set_tr_wd_type(n, npart)
      end function M_SSC_SS_setSpec

      integer function M_SSS_SS_setSpec(n) result(n_M_SSS_SS)
      integer, intent(inout) :: n
      n = n + 1
      n_M_SSS_SS = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 75.d0)
      call set_trpdens(n, DENS_SEAS)
      call set_trradius(n, DG_SSS * .5d-6)
      call set_fq_aer(n, SOLU_SSS)
      call set_tr_wd_type(n, npart)
      end function M_SSS_SS_setSpec

      integer function M_SSS_SU_setSpec(n) result(n_M_SSS_SU)
      integer, intent(inout) :: n
      n = n + 1
      n_M_SSS_SU = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 75.d0)
      call set_trpdens(n, DENS_SULF)
      call set_trradius(n, DG_SSS * .5d-6)
      call set_fq_aer(n, SOLU_SSS)
      call set_tr_wd_type(n, npart)
      end function M_SSS_SU_setSpec

      integer function M_OCC_SU_setSpec(n) result(n_M_OCC_SU)
      integer, intent(inout) :: n
      n = n + 1
      n_M_OCC_SU = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 96.d0)
      call set_trpdens(n, DENS_SULF)
      call set_trradius(n, DG_OCC * .5d-6)
      call set_fq_aer(n, SOLU_OCC)
      call set_tr_wd_type(n, npart)
      end function M_OCC_SU_setSpec

      integer function M_OCC_OC_setSpec(n) result(n_M_OCC_OC)
      integer, intent(inout) :: n
      n = n + 1
      n_M_OCC_OC = n
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 15.6d0)
      call set_trpdens(n, DENS_OCAR)
      call set_trradius(n, DG_OCC * .5d-6)
      call set_fq_aer(n,  SOLU_OCC)
      call set_tr_wd_type(n, npart)
      end function M_OCC_OC_setSpec

      integer function N_OCC_1_setSpec(n) result(n_N_OCC_1)
      integer, intent(inout) :: n
      n = n + 1
      n_N_OCC_1 = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 1.d+0)
      call set_trpdens(n, DENS_OCAR)
      call set_trradius(n, DG_OCC * .5d-6)
      call set_fq_aer(n,  SOLU_OCC)
      call set_tr_wd_type(n, npart)
      end function N_OCC_1_setSpec

      integer function M_BC1_SU_setSpec(n) result(n_M_BC1_SU)
      integer, intent(inout) :: n
      n = n + 1
      n_M_BC1_SU = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 96.d0)
      call set_trpdens(n, DENS_SULF)
      call set_trradius(n, DG_BC1 * .5d-6)
      call set_fq_aer(n, SOLU_BC1)
      call set_tr_wd_type(n, npart)
      end function M_BC1_SU_setSpec

      integer function M_BC1_BC_setSpec(n) result(n_M_BC1_BC)
      integer, intent(inout) :: n
      n = n + 1
      n_M_BC1_BC = n
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 12.d0)
      call set_trpdens(n, DENS_BCAR)
      call set_trradius(n, DG_BC1 * .5d-6)
      call set_fq_aer(n,  SOLU_BC1)
      call set_tr_wd_type(n, npart)
      end function M_BC1_BC_setSpec

      integer function N_BC1_1_setSpec(n) result(n_N_BC1_1)
      integer, intent(inout) :: n
      n = n + 1
      n_N_BC1_1 = n
          call set_ntm_power(n, -11)
          call set_ntsurfsrc(n,  0)
          call set_tr_mm(n, 1.d+0)
          call set_trpdens(n, DENS_BCAR)
          call set_trradius(n, DG_BC1 * .5d-6)
          call set_fq_aer(n,  SOLU_BC1)
          call set_tr_wd_type(n, npart)
      end function N_BC1_1_setSpec

      integer function M_BC2_SU_setSpec(n) result(n_M_BC2_SU)
      integer, intent(inout) :: n
      n = n + 1
      n_M_BC2_SU = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 96.d0)
      call set_trpdens(n, DENS_SULF)
      call set_trradius(n, DG_BC2 * .5d-6)
      call set_fq_aer(n, SOLU_BC2)
      call set_tr_wd_type(n, npart)
      end function M_BC2_SU_setSpec

      integer function M_BC2_BC_setSpec(n) result(n_M_BC2_BC)
      integer, intent(inout) :: n
      n = n + 1
      n_M_BC2_BC = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 12.d0)
      call set_trpdens(n, DENS_BCAR)
      call set_trradius(n, DG_BC2 * .5d-6)
      call set_fq_aer(n, SOLU_BC2)
      call set_tr_wd_type(n, npart)
      end function M_BC2_BC_setSpec

      integer function N_BC2_1_setSpec(n) result(n_N_BC2_1)
      integer, intent(inout) :: n
      n = n + 1
      n_N_BC2_1 = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 1.d+0)
      call set_trpdens(n, DENS_BCAR)
      call set_trradius(n, DG_BC2 * .5d-6)
      call set_fq_aer(n, SOLU_BC2)
      call set_tr_wd_type(n, npart)
      end function N_BC2_1_setSpec

      integer function M_BC3_SU_setSpec(n) result(n_M_BC3_SU)
      integer, intent(inout) :: n
      n = n + 1
      n_M_BC3_SU = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 96.d0)
      call set_trpdens(n, DENS_SULF)
      call set_trradius(n, DG_BC3 * .5d-6)
      call set_fq_aer(n, SOLU_BC3)
      call set_tr_wd_type(n, npart)
      end function M_BC3_SU_setSpec

      integer function M_BC3_BC_setSpec(n) result(n_M_BC3_BC)
      integer, intent(inout) :: n
      n = n + 1
      n_M_BC3_BC = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 12.d0)
      call set_trpdens(n, DENS_BCAR)
      call set_trradius(n, DG_BC3 * .5d-6)
      call set_fq_aer(n, SOLU_BC3)
      call set_tr_wd_type(n, npart)
      end function M_BC3_BC_setSpec

      integer function N_BC3_1_setSpec(n) result(n_N_BC3_1)
      integer, intent(inout) :: n
      n = n + 1
      n_N_BC3_1 = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 1.d+0)
      call set_trpdens(n, DENS_BCAR)
      call set_trradius(n, DG_BC3 * .5d-6)
      call set_fq_aer(n, SOLU_BC3)
      call set_tr_wd_type(n, npart)
      end function N_BC3_1_setSpec

      integer function M_DBC_SU_setSpec(n) result(n_M_DBC_SU)
      integer, intent(inout) :: n
      n = n + 1
      n_M_DBC_SU = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 96.d0)
      call set_trpdens(n, DENS_SULF)
      call set_trradius(n, DG_DBC * .5d-6)
      call set_fq_aer(n, SOLU_DBC)
      call set_tr_wd_type(n, npart)
      end function M_DBC_SU_setSpec

      integer function M_DBC_BC_setSpec(n) result(n_M_DBC_BC)
      integer, intent(inout) :: n
      n = n + 1
      n_M_DBC_BC = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 96.d0)
      call set_trpdens(n, DENS_BCAR)
      call set_trradius(n, DG_DBC * .5d-6)
      call set_fq_aer(n,  SOLU_DBC)
      call set_tr_wd_type(n, npart)
      end function M_DBC_BC_setSpec

      integer function M_DBC_DU_setSpec(n) result(n_M_DBC_DU)
      integer, intent(inout) :: n
      n = n + 1
      n_M_DBC_DU = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 1.d+0)
      call set_trpdens(n, DENS_DUST)
      call set_trradius(n, DG_DBC * .5d-6)
      call set_fq_aer(n,  SOLU_DBC)
      call set_tr_wd_type(n, npart)
      end function M_DBC_DU_setSpec

      integer function N_DBC_1_setSpec(n) result(n_N_DBC_1)
      integer, intent(inout) :: n
      n = n + 1
      n_N_DBC_1 = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 1.d+0)
      call set_trpdens(n, DENS_DUST)
      call set_trradius(n, DG_DBC * .5d-6)
      call set_fq_aer(n,  SOLU_DBC)
      call set_tr_wd_type(n, npart)
      end function N_DBC_1_setSpec

      integer function M_BOC_SU_setSpec(n) result(n_M_BOC_SU)
      integer, intent(inout) :: n
      n = n + 1
      n_M_BOC_SU = n
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 96.d0)
      call set_trpdens(n, DENS_SULF)
      call set_trradius(n, DG_BOC * .5d-6)
      call set_fq_aer(n, SOLU_BOC)
      call set_tr_wd_type(n, npart)
      end function M_BOC_SU_setSpec

      integer function M_BOC_BC_setSpec(n) result(n_M_BOC_BC)
      integer, intent(inout) :: n
      n = n + 1
      n_M_BOC_BC = n
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 12.d0)
      call set_trpdens(n, DENS_BCAR)
      call set_trradius(n, DG_BOC * .5d-6)
      call set_fq_aer(n,  SOLU_BOC)
      call set_tr_wd_type(n, npart)
      end function M_BOC_BC_setSpec

      integer function M_BOC_OC_setSpec(n) result(n_M_BOC_OC)
      integer, intent(inout) :: n
      n = n + 1
      n_M_BOC_OC = n
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 15.6d0)
      call set_trpdens(n, DENS_OCAR)
      call set_trradius(n, DG_BOC * .5d-6)
      call set_fq_aer(n,  SOLU_BOC)
      call set_tr_wd_type(n, npart)
      end function M_BOC_OC_setSpec

      integer function N_BOC_1_setSpec(n) result(n_N_BOC_1)
      integer, intent(inout) :: n
      n = n + 1
      n_N_BOC_1 = n
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 1.d+0)
      call set_trpdens(n, DENS_BCAR)
      call set_trradius(n, DG_BOC * .5d-6)
      call set_fq_aer(n,  SOLU_BOC)
      call set_tr_wd_type(n, npart)
      end function N_BOC_1_setSpec

      integer function M_BCS_SU_setSpec(n) result(n_M_BCS_SU)
      integer, intent(inout) :: n
      n = n + 1
      n_M_BCS_SU = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 96.d0)
      call set_trpdens(n, DENS_SULF)
      call set_trradius(n, DG_BCS * .5d-6)
      call set_fq_aer(n, SOLU_BCS)
      call set_tr_wd_type(n, npart)
      end function M_BCS_SU_setSpec

      integer function M_BCS_BC_setSpec(n) result(n_M_BCS_BC)
      integer, intent(inout) :: n
      n = n + 1
      n_M_BCS_BC = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 12.d0)
      call set_trpdens(n, DENS_BCAR)
      call set_trradius(n, DG_BCS * .5d-6)
      call set_fq_aer(n, SOLU_BCS)
      call set_tr_wd_type(n, npart)
      end function M_BCS_BC_setSpec

      integer function N_BCS_1_setSpec(n) result(n_N_BCS_1)
      integer, intent(inout) :: n
      n = n + 1
      n_N_BCS_1 = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 1.d+0)
      call set_trpdens(n, DENS_BCAR)
      call set_trradius(n, DG_BCS * .5d-6)
      call set_fq_aer(n, SOLU_BCS)
      call set_tr_wd_type(n, npart)
      end function N_BCS_1_setSpec

      integer function M_MXX_SU_setSpec(n) result(n_M_MXX_SU)
      integer, intent(inout) :: n
      n = n + 1
      n_M_MXX_SU = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 96.d0)
      call set_trpdens(n, DENS_SULF)
      call set_trradius(n, DG_MXX * .5d-6)
      call set_fq_aer(n, SOLU_MXX)
      call set_tr_wd_type(n, npart)
      end function M_MXX_SU_setSpec

      integer function M_MXX_BC_setSpec(n) result(n_M_MXX_BC)
      integer, intent(inout) :: n
      n = n + 1
      n_M_MXX_BC = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 12.d0)
      call set_trpdens(n, DENS_BCAR)
      call set_trradius(n, DG_MXX * .5d-6)
      call set_fq_aer(n, SOLU_MXX)
      call set_tr_wd_type(n, npart)
      end function M_MXX_BC_setSpec

      integer function M_MXX_OC_setSpec(n) result(n_M_MXX_OC)
      integer, intent(inout) :: n
      n = n + 1
      n_M_MXX_OC = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 15.6d0)
      call set_trpdens(n, DENS_OCAR)
      call set_trradius(n, DG_MXX * .5d-6)
      call set_fq_aer(n, SOLU_MXX)
      call set_tr_wd_type(n, npart)
      end function M_MXX_OC_setSpec

      integer function M_MXX_DU_setSpec(n) result(n_M_MXX_DU)
      integer, intent(inout) :: n
      n = n + 1
      n_M_MXX_DU = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 1.d+0)
      call set_trpdens(n, DENS_DUST)
      call set_trradius(n, DG_MXX * .5d-6)
      call set_fq_aer(n, SOLU_MXX)
      call set_tr_wd_type(n, npart)
      end function M_MXX_DU_setSpec

      integer function M_MXX_SS_setSpec(n) result(n_M_MXX_SS)
      integer, intent(inout) :: n
      n = n + 1
      n_M_MXX_SS = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 75.d0)
      call set_trpdens(n, DENS_SEAS)
      call set_trradius(n, DG_MXX * .5d-6)
      call set_fq_aer(n, SOLU_MXX)
      call set_tr_wd_type(n, npart)
      end function M_MXX_SS_setSpec

      integer function N_MXX_1_setSpec(n) result(n_N_MXX_1)
      integer, intent(inout) :: n
      n = n + 1
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
      end function N_MXX_1_setSpec

      integer function M_OCS_SU_setSpec(n) result(n_M_OCS_SU)
      integer, intent(inout) :: n
      n = n + 1
      n_M_OCS_SU = n
          call set_ntm_power(n, -11)
          call set_ntsurfsrc(n,  0)
          call set_tr_mm(n, 96.d0)
          call set_trpdens(n, DENS_SULF)
          call set_trradius(n, DG_OCS * .5d-6)
          call set_fq_aer(n, SOLU_OCS)
          call set_tr_wd_type(n, npart)
      end function M_OCS_SU_setSpec

      integer function M_OCS_OC_setSpec(n) result(n_M_OCS_OC)
      integer, intent(inout) :: n
      n = n + 1
      n_M_OCS_OC = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 15.6d0)
      call set_trpdens(n, DENS_OCAR)
      call set_trradius(n, DG_OCS * .5d-6)
      call set_fq_aer(n, SOLU_OCS)
      call set_tr_wd_type(n, npart)
      end function M_OCS_OC_setSpec

      integer function N_OCS_1_setSpec(n) result(n_N_OCS_1)
      integer, intent(inout) :: n
      n = n + 1
      n_N_OCS_1 = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 1.d+0)
      call set_trpdens(n, DENS_OCAR)
      call set_trradius(n, DG_OCS * .5d-6)
      call set_fq_aer(n, SOLU_OCS)
      call set_tr_wd_type(n, npart)
      end function N_OCS_1_setSpec

#endif /* TRACERS_AMP */

#ifdef TRACERS_AEROSOLS_VBS
      integer function VBS_setSpec(name, n, index, type) result(label)
      character(len=*), intent(in) :: name
      integer, intent(inout) :: n
      integer, intent(in) :: index
      character(len=4), intent(in) :: type

      n = n + 1
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

#ifdef TRACERS_TOMAS
      integer function TOMAS_H2SO4_setSpec(n) result(n_H2SO4)
      integer, intent(inout) :: n
      n = n + 1
      n_H2SO4 = n
      call set_ntm_power(n, -11)
      call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 98.d0)
      call set_trpdens(n, 1.78d0)
      call set_fq_aer(n, 1.d0)
      call set_tr_wd_type(n, nGas)
      end function TOMAS_H2SO4_setSpec

      integer function TOMAS_SOAgas_setSpec(n) result(n_SOAgas)
      integer, intent(inout) :: n
      n = n + 1
      n_SOAgas = n
      call set_ntm_power(n, -11)
!     call set_ntsurfsrc(n,  0)
      call set_tr_mm(n, 120.10d0) ! i.e. 10 carbons
#ifdef TRACERS_DRYDEP
      call set_HSTAR(n,  0.D0)  !no dry dep
#endif
      end function TOMAS_SOAgas_setSpec

      function TOMAS_setSpec(func, n, nbins) result (indices)
      interface
        integer function func(n, bin)
        integer, intent(inout) :: n
        integer, intent(in) :: bin
        end function func
      end interface
      integer, intent(inout) :: n
      integer, intent(in) :: nbins
      integer :: indices(nbins)

      integer :: bin

      do bin = 1, nbins
        indices(bin) = func(n, bin)
      end do
      
      end function TOMAS_setSpec

      integer function TOMAS_ANUM_setSpec(n, bin) result(n_ANUM)
      integer, intent(inout) :: n
      integer, intent(in) :: bin
      n = n + 1
      n_ANUM = n  
      
      TOMAS_dens=1.5d3
      TOMAS_radius=(sqrt(xk(bin)*xk(bin+1))/TOMAS_dens/pi/4.*3.)
     *     **(1./3.)  
      
      if(bin.le.5)call set_ntm_power(n, 10)
      if(bin.gt.5)call set_ntm_power(n, 8)
      if(bin.eq.1)  call set_ntsurfsrc(n,  3) ! 1 for SO4,2 for EC, 3 for OC (4 for SS and 5 for DU)
      call set_tr_mm(n, 1.d+0)
      call set_trpdens(n, TOMAS_dens)
      call set_trradius(n, TOMAS_radius)
      call set_fq_aer(n, 1.d0)  !not used in wet deposition
      call set_tr_wd_type(n, npart)
      end function TOMAS_ANUM_setSpec

      integer function TOMAS_ASO4_setSpec(n,bin) result(n_ASO4)
      integer, intent(inout) :: n
      integer, intent(in) :: bin
      n = n + 1
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

      integer function TOMAS_ANACL_setSpec(n,bin) result(n_ANACL)
      integer, intent(inout) :: n
      integer, intent(in) :: bin
      n = n + 1
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

      integer function TOMAS_AECOB_setSpec(n,bin) result(n_AECOB)
      integer, intent(inout) :: n
      integer, intent(in) :: bin
      n = n + 1
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

      integer function TOMAS_AECIL_setSpec(n,bin) result(n_AECIL)
      integer, intent(inout) :: n
      integer, intent(in) :: bin
      n = n + 1
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

      integer function TOMAS_AOCOB_setSpec(n,bin) result(n_AOCOB)
      integer, intent(inout) :: n
      integer, intent(in) :: bin
      n = n + 1
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

      if (bin==1) call sync_param("OCB_om2oc",om2oc(n_AOCOB))
      end function TOMAS_AOCOB_setSpec

      integer function TOMAS_AOCIL_setSpec(n,bin) result(n_AOCIL)
      integer, intent(inout) :: n
      integer, intent(in) :: bin
      n = n + 1
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

      integer function TOMAS_ADUST_setSpec(n,bin) result(n_ADUST)
      integer, intent(inout) :: n
      integer, intent(in) :: bin
      n = n + 1
      n_ADUST = n  
      call set_ntsurfsrc(n,  0)
      call set_ntisurfsrc(n, 1)
      if(bin.le.10) TOMAS_dens= 2.5d3 !clay 
      if(bin.gt.10) TOMAS_dens= 2.65d3 !silt
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

      integer function TOMAS_AH2O_setSpec(n,bin) result(n_AH2O)
      integer, intent(inout) :: n
      integer, intent(in) :: bin
      n = n + 1
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
