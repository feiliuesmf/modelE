      SUBROUTINE initTracerGriddedData()
!@sum init_tracer initializes trace gas attributes and diagnostics
!@auth J. Lerner
!@calls sync_param, SET_TCON, RDLAND, RDDRYCF
      USE DOMAIN_DECOMP_ATM, only:GRID,GET,AM_I_ROOT,
     &     write_parallel,readt8_parallel
      USE CONSTANT, only: mair,mwat,sday,pi
#ifdef TRACERS_AEROSOLS_SOA
     &                   ,gasc
#endif  /* TRACERS_AEROSOLS_SOA */
      USE RESOLUTION, only : jm,lm
      USE MODEL_COM, only: dtsrc,itime
      USE ATM_COM, only: pmidl00
      USE GEOM, only: axyp,byaxyp
      USE ATM_COM, only: am  ! Air mass of each box (kg/m^2)
      USE TRACER_COM
#ifdef TRACERS_ON
      USE TRDIAG_COM
#endif
      USE Dictionary_mod
#ifdef TRACERS_SPECIAL_Lerner
      USE TRACERS_MPchem_COM, only: n_MPtable,tcscale
!@dbparam dsol describes portion of solar cycle being modeled for linoz
!@+      +1.0 = solar max, 0.0 = neutral, -1.0 = solar min
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
     &   ,fracSiltPDFscheme
#endif
#ifdef TRACERS_QUARZHEM
     &   ,DenHema,DenQuarz,FreeFe,FrHeQu
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
     &gw_llim,sw_lim,exclude_us_eu,nn_or_zon,ice_age,nday_ch4,max_days,
     &ns_wet,nra_ch4
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
     & DG_DS1, DG_DS2, DG_SSA, DG_SSC, DG_ACC,
     & DG_SSS, DG_OCC, DG_BC1, DG_BC2, DG_BC3,
     & DG_DBC, DG_BOC, DG_BCS, DG_OCS, DG_MXX,

     & SOLU_DD1, SOLU_DD2, SOLU_AKK, SOLU_ACC,
     & SOLU_DS1, SOLU_DS2, SOLU_SSA, SOLU_SSC,
     & SOLU_SSS, SOLU_OCC, SOLU_BC1, SOLU_BC2, SOLU_BC3,
     & SOLU_DBC, SOLU_BOC, SOLU_BCS, SOLU_OCS, SOLU_MXX

      USE AERO_ACTV, only: DENS_SULF, DENS_DUST,
     &          DENS_SEAS, DENS_BCAR, DENS_OCAR
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

#if (defined TRACERS_OCEAN) && !defined(TRACERS_OCEAN_INDEP)
! atmosphere copies atmosphere-declared tracer info to ocean
! so that the ocean can "inherit" it without referencing atm. code
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
      implicit none
      integer :: l,k,n,kr,m,ns
#ifdef TRACERS_SPECIAL_O18
      real*8 fracls
#endif
#ifdef TRACERS_TOMAS
      integer :: bin
      real*8 :: TOMAS_dens,TOMAS_radius
#endif
#if (defined TRACERS_WATER) || (defined TRACERS_DRYDEP)
!@param convert_HSTAR converts from mole/Joule to mole/(L*atm)
      real*8, parameter :: convert_HSTAR = 1.01325d2
#endif
#ifdef TRACERS_SPECIAL_Shindell
!@var PRES local nominal pressure for vertical interpolations
!@var iu_data unit number
!@var title header read in from file
      REAL*8, DIMENSION(LM) :: PRES
      integer iu_data,i,j,nq
      character*80 title
      character(len=300) :: out_line
      real*8, dimension(6) :: temp_ghg
      integer :: temp_year
#endif /* TRACERS_SPECIAL_Shindell */

#if defined(TRACERS_GASEXCH_ocean) && defined(TRACERS_GASEXCH_ocean_CFC)
      integer i, iu_data
#endif
#if (!defined(TRACERS_GASEXCH_ocean_CO2)) && defined(TRACERS_GASEXCH_land_CO2)
      real*8 :: atmCO2 = 280.d0
#endif

! temp storage for new tracer interfaces
      integer :: values(ntm)
      integer :: val

      INTEGER J_0, J_1, I_0, I_1

#ifdef TRACERS_TOMAS
      CALL initbounds()
      call readbinact ("binact10_12.dat",binact10) 
      call readbinact ("binact02_12.dat",binact02) 
      call readfraction("fraction10_12.dat",fraction10) 
      call readfraction("fraction02_12.dat",fraction02) 
      call readmielut ! aerosol radiation lookup table
#endif
C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT=J_0,       J_STOP=J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

      do n=1,ntm

      select case (trname(n))

#ifdef TRACERS_ON
      case ('N2O')
#ifdef TRACERS_SPECIAL_Shindell
          call openunit('N2O_IC',iu_data,.true.,.true.)
          CALL READT8_PARALLEL(grid,iu_data,NAMEUNIT(iu_data),N2OICIN,0)
          call closeunit(iu_data)
          do j=J_0,J_1  ; do i=I_0,I_1
           N2OICINL(:)=N2OICIN(i,j,:) ! now in PPPM
           CALL LOGPINT(LCOalt,PCOalt,N2OICINL,LM,PRES,N2OICL,.true.)
           N2OICX(I,J,:)=N2OICL(:)*am(:,i,j)*axyp(i,j)
          end do     ; end do
#endif
      case ('CH4')
#ifdef TRACERS_SPECIAL_Shindell
C**** determine initial CH4 distribution if set from rundeck
C**** This is only effective with a complete restart.
          call sync_param("ch4_init_sh",ch4_init_sh)
          call sync_param("ch4_init_nh",ch4_init_nh)
          call sync_param("fix_CH4_chemistry",fix_CH4_chemistry)
          call sync_param("scale_ch4_IC_file",scale_ch4_IC_file)
C         Interpolate CH4 altitude-dependence to model resolution:
          CALL LOGPINT(LCH4alt,PCH4alt,CH4altINT,LM,PRES,CH4altT,.true.)
          CALL LOGPINT(LCH4alt,PCH4alt,CH4altINX,LM,PRES,CH4altX,.true.)
          if(fix_CH4_chemistry.eq.-1)then
            call openunit('CH4_IC',iu_data,.true.,.true.)
            CALL READT8_PARALLEL(grid,iu_data,NAMEUNIT(iu_data),
     &           CH4ICIN,0)
            call closeunit(iu_data)
            do j=J_0,J_1  ; do i=I_0,I_1
             CH4ICINL(:)=CH4ICIN(I,J,:)! now in PPPM
             CALL LOGPINT(LCOalt,PCOalt,CH4ICINL,LM,PRES,CH4ICL,.true.)
             CH4ICX(I,J,:)=CH4ICL(:)*scale_ch4_IC_file*am(:,i,j)*
     *            axyp(i,j)
            end do     ; end do
          end if
#endif /* TRACERS_SPECIAL_Shindell */
#ifdef TRACERS_SPECIAL_Shindell
      case ('Ox')
          call openunit('Ox_IC',iu_data,.true.,.true.)
          CALL READT8_PARALLEL(grid,iu_data,NAMEUNIT(iu_data),OxICIN,0)
          call closeunit(iu_data)
          do j=J_0,J_1  ; do i=I_0,I_1
           OxICINL(:)=OxICIN(I,J,:)! now in PPPM
           CALL LOGPINT(LCOalt,PCOalt,OxICINL,LM,PRES,OxICL,.true.)
           OxIC(I,J,:)=OxICL(:)*am(:,i,j)*axyp(i,j)
          end do     ; end do
          call check_aircraft_sectors ! special 3D source case

      case ('CFC')
          if(AM_I_ROOT( ))then
C          check on GHG file's 1995 value for CFCs:
           call openunit('GHG',iu_data,.false.,.true.)
           do i=1,5; read(iu_data,'(a80)') title; enddo
           temp_year=0
           do while(temp_year <= 1995)
             read(iu_data,*,end=101) temp_year,(temp_ghg(j),j=1,6)
             if(temp_year==1995)then
               temp_ghg(1)=cfc_rad95*0.95d0
               temp_ghg(2)=cfc_rad95*1.05d0
               temp_ghg(3)=(temp_ghg(4)+temp_ghg(5))*1.d-9
               if(temp_ghg(3) < temp_ghg(1) .or.
     &         temp_ghg(3) > temp_ghg(2))then
                 call stop_model('please check on cfc_rad95 2',255)
               endif
             endif
           enddo
 101       continue
           if(temp_year<1995)
     &     call stop_model('please check on cfc_rad95 1',255)
           call closeunit(iu_data)
          endif
C          read the CFC initial conditions:
          call openunit('CFC_IC',iu_data,.true.,.true.)
          CALL READT8_PARALLEL(grid,iu_data,NAMEUNIT(iu_data),CFCICIN,0)
          call closeunit(iu_data)
          do j=J_0,J_1  ; do i=I_0,I_1
           CFCICINL(:)=CFCICIN(I,J,:)! now in PPPM
           CALL LOGPINT(LCOalt,PCOalt,CFCICINL,LM,PRES,CFCICL,.true.)
           CFCIC(I,J,:)=CFCICL(:)*am(:,i,j)*axyp(i,j)
          end do     ; end do

      case ('CO')
          call openunit('CO_IC',iu_data,.true.,.true.)
          CALL READT8_PARALLEL(grid,iu_data,NAMEUNIT(iu_data),COICIN,0)
          call closeunit(iu_data)
          do j=J_0,J_1  ; do i=I_0,I_1
           COICINL(:)=COICIN(I,J,:)! now in PPPM
           CALL LOGPINT(LCOalt,PCOalt,COICINL,LM,PRES,COICL,.true.)
           COIC(I,J,:)=COICL(:)*am(:,i,j)*axyp(i,j)
          end do     ; end do
#endif  /* TRACERS_SPECIAL_Shindell */

#endif /* TRACERS_ON */
      end select
      end do

      return
      end subroutine initTracerGriddedData
