#include "rundeck_opts.h"
!@sum  TRACERS_DRV: tracer-dependent routines for air/water mass
!@+    and ocean tracers
!@+    Routines included:
!@+      Those that MUST EXIST for all tracers:
!@+        Diagnostic specs: init_tracer
!@+        Tracer initialisation + sources: tracer_ic, set_tracer_source
!@+        Entry points: daily_tracer
!@auth Jean Lerner/Gavin Schmidt
!@ver  1.0

      SUBROUTINE init_tracer
!@sum init_tracer initializes trace gas attributes and diagnostics
!@auth J. Lerner
!@calls sync_param, SET_TCON, RDLAND, RDDRYCF
      USE DOMAIN_DECOMP, only : GRID, GET
      USE CONSTANT, only: mair,mwat,sday
      USE MODEL_COM, only: dtsrc,byim,ptop,psf,sig,lm,jm,itime
      USE DIAG_COM, only: ia_src,ia_12hr,ir_log2,npts,ia_rad
      USE TRACER_COM
#ifdef TRACERS_ON
      USE TRDIAG_COM
#endif
      USE PARAM
#ifdef TRACERS_SPECIAL_Lerner
      USE TRACERS_MPchem_COM, only: n_MPtable,tcscale
!@dbparam dsol describes portion of solar cycle being modeled for linoz
!@+      +1.0 = solar max, 0.0 = neutral, -1.0 = solar min
      USE LINOZ_CHEM_COM, only: dsol
#endif
#ifdef TRACERS_WATER
      USE LANDICE_COM, only : trli0    ! should these be in tracer_com?
      USE SEAICE_COM, only : trsi0
#endif
#ifdef TRACERS_SPECIAL_Shindell
      USE TRCHEM_Shindell_COM,only:COaltIN,LCOalt,PCOalt,COalt,
     & mass2vol,bymass2vol,CH4altINT,CH4altINX,LCH4alt,PCH4alt,
     &     CH4altX,CH4altT,ch4_init_sh,ch4_init_nh,
     &     OxICIN,OxIC,OxICINL,OxICL,corrOxIN,corrOx,LcorrOx,PcorrOx
     &     ,pfix_CH4_N,pfix_CH4_S,fix_CH4_chemistry,which_trop,
     &     PI_run,PIratio_N,PIratio_CO_T,PIratio_CO_S,PIratio_other,
     &     PIratio_indus,PIratio_bburn
#ifdef SHINDELL_STRAT_CHEM
     &     ,BrOxaltIN,ClOxaltIN,ClONO2altIN,HClaltIN,BrOxalt,
     &     ClOxalt,ClONO2alt,HClalt,N2OICIN,N2OICX,N2OICINL,N2OICL,
     &     CFCICIN,CFCIC,CFCICINL,CFCICL
#endif
#endif
#ifdef TRACERS_COSMO
      USE COSMO_SOURCES, only: be7_src_param
#endif
      USE FILEMANAGER, only: openunit,closeunit
      implicit none
      integer :: l,k,n,ntemp,n2,ltop,g
      character*20 sum_unit(ntm),inst_unit(ntm)   ! for conservation
      character*10 CMR
#ifdef TRACERS_ON
      logical :: qcon(KTCON-1), qsum(KTCON-1), T=.TRUE. , F=.FALSE.
#endif
      character*50 :: unit_string
#ifdef TRACERS_WATER
      real*8 fracls
#endif
#if (defined TRACERS_WATER) || (defined TRACERS_DRYDEP)
!@param convert_HSTAR converts from mole/Joule to mole/(L*atm)
      real*8, parameter :: convert_HSTAR = 1.01325d2
#endif
#ifdef TRACERS_SPECIAL_Shindell
!@var PRES local nominal pressure for vertical interpolations
!@var iu_data unit number
!@var title header read in from file
      REAL*8, DIMENSION(LM) :: PRES, tempOx2
      REAL*8, DIMENSION(LcorrOx) :: tempOx1
      integer iu_data,m,i,j,nq
      character*80 title
#ifdef regional_Ox_tracers
!@var Ox_a_tracer logical is true if Ox is one of the tracers
      logical Ox_a_tracer
#endif
#endif
      INTEGER J_0, J_1
C****                
C**** Extract useful local domain parameters from "grid"
C****                   
      CALL GET(grid, J_STRT=J_0,       J_STOP=J_1)

C**** Set defaults for tracer attributes (all dimensioned ntm)
      itime_tr0=itime
      t_qlimit = .true.
      trdecay = 0.
#ifdef TRACERS_ON
      needtrs = .false.
      ntsurfsrc = 0
      trpdens = 0.
      trradius = 0.
#ifdef TRACERS_SPECIAL_Lerner
      n_MPtable = 0
      tcscale = 0.
#endif
#endif
      tr_wd_TYPE = nGas       !other options are nPART or nWATER
      tr_RKD = 0.
      tr_DHD = 0.
      fq_aer = 0.
#ifdef TRACERS_WATER
      trli0 = 0.
      trsi0 = 0.
      tr_H2ObyCH4 = 0.
      dowetdep = .false.
#endif
#ifdef TRACERS_DRYDEP
      dodrydep = .false.
      F0 = 0.
      HSTAR = 0.   ! tr_RKD * convert_HSTAR
#endif
#if (defined TRACERS_WATER) || (defined TRACERS_OCEAN)
      trw0 = 0.
      ntrocn = 0
#endif
#ifdef TRACERS_OCEAN
      trglac = 0.
#endif
#ifdef TRACERS_SPECIAL_Shindell
      call sync_param("which_trop",which_trop)
      call sync_param("PI_run",PI_run)
      call sync_param("PIratio_N",PIratio_N)
      call sync_param("PIratio_CO_T",PIratio_CO_T)
      call sync_param("PIratio_CO_S",PIratio_CO_S)
      call sync_param("PIratio_other",PIratio_other)
      call sync_param("PIratio_indus",PIratio_indus)
      call sync_param("PIratio_bburn",PIratio_bburn)
      PRES(:)=SIG(:)*(PSF-PTOP)+PTOP
#ifdef regional_Ox_tracers
      Ox_a_tracer=.false.
      do n2=1,ntm; if(trname(n2).eq.'Ox')Ox_a_tracer=.true.; end do
      if(.not.Ox_a_tracer)call stop_model(
     &'Ox must be tracer if using regional Ox tracers; define 1st',255)
      regOx_t=0.d0; regOx_b=0.d0; regOx_n=0.d0
      regOx_s=0.d0; regOx_e=0.d0; regOx_w=0.d0
      ntemp=0
#endif
C**** initialise source arrays
       do l=1,lm
       do j=1,jm
       do i=1,im
       oh_live(i,j,l)=0.0d0
       no3_live(i,j,l)=0.0d0
       end do
       end do
       end do
#endif

C**** Define a max layer for some optionally trop/strat tracers
      LTOP = LM

#ifdef TRACERS_COSMO
C**** get rundeck parameter for cosmogenic source factor
      call sync_param("be7_src_param", be7_src_param)
#endif

C**** Define individual tracer characteristics
      do n=1,ntm
      select case (trname(n))

#ifdef TRACERS_ON
      case ('Air')
      n_Air = n
          ntm_power(n) = -2
          tr_mm(n) = mair

      case ('SF6')
      n_SF6 = n
          ntm_power(n) = -14
          tr_mm(n) = 146.01d0
          ntsurfsrc(n) = 1

      case ('Rn222')
      n_Rn222 = n
          ntm_power(n) = -21
          tr_mm(n) = 222.d0
          trdecay(n) = 2.1d-6
          ntsurfsrc(n) = 1

      case ('CO2')
      n_CO2 = n
          ntm_power(n) = -6
          tr_mm(n) = 44.d0
          t_qlimit(n) = .false.
          ntsurfsrc(n) = 6

      case ('N2O')
      n_N2O = n
          ntm_power(n) = -9
          tr_mm(n) = 44.d0
#ifdef TRACERS_SPECIAL_Lerner
          ntsurfsrc(n) = 1
          n_MPtable(n) = 1
          tcscale(n_MPtable(n)) = 1.
#endif
#ifdef SHINDELL_STRAT_CHEM
          call openunit('N2O_IC',iu_data,.true.,.true.)
          read (iu_data) title,N2OICIN
          call closeunit(iu_data)
          write(6,*) title,' read from N2O_IC'
          do j=J_0,J_1  ; do i=1,im
           N2OICINL(:)=N2OICIN(I,J,:)
           CALL LOGPINT(LCOalt,PCOalt,N2OICINL,LM,PRES,N2OICL,.true.)
           N2OICX(I,J,:)=N2OICL(:)
          end do     ; end do
#endif

      case ('CFC11')   !!! should start April 1
      n_CFC11 = n
          ntm_power(n) = -12
          tr_mm(n) = 137.4d0
          ntsurfsrc(n) = 1
#ifdef TRACERS_SPECIAL_Lerner
          n_MPtable(n) = 2
          tcscale(n_MPtable(n)) = 1.
#endif

      case ('14CO2')   !!! should start 10/16
      n_14CO2 = n
          ntm_power(n) = -18
          tr_mm(n) = 46.d0
          ntsurfsrc(n) = 1

      case ('CH4')
      n_CH4 = n
          tr_mm(n) = 16.d0
#ifdef TRACERS_SPECIAL_Lerner
          ntsurfsrc(n) = 14
          ntm_power(n) = -9
          n_MPtable(n) = 3
          tcscale(n_MPtable(n)) = 1.
#endif
#ifdef TRACERS_SPECIAL_Shindell
C**** determine initial CH4 distribution if set from rundeck
C**** This is only effective with a complete restart.
          call sync_param("ch4_init_sh",ch4_init_sh)
          call sync_param("ch4_init_nh",ch4_init_nh)
          call sync_param("pfix_CH4_N",pfix_CH4_N)
          call sync_param("pfix_CH4_S",pfix_CH4_S)
          call sync_param("fix_CH4_chemistry",fix_CH4_chemistry)
#ifdef EDGAR_HYDE_SOURCES
          ntsurfsrc(n) = 15
#else
          ntsurfsrc(n) = 14
#endif
          ntm_power(n) = -8
C         Interpolate CH4 altitude-dependence to model resolution:
          CALL LOGPINT(LCH4alt,PCH4alt,CH4altINT,LM,PRES,CH4altT,.true.)
          CALL LOGPINT(LCH4alt,PCH4alt,CH4altINX,LM,PRES,CH4altX,.true.)
#endif

      case ('O3')
      n_O3 = n
          ntm_power(n) = -8
          tr_mm(n) = 48.d0
#ifdef TRACERS_SPECIAL_Lerner
C**** Get solar variability coefficient from namelist if it exits
          dsol = 0.
          call sync_param("dsol",dsol)
#endif
#endif
      case ('Water')
      n_Water = n
#if (defined TRACERS_WATER) || (defined TRACERS_OCEAN)
          trw0(n) = 1.
          ntrocn(n)= 0
#endif
#ifdef TRACERS_WATER
          ntm_power(n) = -4
          tr_mm(n) = mwat
          needtrs(n) = .true.
          tr_wd_TYPE(n) = nWater
          trli0(n) = 1.
          trsi0(n) = 1.
          tr_H2ObyCH4(n) = 1.
#endif
#ifdef TRACERS_OCEAN
          trglac(n) = 1.
#endif

#ifdef TRACERS_ON
#ifdef TRACERS_SPECIAL_O18
      case ('H2O18')
      n_H2O18 = n
          ntm_power(n) = -7
          tr_mm(n) = 20.
          needtrs(n) = .true.
          tr_wd_TYPE(n) = nWater
          trw0(n) = 2.228d-3   ! SMOW mass ratio of water molecules
          trli0(n) = 0.980d0*trw0(n)  ! d=-20
          trsi0(n) = fracls(trname(n))*trw0(n)
          tr_H2ObyCH4(n) = trw0(n)*1.023d0 ! d=+23 (ie. no frac from O2)
          ntrocn(n) = -3
#ifdef TRACERS_OCEAN
          trglac(n) = trw0(n)*0.98d0   ! d=-20
#endif

      case ('HDO')
      n_HDO = n
          ntm_power(n) = -8
          tr_mm(n) = 19.
          needtrs(n) = .true.
          tr_wd_TYPE(n) = nWater
          trw0(n) = 3.29d-4    ! SMOW mass ratio of water molecules
          trli0(n) = 0.830d0*trw0(n)  ! d=-170
          trsi0(n) = fracls(trname(n))*trw0(n)
          tr_H2ObyCH4(n) = trw0(n)*0.93d0  ! d=-70
          ntrocn(n) = -4
#ifdef TRACERS_OCEAN
          trglac(n) = trw0(n)*0.84d0   ! d=-160
#endif

      case ('HTO')
      n_HTO = n
          ntm_power(n) = -18
          tr_mm(n) = 20.
          needtrs(n) = .true.
          tr_wd_TYPE(n) = nWater
          trw0(n) = 0. !2.22d-18   ! SMOW mass ratio of water molecules
          trli0(n) = 0.
          trsi0(n) = 0.
          tr_H2ObyCH4(n) = 0.
          trdecay(n) = 1.77d-9      ! =5.59d-2 /yr
          ntrocn(n) = -18
#ifdef TRACERS_OCEAN
          trglac(n) = 0.
#endif
#endif
#ifdef TRACERS_SPECIAL_Shindell
      case ('Ox')
      n_Ox = n
          call openunit('Ox_IC',iu_data,.true.,.true.)
          read (iu_data) title,OxICIN
          call closeunit(iu_data)
          write(6,*) title,' read from OxIC'
          do j=J_0,J_1  ; do i=1,im
           OxICINL(:)=OxICIN(I,J,:)
           CALL LOGPINT(LCOalt,PCOalt,OxICINL,LM,PRES,OxICL,.true.)
           OxIC(I,J,:)=OxICL(:)
          end do     ; end do
c         read stratospheric correction from files:
          call openunit('Ox_corr',iu_data,.true.,.true.)
          read (iu_data) title,corrOxIN
          call closeunit(iu_data)
          write(6,*) title,' read from Ox_corr'
          DO m=1,12; DO j=1,jm
           tempOx1(:)=CorrOxIN(J,:,M)
           CALL LOGPINT(LcorrOX,PcorrOx,tempOx1,LM,PRES,
     &                  tempOx2,.true.)
           CorrOx(J,:,M)=tempOx2(:)
          END DO   ; END DO
C         Only alter Ox between 250 and 30 hPa:
          DO L=1,LM
            IF(PRES(L).lt.30.d0.or.PRES(L).gt.250.d0)
     &      corrOx(:,L,:)=1.0d0
          END DO
          ntm_power(n) = -8
          tr_mm(n) = 48.d0
#ifdef TRACERS_DRYDEP
          F0(n) = 1.8d0
          HSTAR(n) = 1.d-2
#endif

      case ('NOx')
      n_NOx = n
          ntm_power(n) = -11
#ifdef EDGAR_HYDE_SOURCES
          ntsurfsrc(n) = 7
#else
          ntsurfsrc(n) = 3 ! fossil fuel, biomass burning, soil
#endif
          tr_mm(n) = 14.01d0
#ifdef TRACERS_DRYDEP
          F0(n) = 1.d-1
          HSTAR(n) = 1.d-2
#endif

      case ('N2O5')
      n_N2O5 = n
          ntm_power(n) = -12
          tr_mm(n) = 108.02d0

#ifdef SHINDELL_STRAT_CHEM
      case ('ClOx')
      n_ClOx = n
          ntm_power(n) = -11
          tr_mm(n) = 51.5d0
C         Interpolate ClOx altitude-dependence to model resolution:
          CALL LOGPINT(LCOalt,PCOalt,ClOxaltIN,LM,PRES,ClOxalt,.true.)

      case ('BrOx')
      n_BrOx = n
          ntm_power(n) = -14
          tr_mm(n) = 95.9d0
C         Interpolate BrOx altitude-dependence to model resolution:
          CALL LOGPINT(LCOalt,PCOalt,BrOxaltIN,LM,PRES,BrOxalt,.true.)

      case ('HCl')
      n_HCl = n
          ntm_power(n) = -10
          tr_mm(n) = 36.5d0
C         Interpolate HCl altitude-dependence to model resolution:
          CALL LOGPINT(LCOalt,PCOalt,HClaltIN,LM,PRES,HClalt,.true.)

      case ('ClONO2')
      n_ClONO2 = n
          ntm_power(n) = -11
          tr_mm(n) = 97.5d0
C         Interpolate ClONO2 altitude-dependence to model resolution:
          CALL
     &    LOGPINT(LCOalt,PCOalt,ClONO2altIN,LM,PRES,ClONO2alt,.true.)
#endif

      case ('HOCl')
      n_HOCl = n
          ntm_power(n) = -14
          tr_mm(n) = 52.5d0

      case ('HBr')
      n_HBr = n
          ntm_power(n) = -14
          tr_mm(n) = 80.9d0

      case ('HOBr')
      n_HOBr = n
          ntm_power(n) = -14
          tr_mm(n) = 96.9d0

      case ('BrONO2')
      n_BrONO2 = n
          ntm_power(n) = -11
          tr_mm(n) = 141.9d0

      case ('CFC')
      n_CFC = n
          ntm_power(n) = -12
          tr_mm(n) = 137.4d0 !CFC11
#ifdef SHINDELL_STRAT_CHEM
          call openunit('CFC_IC',iu_data,.true.,.true.)
          read (iu_data) title,CFCICIN
          call closeunit(iu_data)
          write(6,*) title,' read from CFC_IC'
          do j=J_0,J_1  ; do i=1,im
           CFCICINL(:)=CFCICIN(I,J,:)
           CALL LOGPINT(LCOalt,PCOalt,CFCICINL,LM,PRES,CFCICL,.true.)
           CFCIC(I,J,:)=CFCICL(:)
          end do     ; end do
#endif

      case ('HNO3')
      n_HNO3 = n
          ntm_power(n) = -11
          tr_mm(n) = 63.018d0
          tr_RKD(n) = 2.073d3 ! in mole/J = 2.1d5 mole/(L atm)
#ifdef TRACERS_DRYDEP
          HSTAR(n)=1.d14
#endif

      case ('H2O2')
      n_H2O2 = n
          ntm_power(n) = -11
          tr_mm(n) = 34.016d0
          tr_RKD(n) = 9.869d2    ! in mole/J = 1.d5 mole/(L atm)
          tr_DHD(n) = -5.52288d4 ! in J/mole = -13.2 kcal/mole.
#ifdef TRACERS_DRYDEP
          HSTAR(n)=tr_RKD(n)*convert_HSTAR
          F0(n) = 1.d0
#endif

      case ('CH3OOH')
      n_CH3OOH = n
          ntm_power(n) = -11
          tr_mm(n) = 48.042d0
#ifdef TRACERS_DRYDEP
          HSTAR(n) = 3.d2
#endif

      case ('HCHO')
      n_HCHO = n
          ntm_power(n) = -11
          tr_mm(n) = 30.026d0
          tr_RKD(n) = 6.218d1 ! mole/J = 6.3d3 mole/(L atm)
#ifdef TRACERS_DRYDEP
          HSTAR(n)=6.d3
#endif

      case ('HO2NO2')
      n_HO2NO2 = n
          ntm_power(n) = -12
          tr_mm(n) = 79.018d0

      case ('CO')
      n_CO = n
          ntm_power(n) = -8
#ifdef EDGAR_HYDE_SOURCES
          ntsurfsrc(n) = 6
#else
          ntsurfsrc(n) = 2 ! industrial + biomass burning
#endif
          tr_mm(n) = 28.01d0
C         Interpolate CO altitude-dependence to model resolution:
          CALL LOGPINT(LCOalt,PCOalt,COaltIN,LM,PRES,COalt,.true.)

      case ('PAN')
      n_PAN = n
          ntm_power(n) = -11
          tr_mm(n) = 121.054d0   ! assuming CH3COOONO2 = PAN
#ifdef TRACERS_DRYDEP
          HSTAR(n) = 3.6d0
#endif

      case ('Isoprene')
      n_Isoprene = n
          ntm_power(n) = -11
          ntsurfsrc(n) = 1   ! (vegetation?)
          tr_mm(n) = 60.05d0 ! i.e. 5 carbons
#ifdef TRACERS_DRYDEP
          HSTAR(n) = 1.3d-2
#endif

      case ('AlkylNit')
      n_AlkylNit = n
          ntm_power(n) = -11
          tr_mm(n) = mair !unknown molecular weight, so use air and make
                          ! note in the diagnostics write-out...

      case ('Alkenes')
      n_Alkenes = n
          ntm_power(n) = -10
#ifdef EDGAR_HYDE_SOURCES          
          ntsurfsrc(n) = 8
#else          
          ntsurfsrc(n) = 3   ! industrial + biomass burning + vegetation
#endif

          tr_mm(n) = 46.59d0 ! i.e. carbon mass, weighted as:
C 68.6% indust. x (41.2% Propene(3C) + 58.8% other Alkenes/ynes(4.8C)) +
C 31.2% biomass burn x (47.9% Propene + 52.1% other Alkenes/ynes)= 3.88C
C This number wasn't adjusted when the vegetation was added.

      case ('Paraffin')
      n_Paraffin = n
          ntm_power(n) = -10
#ifdef EDGAR_HYDE_SOURCES          
          ntsurfsrc(n) = 8
#else          
          ntsurfsrc(n) = 3   ! industrial + biomass burning + vegetation
#endif
          tr_mm(n) = 59.50d0 ! i.e. carbon mass, weighted as:
C 94.2%indust. x (12.% Ethane(2C) + 11.1% Propane(3C) + 20.5% Butane(4C)
C + 18.% Pentane(5C) + 34.% higher Alkanes(7.5C) + 4.4% Ketones(4.6C)) +
C 5.8% biomass burning x (51.9% Ethane(2C) + 15.1% Propane(3C)
C + 4.5% Butane(4C) + 10.8% Pentane(5C) + 12.6% higher alkanes(8C)
C + 5.1% Ketones(3.6C)) = 4.95 C
C This number wasn't adjusted when the vegetation source was added.

#ifdef regional_Ox_tracers
      case ('OxREG1')
      n_OxREG1 = n
          ntemp= ntemp+1
          ntm_power(n) = -8
          tr_mm(n) = 48.d0
#ifdef TRACERS_DRYDEP
          F0(n) = 1.8d0
          HSTAR(n) = 1.d-2
#endif
          regOx_n(ntemp)=90.d0  !deg
          regOx_s(ntemp)=30.d0  !deg
          regOx_w(ntemp)=-180.d0!deg
          regOx_e(ntemp)= 180.d0!deg
          regOx_t(ntemp)=710.d0 !hPa
          regOx_b(ntemp)=2000.d0 ! intentionally large

      case ('OxREG2')
      n_OxREG2 = n
          ntemp= ntemp+1
          ntm_power(n) = -8
          tr_mm(n) = 48.d0
#ifdef TRACERS_DRYDEP
          F0(n) = 1.8d0
          HSTAR(n) = 1.d-2
#endif
          regOx_n(ntemp)=30.d0
          regOx_s(ntemp)=-30.d0
          regOx_w(ntemp)=-180.d0
          regOx_e(ntemp)= 180.d0
          regOx_t(ntemp)=710.d0
          regOx_b(ntemp)=2000.d0 ! intentionally large

      case ('OxREG3')
      n_OxREG3 = n
          ntemp= ntemp+1
          ntm_power(n) = -8
          tr_mm(n) = 48.d0
#ifdef TRACERS_DRYDEP
          F0(n) = 1.8d0
          HSTAR(n) = 1.d-2
#endif
          regOx_n(ntemp)=-30.d0
          regOx_s(ntemp)=-90.d0
          regOx_w(ntemp)=-180.d0
          regOx_e(ntemp)= 180.d0
          regOx_t(ntemp)=710.d0
          regOx_b(ntemp)=2000.d0 ! intentionally large

      case ('OxREG4')
      n_OxREG4 = n
          ntemp= ntemp+1
          ntm_power(n) = -8
          tr_mm(n) = 48.d0
#ifdef TRACERS_DRYDEP
          F0(n) = 1.8d0
          HSTAR(n) = 1.d-2
#endif
          regOx_n(ntemp)=90.d0
          regOx_s(ntemp)=30.d0
          regOx_w(ntemp)=-180.d0
          regOx_e(ntemp)= 180.d0
          regOx_t(ntemp)=245.d0
          regOx_b(ntemp)=710.d0

      case ('OxREG5')
      n_OxREG5 = n
          ntemp= ntemp+1
          ntm_power(n) = -8
          tr_mm(n) = 48.d0
#ifdef TRACERS_DRYDEP
          F0(n) = 1.8d0
          HSTAR(n) = 1.d-2
#endif
          regOx_n(ntemp)=30.d0
          regOx_s(ntemp)=-30.d0
          regOx_w(ntemp)=-180.d0
          regOx_e(ntemp)= 180.d0
          regOx_t(ntemp)=150.d0
          regOx_b(ntemp)=710.d0

      case ('OxREG6')
      n_OxREG6 = n
          ntemp= ntemp+1
          ntm_power(n) = -8
          tr_mm(n) = 48.d0
#ifdef TRACERS_DRYDEP
          F0(n) = 1.8d0
          HSTAR(n) = 1.d-2
#endif
          regOx_n(ntemp)=-30.d0
          regOx_s(ntemp)=-90.d0
          regOx_w(ntemp)=-180.d0
          regOx_e(ntemp)= 180.d0
          regOx_t(ntemp)=245.d0
          regOx_b(ntemp)=710.d0
#endif
#endif

      case ('DMS')
      n_DMS = n
          ntm_power(n) = -12
! the ocean source of DMS is actually interactive and therefore should
! not count for ntsurfsrc....
          ntsurfsrc(n) = 0   ! ocean DMS concentration
          ntisurfsrc(n)=1.
          tr_mm(n) = 62.
          needtrs(n)=.true.

      case ('MSA')
      n_MSA = n
          ntm_power(n) = -13
          tr_mm(n) = 96.  !(H2O2 34;SO2 64)
          trpdens(n)=1.7d3   !kg/m3 this is sulfate value
          trradius(n)=5.d-7 !m (SO4 3;BC 1;OC 3)
          fq_aer(n)=1.0   !fraction of aerosol that dissolves
          tr_wd_TYPE(n) = nPART

      case ('SO2')
      n_SO2 = n
          ntm_power(n) = -11
          ntsurfsrc(n) = 1   !Industrial
          tr_mm(n) = 64.
          tr_RKD(n) =0.0118d0 !mole/J or  1.2  M/atm
          tr_DHD(n) =-2.62d4! in J/mole= -6.27 kcal/mol
          tr_wd_TYPE(n) = nGAS
#ifdef TRACERS_DRYDEP
c         HSTAR(n)=tr_RKD(n)*convert_HSTAR
          HSTAR(N)=1.D5
#endif
      case ('SO4')
      n_SO4 = n
          ntm_power(n) = -11
          ntsurfsrc(n) = 1   ! EDGAR
          tr_mm(n) = 96.
          trpdens(n)=1.7d3   !kg/m3 this is sulfate value
          trradius(n)=3.d-7 !m
          fq_aer(n)=1.   !fraction of aerosol that dissolves
          tr_wd_TYPE(n) = nPART

#ifdef TRACERS_HETCHEM
      case ('SO4_d1')
      n_SO4_d1 = n
          ntm_power(n) = -11
          ntsurfsrc(n) = 0
          tr_mm(n) = 96.   !!!! Sulfat
          trpdens(n)=2.5d3   !kg/m3 this is clay density
          trradius(n)=0.75D-06 !m
          fq_aer(n)=1.   !fraction of aerosol that dissolves
          tr_wd_TYPE(n) = nPART
      case ('SO4_d2')
      n_SO4_d2 = n
          ntm_power(n) = -11
          ntsurfsrc(n) = 0
          tr_mm(n) = 96.
          trpdens(n)=2.65d3   !kg/m3 this is silt1 value
          trradius(n)=2.2D-06 !m
          fq_aer(n)=1.   !fraction of aerosol that dissolves
          tr_wd_TYPE(n) = nPART
      case ('SO4_d3')
      n_SO4_d3 = n
          ntm_power(n) = -11
          ntsurfsrc(n) = 0
          tr_mm(n) = 96.
          trpdens(n)=2.65d3   !this is silt2 value
          trradius(n)=4.4D-06 !m this is silt2 value
          fq_aer(n)=1.   !fraction of aerosol that dissolves
          tr_wd_TYPE(n) = nPART
      case ('SO4_d4')
      n_SO4_d4 = n
          ntm_power(n) = -11
          ntsurfsrc(n) = 0
          tr_mm(n) = 96.
          trpdens(n)=2.65d3   !kg/m3 this is silt3 value
          trradius(n)=6.7D-06 !m this is silt3 value
          fq_aer(n)=1.   !fraction of aerosol that dissolves
          tr_wd_TYPE(n) = nPART
      case ('SO4_s1')
      n_SO4_s1 = n
          ntm_power(n) = -11
          ntsurfsrc(n) = 0
          tr_mm(n) = 96.   !!!! Sulfat
          trpdens(n)=1.7d3   !kg/m3 
          trradius(n)=3.d-7 !m
          fq_aer(n)=1.   !fraction of aerosol that dissolves
          tr_wd_TYPE(n) = nPART
      case ('SO4_s2')
      n_SO4_s2 = n
          ntm_power(n) = -11
          ntsurfsrc(n) = 0
          tr_mm(n) = 96.   !!!! Sulfat
          trpdens(n)=1.7d3   !kg/m3 this is clay density
          trradius(n)=3.D-7 !m
          fq_aer(n)=1.   !fraction of aerosol that dissolves
          tr_wd_TYPE(n) = nPART
#endif

      case ('BCII')  !Insoluble industrial BC
      n_BCII = n
          ntm_power(n) = -12
          ntsurfsrc(n) = 1    
          tr_mm(n) = 12.
          trpdens(n)=1.3d3   !kg/m3 
          trradius(n)=1.d-7 !m
          fq_aer(n)=0.0   !fraction of aerosol that dissolves
          tr_wd_TYPE(n) = nPART
      case ('BCIA') !Soluble (aged) industrial BC
      n_BCIA = n
          ntm_power(n) = -12
          ntsurfsrc(n) = 0    
          tr_mm(n) = 12.
          trpdens(n)=1.3d3   !kg/m3 
          trradius(n)=1.d-7 !m
          fq_aer(n)=1.   !fraction of aerosol that dissolves
          tr_wd_TYPE(n) = nPART
      case ('BCB') !Biomass BC
      n_BCB = n
          ntm_power(n) = -12
          ntsurfsrc(n) = 0    
          tr_mm(n) = 12.
          trpdens(n)=1.3d3   !kg/m3 
          trradius(n)=1.d-7 !m
          fq_aer(n)=0.6   !fraction of aerosol that dissolves
          tr_wd_TYPE(n) = nPART
      case ('OCII') !Insoluble industrial organic mass
      n_OCII = n
          ntm_power(n) = -11
          ntsurfsrc(n) = 2 !Terpene and industrial emissions    
          tr_mm(n) = 15.6
          trpdens(n)=1.5d3   !kg/m3 
          trradius(n)=3.d-7 !m
          fq_aer(n)=0.0   !fraction of aerosol that dissolves
          tr_wd_TYPE(n) = nPART
      case ('OCIA') !Aged industrial organic mass
      n_OCIA = n
          ntm_power(n) = -11
          ntsurfsrc(n) = 0    
          tr_mm(n) = 15.6
          trpdens(n)=1.5d3   !kg/m3 
          trradius(n)=3.d-7 !m
          fq_aer(n)=1.   !fraction of aerosol that dissolves
          tr_wd_TYPE(n) = nPART
      case ('OCB') !Insoluble industrial organic mass
      n_OCB = n
          ntm_power(n) = -11
          ntsurfsrc(n) = 0    
          tr_mm(n) = 15.6
          trpdens(n)=1.5d3   !kg/m3 
          trradius(n)=3.d-7 !m
          fq_aer(n)=0.8   !fraction of aerosol that dissolves
          tr_wd_TYPE(n) = nPART

      case ('Be7')
      n_Be7 = n
          ntm_power(n) = -23        ! power of ten for tracer
          tr_mm(n) = 7.d0
          trdecay(n) =  1.51d-7
          trpdens(n) = 1.7d3    !kg/m3 this is SO4 value
          trradius(n) = 1.d-7  !appropriate for stratosphere
          fq_aer(n)=1.   !fraction of aerosol that dissolves
          tr_wd_TYPE(n) = nPART ! same as SO4

      case ('Be10')
      n_Be10 = n
          ntm_power(n) = -23
          tr_mm(n) = 10.d0
          trpdens(n) = 1.7d3   !kg/m3 this is SO4 value
          trradius(n) = 1.d-7  !appropriate for stratosphere
          fq_aer(n)=1.   !fraction of aerosol that dissolves
          tr_wd_TYPE(n) = nPART ! same as SO4

      case ('Pb210')
          n_Pb210 = n
          ntm_power(n) = -22
          tr_mm(n) = 210.d0
          trdecay(n) = 9.85d-10
          trpdens(n) = 1.7d3    !kg/m3 this is SO4 value
          trradius(n) = 3.d-7  !again S04 value
          fq_aer(n)=1.   !fraction of aerosol that dissolves
          tr_wd_TYPE(n) = nPART ! same as SO4

      case ('H2O2_s')
      n_H2O2_s = n
          ntm_power(n) = -10
          tr_mm(n) = 34.016
          tr_RKD(n) = 986.9 
          tr_DHD(n) = -5.52288d4 ! in J/mole = -13.2 kcal/mole.
          tr_wd_TYPE(n) = nGAS
#ifdef TRACERS_DRYDEP
          HSTAR(n)=tr_RKD(n)*convert_HSTAR
          F0(n) = 1.d0
#endif
      case ('seasalt1')
      n_seasalt1 = n
          ntsurfsrc(n) = 0   ! ocean bubbles
          ntisurfsrc(n)=1.
          ntm_power(n) = -10
          tr_mm(n) = 75. !Na x 3.256
          trpdens(n)=2.2d3  !kg/m3 This is for non-hydrated
          trradius(n)=4.4d-7 ! This is non-hydrated
          fq_aer(n)=1.0   !fraction of aerosol that dissolves
          tr_wd_TYPE(n) = nPART
      case ('seasalt2')
      n_seasalt2 = n
          ntsurfsrc(n) = 0   ! ocean bubbles
          ntisurfsrc(n)=1.
          ntm_power(n) = -9
          tr_mm(n) = 75. !Na x 3.256
          trpdens(n)=2.2d3  !kg/m3 This is for non-hydrated
          trradius(n)=5.0d-6                 ! This is non-hydrated
#ifdef TRACERS_AEROSOLS_Koch
          if (imAER.eq.0) trradius(n)=1.7d-6 ! This is non-hydrated
#endif
          fq_aer(n)=1.0   !fraction of aerosol that dissolves
          tr_wd_TYPE(n) = nPART

#ifdef TRACERS_DUST
      CASE('Clay')
      n_clay=n
          ntm_power(n)=-9
          trpdens(n)=2.5d3
#ifdef TRACERS_DRYDEP
          trradius(n)=0.46D-06
#endif
          fq_aer(n)=0.
          tr_wd_TYPE(n)=nPART
          tr_mm(n) = 1.
      CASE('Silt1')
      n_silt1=n
          ntm_power(n)=-9
          trpdens(n)=2.65d3
#ifdef TRACERS_DRYDEP
          trradius(n)=1.47D-06
#endif
          fq_aer(n)=0.
          tr_wd_TYPE(n)=nPART
          tr_mm(n) = 1.
      CASE('Silt2')
      n_silt2=n
          ntm_power(n)=-9
          trpdens(n)=2.65d3
#ifdef TRACERS_DRYDEP
          trradius(n)=2.94D-06
#endif
          fq_aer(n)=0.
          tr_wd_TYPE(n)=nPART
          tr_mm(n) = 1.
      CASE('Silt3')
      n_silt3=n
          ntm_power(n)=-9
          trpdens(n)=2.65d3
#ifdef TRACERS_DRYDEP
          trradius(n)=5.88D-06
#endif
          fq_aer(n)=0.
          tr_wd_TYPE(n)=nPART
          tr_mm(n) = 1.
      CASE('Silt4')
      n_silt4=n
          ntm_power(n)=-9   
          trpdens(n)=2.65d3
#ifdef TRACERS_DRYDEP
          trradius(n)=11.77D-06
#endif
          fq_aer(n)=0.
          tr_wd_TYPE(n)=nPART
          tr_mm(n) = 1.
#else
#ifdef TRACERS_MINERALS
      CASE('ClayIlli')          ! http://webmineral.com/data/Illite.shtml
      n_clayilli=n
          ntm_power(n)=-9
          trpdens(n)=2.61D3
#ifdef TRACERS_DRYDEP
          trradius(n)=0.46D-06
#endif
          fq_aer(n)=0.
          tr_wd_TYPE(n)=nPART
          tr_mm(n) = 1.
      CASE('ClayKaol')       ! http://www.webmineral.com/data/Kaolinite.shtml
      n_claykaol=n
          ntm_power(n)=-9
          trpdens(n)=2.63D3
#ifdef TRACERS_DRYDEP
          trradius(n)=0.46D-06
#endif
          fq_aer(n)=0.
          tr_wd_TYPE(n)=nPART
          tr_mm(n) = 1.
      CASE('ClaySmec')       ! http://www.webmineral.com/data/Rectorite.shtml
      n_claysmec=n
          ntm_power(n)=-9
          trpdens(n)=2.35D3     ! for Montmorillonite
#ifdef TRACERS_DRYDEP
          trradius(n)=0.46D-06
#endif
          fq_aer(n)=0.
          tr_wd_TYPE(n)=nPART
          tr_mm(n) = 1.
      CASE('ClayCalc')       ! http://www.webmineral.com/data/Calcite.shtml
      n_claycalc=n
          ntm_power(n)=-9
          trpdens(n)=2.71D3
#ifdef TRACERS_DRYDEP
          trradius(n)=0.46D-06
#endif
          fq_aer(n)=0.
          tr_wd_TYPE(n)=nPART
          tr_mm(n) = 1.
      CASE('ClayQuar')       ! http://www.webmineral.com/data/Quartz.shtml
      n_clayquar=n
          ntm_power(n)=-9
          trpdens(n)=DenQuarz   ! 2.62D3
#ifdef TRACERS_DRYDEP
          trradius(n)=0.46D-06
#endif
          fq_aer(n)=0.
          tr_wd_TYPE(n)=nPART
          tr_mm(n) = 1.
      CASE('Sil1Quar')       ! http://www.webmineral.com/data/Quartz.shtml
      n_sil1quar=n
          ntm_power(n)=-9
          trpdens(n)=DenQuarz   ! 2.62D3
#ifdef TRACERS_DRYDEP
          trradius(n)=1.47D-06
#endif
          fq_aer(n)=0.
          tr_wd_TYPE(n)=nPART
          tr_mm(n) = 1.
      CASE('Sil1Feld')       ! http://www.mindat.org/min-1624.html
      n_sil1feld=n
          ntm_power(n)=-9
          trpdens(n)=2.65D3     ! assumed, varies strongly among types
#ifdef TRACERS_DRYDEP
          trradius(n)=1.47D-06
#endif
          fq_aer(n)=0.
          tr_wd_TYPE(n)=nPART
          tr_mm(n) = 1.
      CASE('Sil1Calc')       ! http://www.webmineral.com/data/Calcite.shtml
      n_sil1calc=n
          ntm_power(n)=-9
          trpdens(n)=2.71D3
#ifdef TRACERS_DRYDEP
          trradius(n)=1.47D-06
#endif
          fq_aer(n)=0.
          tr_wd_TYPE(n)=nPART
          tr_mm(n) = 1.
      CASE('Sil1Hema')       ! http://www.webmineral.com/data/Hematite.shtml
      n_sil1hema=n
          ntm_power(n)=-9
          trpdens(n)=DenHema    ! 5.3D3
#ifdef TRACERS_DRYDEP
          trradius(n)=1.47D-06
#endif
          fq_aer(n)=0.
          tr_wd_TYPE(n)=nPART
          tr_mm(n) = 1.
      CASE('Sil1Gyps')       ! http://www.webmineral.com/data/Gypsum.shtml
      n_sil1gyps=n
          ntm_power(n)=-9
          trpdens(n)=2.3D3
#ifdef TRACERS_DRYDEP
          trradius(n)=1.47D-06
#endif
          fq_aer(n)=0.
          tr_wd_TYPE(n)=nPART
          tr_mm(n) = 1.
      CASE('Sil2Quar')       ! http://www.webmineral.com/data/Quartz.shtml
      n_sil2quar=n
          ntm_power(n)=-9
          trpdens(n)=DenQuarz   ! 2.62D3
#ifdef TRACERS_DRYDEP
          trradius(n)=2.94D-06
#endif
          fq_aer(n)=0.
          tr_wd_TYPE(n)=nPART
          tr_mm(n) = 1.
      CASE('Sil2Feld')       ! http://www.mindat.org/min-1624.html
      n_sil2feld=n
          ntm_power(n)=-9
          trpdens(n)=2.65D3     ! assumed, varies strongly among types
#ifdef TRACERS_DRYDEP
          trradius(n)=2.94D-06
#endif
          fq_aer(n)=0.
          tr_wd_TYPE(n)=nPART
          tr_mm(n) = 1.
      CASE('Sil2Calc')       ! http://www.webmineral.com/data/Calcite.shtml
      n_sil2calc=n
          ntm_power(n)=-9
          trpdens(n)=2.71D3
#ifdef TRACERS_DRYDEP
          trradius(n)=2.94D-06
#endif
          fq_aer(n)=0.
          tr_wd_TYPE(n)=nPART
          tr_mm(n) = 1.
      CASE('Sil2Hema')       ! http://www.webmineral.com/data/Hematite.shtml
      n_sil2hema=n
          ntm_power(n)=-9
          trpdens(n)=DenHema    ! 5.3D3
#ifdef TRACERS_DRYDEP
          trradius(n)=2.94D-06
#endif
          fq_aer(n)=0.
          tr_wd_TYPE(n)=nPART
          tr_mm(n) = 1.
      CASE('Sil2Gyps')       ! http://www.webmineral.com/data/Gypsum.shtml
      n_sil2gyps=n
          ntm_power(n)=-9
          trpdens(n)=2.3D3
#ifdef TRACERS_DRYDEP
          trradius(n)=2.94D-06
#endif
          fq_aer(n)=0.
          tr_wd_TYPE(n)=nPART
          tr_mm(n) = 1.
      CASE('Sil3Quar')       ! http://www.webmineral.com/data/Quartz.shtml
      n_sil3quar=n
          ntm_power(n)=-9
          trpdens(n)=DenQuarz   ! 2.62D3
#ifdef TRACERS_DRYDEP
          trradius(n)=5.88D-06
#endif
          fq_aer(n)=0.
          tr_wd_TYPE(n)=nPART
          tr_mm(n) = 1.
      CASE('Sil3Feld')       ! http://www.mindat.org/min-1624.html
      n_sil3feld=n
          ntm_power(n)=-9
          trpdens(n)=2.65D3     ! assumed, varies strongly among types
#ifdef TRACERS_DRYDEP
          trradius(n)=5.88D-06
#endif
          fq_aer(n)=0.
          tr_wd_TYPE(n)=nPART
          tr_mm(n) = 1.
      CASE('Sil3Calc')       ! http://www.webmineral.com/data/Calcite.shtml
      n_sil3calc=n
          ntm_power(n)=-9
          trpdens(n)=2.71D3
#ifdef TRACERS_DRYDEP
          trradius(n)=5.88D-06
#endif
          fq_aer(n)=0.
          tr_wd_TYPE(n)=nPART
          tr_mm(n) = 1.
      CASE('Sil3Hema')       ! http://www.webmineral.com/data/Hematite.shtml
      n_sil3hema=n
          ntm_power(n)=-9
          trpdens(n)=DenHema    ! 5.3D3
#ifdef TRACERS_DRYDEP
          trradius(n)=5.88D-06
#endif
          fq_aer(n)=0.
          tr_wd_TYPE(n)=nPART
          tr_mm(n) = 1.
      CASE('Sil3Gyps')       ! http://www.webmineral.com/data/Gypsum.shtml
      n_sil3gyps=n
          ntm_power(n)=-9
          trpdens(n)=2.3D3
#ifdef TRACERS_DRYDEP
          trradius(n)=5.88D-06
#endif
          fq_aer(n)=0.
          tr_wd_TYPE(n)=nPART
          tr_mm(n) = 1.
#endif
#ifdef TRACERS_QUARZHEM
      CASE('Sil1QuHe')
      n_sil1quhe=n
          ntm_power(n)=-9
          trpdens(n)=(1-FrHeQu)*DenQuarz+FrHeQu*DenHema
#ifdef TRACERS_DRYDEP
          trradius(n)=1.47D-06
#endif
          fq_aer(n)=0.
          tr_wd_TYPE(n)=nPART
          tr_mm(n) = 1.
      CASE('Sil2QuHe')
      n_sil2quhe=n
          ntm_power(n)=-9
          trpdens(n)=(1-FrHeQu)*DenQuarz+FrHeQu*DenHema
#ifdef TRACERS_DRYDEP
          trradius(n)=2.94D-06
#endif
          fq_aer(n)=0.
          tr_wd_TYPE(n)=nPART
          tr_mm(n) = 1.
      CASE('Sil3QuHe')
      n_sil3quhe=n
          ntm_power(n)=-9
          trpdens(n)=(1-FrHeQu)*DenQuarz+FrHeQu*DenHema
#ifdef TRACERS_DRYDEP
          trradius(n)=5.88D-06
#endif
          fq_aer(n)=0.
          tr_wd_TYPE(n)=nPART
          tr_mm(n) = 1.
#endif
#endif

#endif
      end select

#if (defined TRACERS_WATER) 
C**** Tracers that are soluble or are scavenged or are water => wet dep
      if (tr_wd_TYPE(n).eq.nWater.or.tr_wd_TYPE(n) .EQ. nPART .or.
     *  tr_RKD(n).gt.0) then
        dowetdep(n)=.true.
      end if
#endif
#ifdef TRACERS_DRYDEP
C**** If tracers are particles or have non-zero HSTAR or F0 do dry dep:
C**** Any tracers that dry deposits needs the surface concentration:
      if(HSTAR(n).GT.0..OR.F0(n).GT.0..OR.tr_wd_TYPE(n).eq.nPART) then
        dodrydep(n)=.true.
        needtrs(n)=.true.
#ifdef TRACERS_WATER
        if (tr_wd_TYPE(n).eq.nWATER) call stop_model
     &       ('A water tracer should not undergo dry deposition.',255)
#endif
      end if
#endif
#ifdef TRACERS_SPECIAL_Shindell
C**** Define the conversion from mass to volume units here so it is not
C**** done each hour:
      mass2vol(n)  =mair/TR_MM(n)
      bymass2vol(n)=TR_MM(n)/mair
#endif

      end do

C**** DIAGNOSTIC DEFINTIONS 
C**** Please note that short names for diags i.e. sname_jls are used
C**** in special ways and MUST NOT contain spaces, commas or % signs.
C**** Underscores and minus signs are allowed.

#ifdef TRACERS_ON
C**** Tracer sources and sinks
C**** Defaults for jls (sources, sinks, etc.)
C**** These need to be 'hand coded' depending on circumstances
      do k=1,ktajls  ! max number of sources and sinks
        jgrid_jls(k) = 1
        jwt_jls(k) = 1
        ia_jls(k) = ia_src
        scale_jls(k) = 1./DTsrc
      end do
      jls_grav=0
#ifdef TRACERS_WATER
C**** set defaults for some precip/wet-dep related diags
      jls_prec(:,:)=0    
#endif

C****
C**** Set some diags that are the same regardless
      call set_generic_tracer_diags
C****

      k = 0
      do n=1,ntm
      select case (trname(n))
      case ('SF6')
        k = k + 1
        jls_source(1,n) = k
        sname_jls(k) = 'Layer_1_source_of_'//trname(n)
        lname_jls(k) = 'SF6 CFC-GRID SOURCE, LAYER 1'
        jls_ltop(k) = 1
        jls_power(k) = -3.
        units_jls(k) = unit_string(jls_power(k),'kg/s')

      case ('Rn222')
        k = k + 1
        jls_decay(n) = k   ! special array for all radioactive sinks
        sname_jls(k) = 'Decay_of_'//trname(n)
        lname_jls(k) = 'LOSS OF RADON-222 BY DECAY'
        jls_ltop(k) = lm
        jls_power(k) = -26
        units_jls(k) = unit_string(jls_power(k),'kg/s/mb/m^2')
        jwt_jls(k)=3

        k = k + 1
        jls_source(1,n) = k
        sname_jls(k) = 'Ground_Source_of_'//trname(n)
        lname_jls(k) = 'RADON-222 SOURCE, LAYER 1'
        jls_ltop(k) = 1
        jls_power(k) = -10.
        units_jls(k) = unit_string(jls_power(k),'kg/s')

! keep AIJ and AJL CO2 sources in same order !!
      case ('CO2')
        k = k + 1
        jls_source(1,n) = k
        sname_jls(k) = 'Fossil_fuel_source_'//trname(n)
        lname_jls(k) = 'CO2 Fossil fuel source (Marland)'
        jls_ltop(k) = 1
        jls_power(k) = 3
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_source(2,n) = k
        sname_jls(k) = 'fertilization_sink_'//trname(n)
        lname_jls(k) = 'CO2 fertilization sink (Friedlingstein)'
        jls_ltop(k) = 1
        jls_power(k) = 3
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_source(3,n) = k
        sname_jls(k) = 'Northern_forest_regrowth_'//trname(n)
        lname_jls(k) = 'CO2 Northern forest regrowth sink'
        jls_ltop(k) = 1
        jls_power(k) = 3
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_source(4,n) = k
        sname_jls(k) = 'Land_Use_Modification_'//trname(n)
        lname_jls(k) = 'CO2 from Land use modification (Houton)'
        jls_ltop(k) = 1
        jls_power(k) = 3
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_source(5,n) = k
        sname_jls(k) = 'Ecosystem_exchange_'//trname(n)
        lname_jls(k) = 'CO2 Ecosystem exchange (Matthews)'
        jls_ltop(k) = 1
        jls_power(k) = 3
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_source(6,n) = k
        sname_jls(k) = 'Ocean_exchange_'//trname(n)
        lname_jls(k) = 'CO2 Ocean exchange'
        jls_ltop(k) = 1
        jls_power(k) = 3
        units_jls(k) = unit_string(jls_power(k),'kg/s')

      case ('N2O')
#ifdef TRACERS_SPECIAL_Shindell
        k = k + 1
        jls_3Dsource(1,n) = k
        sname_jls(k) = 'chemistry_source_of'//trname(n)
        lname_jls(k) = 'CHANGE OF '//trname(n)//' BY CHEMISTRY'
        jls_ltop(k) = LM
        jls_power(k) = -1.
        units_jls(k) = unit_string(jls_power(k),'kg/s')
#endif
#ifdef TRACERS_SPECIAL_Lerner
        k = k + 1
        jls_source(1,n) = k
        sname_jls(k) = 'L1_sink_'//trname(n)
        lname_jls(k) = 'CHANGE OF N20 BY RESETTING TO 462.2d-9, L1'
        jls_ltop(k) = 1
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_3Dsource(1,n) = k
        sname_jls(k) = 'Stratos_chem_change_'//trname(n)
        lname_jls(k) = 'CHANGE OF N2O BY CHEMISTRY IN STRATOS'
        jls_ltop(k) = lm
        jls_power(k) = -1
        units_jls(k) = unit_string(jls_power(k),'kg/s')
#endif

      case ('CFC11')   !!! should start April 1
        k = k + 1
        jls_source(1,n) = k
        sname_jls(k) = 'L1_sink_'//trname(n)
        lname_jls(k) = 'CHANGE OF CFC-11 BY SOURCE, L1'
        jls_ltop(k) = 1
        jls_power(k) = -1
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_3Dsource(1,n) = k
        sname_jls(k) = 'Stratos_chem_change_'//trname(n)
        lname_jls(k) = 'CHANGE OF CFC-11 BY CHEMISTRY IN STRATOS'
        jls_ltop(k) = lm
        jls_power(k) = -3
        units_jls(k) = unit_string(jls_power(k),'kg/s')

      case ('14CO2')   !!! should start 10/16
        k = k + 1
        jls_source(1,n) = k
        sname_jls(k) = 'L1_sink_'//trname(n)
        lname_jls(k) = 'CHANGE OF 14CO2 by SINK, L1'
        jls_ltop(k) = 1
        jls_power(k) = -4
        units_jls(k) = unit_string(jls_power(k),'kg/s')

      case ('CH4')
        k = k + 1
        jls_source(6,n) = k
        sname_jls(k) = 'Soil_sink_'//trname(n)
        lname_jls(k) = 'CH4 sink due to soil absorption'
        jls_ltop(k) = 1
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_source(7,n) = k
        sname_jls(k) = 'Termite_source_'//trname(n)
        lname_jls(k) = 'CH4 Termite source'
        jls_ltop(k) = 1
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_source(9,n) = k
        sname_jls(k) = 'Ocean_source_'//trname(n)
        lname_jls(k) = 'CH4 Ocean source'
        jls_ltop(k) = 1
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_source(10,n) = k
        sname_jls(k) = 'Fresh_Water_lake_source_'//trname(n)
        lname_jls(k) = 'CH4 Fresh Water lake source'
        jls_ltop(k) = 1
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_source(11,n) = k
        sname_jls(k) = 'Misc_Ground_source_'//trname(n)
        lname_jls(k) = 'CH4 Misc_Ground source'
        jls_ltop(k) = 1
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_source(14,n) = k
        sname_jls(k) = 'Wetlands+Tundra_source_'//trname(n)
        lname_jls(k) = 'CH4 Wetlands+Tundra source'
        jls_ltop(k) = 1
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),'kg/s')
#ifdef EDGAR_HYDE_SOURCES
        k = k + 1
        jls_source(1,n) = k
        sname_jls(k) = 'EH_Animal_source_'//trname(n)
        lname_jls(k) = 'CH4 EH Animal source'
        jls_ltop(k) = 1
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_source(2,n) = k
        sname_jls(k) = 'EH_Fos_Fuel_Comb_source_'//trname(n)
        lname_jls(k) = 'CH4 EH Fossil Fuel Combustion source'
        jls_ltop(k) = 1
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_source(3,n) = k
        sname_jls(k) = 'EH_Landfill_source_'//trname(n)
        lname_jls(k) = 'CH4 EH Landfills source'
        jls_ltop(k) = 1
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_source(4,n) = k
        sname_jls(k) = 'EH_Fos_Fuel_Prod_source_'//trname(n)
        lname_jls(k) = 'CH4 EH Fossil Fuel Production etc source'
        jls_ltop(k) = 1
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_source(5,n) = k
        sname_jls(k) = 'EH_Ag_waste_burn_source_'//trname(n)
        lname_jls(k) = 'CH4 EH Agricultural Waste Burning source'
        jls_ltop(k) = 1
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_source(8,n) = k
        sname_jls(k) = 'EH_Deforest_source_'//trname(n)
        lname_jls(k) = 'CH4 EH Deforestation source'
        jls_ltop(k) = 1
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_source(12,n) = k
        sname_jls(k) = 'EH_savan_burn_source_'//trname(n)
        lname_jls(k) = 'CH4 EH Savannah Burning source'
        jls_ltop(k) = 1
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_source(13,n) = k
        sname_jls(k) = 'EH_Biofuel_source_'//trname(n)
        lname_jls(k) = 'CH4 EH Biofuel P,T,C source'
        jls_ltop(k) = 1
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_source(15,n) = k
        sname_jls(k) = 'EH_Ag_Land_source_'//trname(n)
        lname_jls(k) = 'CH4 EH Agricultural Land Activities source'
        jls_ltop(k) = 1
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),'kg/s')        
#else
        k = k + 1
        jls_source(1,n) = k
        sname_jls(k) = 'Animal_source_of'//trname(n)
        lname_jls(k) = 'CH4 Animal source'
        jls_ltop(k) = 1
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_source(2,n) = k
        sname_jls(k) = 'Coal_Mine_source_'//trname(n)
        lname_jls(k) = 'CH4 Coal Mine source'
        jls_ltop(k) = 1
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_source(3,n) = k
        sname_jls(k) = 'Gas_Leak_source_'//trname(n)
        lname_jls(k) = 'CH4 Gas Leak source'
        jls_ltop(k) = 1
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_source(4,n) = k
        sname_jls(k) = 'Gas_Venting_source_'//trname(n)
        lname_jls(k) = 'CH4 Gas Venting source'
        jls_ltop(k) = 1
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_source(5,n) = k
        sname_jls(k) = 'Municipal_solid_waste_source_'//trname(n)
        lname_jls(k) = 'CH4 Municipal solid waste source'
        jls_ltop(k) = 1
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_source(8,n) = k
        sname_jls(k) = 'Coal_combustion_source_'//trname(n)
        lname_jls(k) = 'CH4 Coal combustion source'
        jls_ltop(k) = 1
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_source(12,n) = k
        sname_jls(k) = 'Biomass_burning_source_'//trname(n)
        lname_jls(k) = 'CH4 Biomass burning source'
        jls_ltop(k) = 1
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_source(13,n) = k
        sname_jls(k) = 'Rice_Cultivation_source_'//trname(n)
        lname_jls(k) = 'CH4 Rice Cultivation source'
        jls_ltop(k) = 1
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),'kg/s')
#endif
#ifdef TRACERS_SPECIAL_Shindell
        k = k + 1
        jls_3Dsource(nChemistry,n) = k
        sname_jls(k) = 'chemistry_source_of'//trname(n)
        lname_jls(k) = 'CHANGE OF CH4 BY CHEMISTRY'
        jls_ltop(k) = LTOP
        jls_power(k) = 0.
        units_jls(k) = unit_string(jls_power(k),'kg/s')
#ifndef SHINDELL_STRAT_CHEM
        k = k + 1
        jls_3Dsource(nStratwrite,n) = k
        sname_jls(k) = 'strat_overwrite_source_of'//trname(n)
        lname_jls(k) = 'CHANGE OF CH4 BY STRATOSPHERIC OVERWRITE'
        jls_ltop(k) = LM
        jls_power(k) = 0.
        units_jls(k) = unit_string(jls_power(k),'kg/s')
#endif
#else
       k = k + 1
        jls_3Dsource(1,n) = k
        sname_jls(k) = 'Tropos_Chem_change_'//trname(n)
        lname_jls(k) = 'CHANGE OF CH4 BY CHEMISTRY IN TROPOSPHERE'
        jls_ltop(k) = lm
        jls_power(k) = -1
        units_jls(k) = unit_string(jls_power(k),'kg/s')
       k = k + 1
        jls_3Dsource(2,n) = k
        sname_jls(k) = 'Stratos_Chem_change_'//trname(n)
        lname_jls(k) = 'CHANGE OF CH4 BY CHEMISTRY IN STRATOS'
        jls_ltop(k) = lm
        jls_power(k) = -1
        units_jls(k) = unit_string(jls_power(k),'kg/s')
#endif

      case ('O3')
       k = k + 1
        jls_3Dsource(1,n) = k
        sname_jls(k) = 'Tropos_Chem_change_'//trname(n)
        lname_jls(k) = 'CHANGE OF O3 BY CHEMISTRY IN TROPOSPHERE'
        jls_ltop(k) = lm
        jls_power(k) = 1
        units_jls(k) = unit_string(jls_power(k),'kg/s')
       k = k + 1
        jls_3Dsource(2,n) = k
        sname_jls(k) = 'Stratos_Chem_change_'//trname(n)
        lname_jls(k) = 'CHANGE OF O3 BY CHEMISTRY IN STRATOS'
        jls_ltop(k) = lm
        jls_power(k) = 1
        units_jls(k) = unit_string(jls_power(k),'kg/s')

#ifdef TRACERS_WATER
C**** generic ones for many water tracers
      case ('Water', 'H2O18', 'HDO', 'HTO')
       k = k + 1
        jls_source(1,n)=k
        sname_jls(k) = 'Evap_'//trname(n)
        lname_jls(k) = 'EVAPORATION OF '//trname(n)
        jls_ltop(k) = 1
        jls_power(k) = ntm_power(n)+4
        scale_jls(k) = SDAY*byim/DTsrc
        units_jls(k) = unit_string(jls_power(k),'mm/day')
       k = k + 1
        jls_source(2,n)=k
        sname_jls(k) = 'Ocn_Evap_'//trname(n)
        lname_jls(k) = 'OCEAN EVAP OF '//trname(n)
        jls_ltop(k) = 1
        jls_power(k) = ntm_power(n)+4
        scale_jls(k) = SDAY*byim/DTsrc
        units_jls(k) = unit_string(jls_power(k),'mm/day')
       k = k + 1
        jls_prec(1,n)=k
        sname_jls(k) = 'Precip_'//trname(n)
        lname_jls(k) = 'PRECIPITATION OF '//trname(n)
        jls_ltop(k) = 1
        jls_power(k) = ntm_power(n)+4
        scale_jls(k) = SDAY*byim/DTsrc
        units_jls(k) = unit_string(jls_power(k),'mm/day')
       k = k + 1
        jls_prec(2,n)=k
        sname_jls(k) = 'Ocn_Precip_'//trname(n)
        lname_jls(k) = 'OCEAN PRECIP OF '//trname(n)
        jls_ltop(k) = 1
        jls_power(k) = ntm_power(n)+4
        scale_jls(k) = SDAY*byim/DTsrc
        units_jls(k) = unit_string(jls_power(k),'mm/day')

C**** special one unique to HTO
      if (trname(n).eq."HTO") then
       k = k + 1
        jls_decay(n) = k   ! special array for all radioactive sinks
        sname_jls(k) = 'Decay_of_'//trname(n)
        lname_jls(k) = 'LOSS OF '//TRIM(trname(n))//' BY DECAY'
        jls_ltop(k) = lm
        jls_power(k) = ntm_power(n)+8
        scale_jls(k) = 1./DTsrc
        units_jls(k) = unit_string(jls_power(k),'kg/s')
      end if
#endif

      case ('NOx')
        k = k + 1
        jls_3Dsource(nChemistry,n) = k
        sname_jls(k) = 'chemistry_source_of'//trname(n)
        lname_jls(k) = 'CHANGE OF NOx BY CHEMISTRY'
        jls_ltop(k) = LTOP
        jls_power(k) = -1.
        units_jls(k) = unit_string(jls_power(k),'kg/s')
#ifndef SHINDELL_STRAT_CHEM
        k = k + 1
        jls_3Dsource(nStratwrite,n) = k
        sname_jls(k) = 'strat_overwrite_source_of'//trname(n)
        lname_jls(k) = 'CHANGE OF NOx BY STRATOSPHERIC OVERWRITE'
        jls_ltop(k) = LM
        jls_power(k) = 0.
        units_jls(k) = unit_string(jls_power(k),'kg/s')
#endif
        k = k + 1
        jls_3Dsource(nLightning,n) = k
        sname_jls(k) = 'lightning_source_of'//trname(n)
        lname_jls(k) = 'CHANGE OF NOx BY LIGHTNING'
        jls_ltop(k) = LM
        jls_power(k) = -2.
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_3Dsource(nAircraft,n) = k
        sname_jls(k) = 'aircraft_source_of'//trname(n)
        lname_jls(k) = 'CHANGE OF NOx BY AIRCRAFT'
        jls_ltop(k) = LM
        jls_power(k) = -2.
        units_jls(k) = unit_string(jls_power(k),'kg/s')
#ifdef EDGAR_HYDE_SOURCES
        k = k + 1
        jls_source(1,n) = k
        sname_jls(k) = 'EH_Fos_Fuel_Source_of_'//trname(n)
        lname_jls(k) = 'NOx EH fossil fuel combustion source'
        jls_ltop(k) = 1
        jls_power(k) = -2.
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_source(2,n) = k
        sname_jls(k) = 'EH_Industrial_Source_of_'//trname(n)
        lname_jls(k) = 'NOx EH Industrial Processes source'
        jls_ltop(k) = 1
        jls_power(k) = -2.
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_source(3,n) = k
        sname_jls(k) = 'EH_Ag_Land_source_'//trname(n)
        lname_jls(k) = 'NOx EH Agricultural Land Activitied source'
        jls_ltop(k) = 1
        jls_power(k) = -2.
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_source(4,n) = k
        sname_jls(k) = 'EH_Ag_waste_burn_source_'//trname(n)
        lname_jls(k) = 'NOx EH Agricultural Waste Burning source'
        jls_ltop(k) = 1
        jls_power(k) = -2.
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_source(5,n) = k
        sname_jls(k) = 'EH_Biofuel_source_'//trname(n)
        lname_jls(k) = 'NOx EH Biofuel P,T,C source'
        jls_ltop(k) = 1
        jls_power(k) = -2.
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_source(6,n) = k
        sname_jls(k) = 'EH_Deforest_source_'//trname(n)
        lname_jls(k) = 'NOx EH Deforestation source'
        jls_ltop(k) = 1
        jls_power(k) = -2.
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_source(7,n) = k
        sname_jls(k) = 'EH_Savannah_source_'//trname(n)
        lname_jls(k) = 'NOx EH Savannah Burning source'
        jls_ltop(k) = 1
        jls_power(k) = -2.
        units_jls(k) = unit_string(jls_power(k),'kg/s')
#else
        k = k + 1
        jls_source(1,n) = k
        sname_jls(k) = 'Fossil_Fuel_Source_of_'//trname(n)
        lname_jls(k) = 'NOx fossil fuel source'
        jls_ltop(k) = 1
        jls_power(k) = -2.
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_source(2,n) = k
        sname_jls(k) = 'Biomass_Burning_Source_of_'//trname(n)
        lname_jls(k) = 'NOx biomass burning source'
        jls_ltop(k) = 1
        jls_power(k) = -2.
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_source(3,n) = k
        sname_jls(k) = 'Soil_Source_of_'//trname(n)
        lname_jls(k) = 'NOx soil source'
        jls_ltop(k) = 1
        jls_power(k) = -2.
        units_jls(k) = unit_string(jls_power(k),'kg/s')
#endif

      case ('N2O5','HNO3','H2O2','CH3OOH','HCHO','HO2NO2','PAN',
     &      'AlkylNit','Ox','OxREG1','OxREG2','OxREG3',
     &      'OxREG4','OxREG5','OxREG6')
        k = k + 1
        jls_3Dsource(nChemistry,n) = k
        sname_jls(k) = 'chemistry_source_of'//trname(n)
        lname_jls(k) = 'CHANGE OF '//trname(n)//' BY CHEMISTRY'
        jls_ltop(k) = LTOP
        select case (trname(n))
        case ('Ox')
          jls_power(k) = 1.        
        case ('H2O2','HCHO','OxREG1','OxREG2','OxREG3','OxREG4',
     &        'OxREG5','OxREG6')
          jls_power(k) = 0.
        case default
          jls_power(k) = -1.
        end select
        units_jls(k) = unit_string(jls_power(k),'kg/s')
#ifndef SHINDELL_STRAT_CHEM
        k = k + 1
        jls_3Dsource(nStratwrite,n) = k
        sname_jls(k) = 'strat_overwrite_source_of'//trname(n)
        lname_jls(k) =
     &  'CHANGE OF '//trname(n)//' BY STRATOSPHERIC OVERWRITE'
        jls_ltop(k) = LM
        jls_power(k) = 0.
        units_jls(k) = unit_string(jls_power(k),'kg/s')
#endif

       case ('CO')
        k = k + 1
        jls_3Dsource(nChemistry,n) = k
        sname_jls(k) = 'chemistry_source_of'//trname(n)
        lname_jls(k) = 'CHANGE OF CO BY CHEMISTRY'
        jls_ltop(k) = LTOP
        jls_power(k) = 0.
        units_jls(k) = unit_string(jls_power(k),'kg/s')
#ifndef SHINDELL_STRAT_CHEM
        k = k + 1
        jls_3Dsource(nStratwrite,n) = k
        sname_jls(k) = 'strat_overwrite_source_of'//trname(n)
        lname_jls(k) = 'CHANGE OF CO BY STRATOSPHERIC OVERWRITE'
        jls_ltop(k) = LM
        jls_power(k) = 0.
        units_jls(k) = unit_string(jls_power(k),'kg/s')
#endif
#ifdef EDGAR_HYDE_SOURCES
        k = k + 1
        jls_source(1,n) = k
        sname_jls(k) = 'EH_Fos_Fuel_Source_'//trname(n)
        lname_jls(k) = 'CO EH Fossil Fuel Combustion source'
        jls_ltop(k) = 1
        jls_power(k) = -1
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_source(2,n) = k
        sname_jls(k) = 'EH_Industrial_Source_'//trname(n)
        lname_jls(k) = 'CO EH Industrial Processes source'
        jls_ltop(k) = 1
        jls_power(k) = -1
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_source(3,n) = k
        sname_jls(k) = 'EH_Ag_Waste_Source_'//trname(n)
        lname_jls(k) = 'CO EH Agricultural Waste Burning source'
        jls_ltop(k) = 1
        jls_power(k) = -1
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_source(4,n) = k
        sname_jls(k) = 'EH_Savannah_Source_'//trname(n)
        lname_jls(k) = 'CO EH Savannah Burning source'
        jls_ltop(k) = 1
        jls_power(k) = -1
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_source(5,n) = k
        sname_jls(k) = 'EH_Biofuel_Source_'//trname(n)
        lname_jls(k) = 'CO EH Biofuel P,T,C source'
        jls_ltop(k) = 1
        jls_power(k) = -1
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_source(6,n) = k
        sname_jls(k) = 'EH_Deforest_Source_'//trname(n)
        lname_jls(k) = 'CO EH Deforestation source'
        jls_ltop(k) = 1
        jls_power(k) = -1
        units_jls(k) = unit_string(jls_power(k),'kg/s')
#else
        k = k + 1
        jls_source(1,n) = k
        sname_jls(k) = 'Industrial_Source_of'//trname(n)
        lname_jls(k) = 'CO industrial source'
        jls_ltop(k) = 1
        jls_power(k) = -1
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_source(2,n) = k
        sname_jls(k) = 'Biomass_Burning_Source_of'//trname(n)
        lname_jls(k) = 'CO biomass burning source'
        jls_ltop(k) = 1
        jls_power(k) = -1
        units_jls(k) = unit_string(jls_power(k),'kg/s')
#endif

      case ('Isoprene')
        k = k + 1
        jls_3Dsource(nChemistry,n) = k
        sname_jls(k) = 'chemistry_source_of'//trname(n)
        lname_jls(k) = 'CHANGE OF Isoprene BY CHEMISTRY'
        jls_ltop(k) = LTOP
        jls_power(k) = 1.
        units_jls(k) = unit_string(jls_power(k),'kg/s')
#ifndef SHINDELL_STRAT_CHEM
        k = k + 1
        jls_3Dsource(nStratwrite,n) = k
        sname_jls(k) = 'strat_overwrite_source_of'//trname(n)
        lname_jls(k) = 'CHANGE OF Isoprene BY STRATOSPHERIC OVERWRITE'
        jls_ltop(k) = LM
        jls_power(k) = 0.
        units_jls(k) = unit_string(jls_power(k),'kg/s')
#endif
        k = k + 1
        jls_source(1,n) = k
        sname_jls(k) = 'Vegetation_source_of'//trname(n)
        lname_jls(k) = 'Isoprene vegetation source'
        jls_ltop(k) = 1
        jls_power(k) = -1
        units_jls(k) = unit_string(jls_power(k),'kg/s')

      case ('Alkenes','Paraffin')
        k = k + 1
        jls_3Dsource(nChemistry,n) = k
        sname_jls(k) = 'chem_source_of'//trname(n)
        lname_jls(k) = 'CHANGE OF '//trname(n)//' BY CHEMISTRY'
        jls_ltop(k) = LTOP
        jls_power(k) = -1.
        units_jls(k) = unit_string(jls_power(k),'kg/s')
#ifndef SHINDELL_STRAT_CHEM
        k = k + 1
        jls_3Dsource(nStratwrite,n) = k
        sname_jls(k) = 'strat_overwrite_source_of'//trname(n)
        lname_jls(k) =
     &  'CHANGE OF '//trname(n)//' BY STRATOSPHERIC OVERWRITE'
        jls_ltop(k) = LM
        jls_power(k) = 0.
        units_jls(k) = unit_string(jls_power(k),'kg/s')
#endif
        k = k + 1
        jls_source(3,n) = k
        sname_jls(k) = 'Vegetation_source_of'//trname(n)
        lname_jls(k) = trname(n)//' vegetation source'
        jls_ltop(k) = 1
        jls_power(k) = -2
        units_jls(k) = unit_string(jls_power(k),'kg/s')
#ifdef EDGAR_HYDE_SOURCES
        k = k + 1
        jls_source(1,n) = k
        sname_jls(k) = 'EH_Foss_Fuel_Comb_source_'//trname(n)
        lname_jls(k) = trname(n)//' EH Fossil Fuel Combustion source'
        jls_ltop(k) = 1
        jls_power(k) = -2
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_source(2,n) = k
        sname_jls(k) = 'EH_Industrial_source_'//trname(n)
        lname_jls(k) = trname(n)//' EH Idustrial Processes source'
        jls_ltop(k) = 1
        jls_power(k) = -2
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
!! #3 is (common) vegetation source.  See above.
        jls_source(4,n) = k
        sname_jls(k) = 'EH_Foss_Fuel_Prod_source_'//trname(n)
        lname_jls(k) = trname(n)//' EH Fossil Fuel Production source'
        jls_ltop(k) = 1
        jls_power(k) = -2
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_source(5,n) = k
        sname_jls(k) = 'EH_Ag_Waste_source_'//trname(n)
        lname_jls(k) =
     &  trname(n)//' EH Agricultural Waste Burning source'
        jls_ltop(k) = 1
        jls_power(k) = -2
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_source(6,n) = k
        sname_jls(k) = 'EH_Savanna_source_'//trname(n)
        lname_jls(k) = trname(n)//' EH Savanna Burning source'
        jls_ltop(k) = 1
        jls_power(k) = -2
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_source(7,n) = k
        sname_jls(k) = 'EH_Biofuel_source_'//trname(n)
        lname_jls(k) = trname(n)//' EH Biofuel P,T,C source'
        jls_ltop(k) = 1
        jls_power(k) = -2
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_source(8,n) = k
        sname_jls(k) = 'EH_Deforest_source_'//trname(n)
        lname_jls(k) = trname(n)//' EH Deforestation source'
        jls_ltop(k) = 1
        jls_power(k) = -2
        units_jls(k) = unit_string(jls_power(k),'kg/s')
#else
        k = k + 1
        jls_source(1,n) = k
        sname_jls(k) = 'Industrial_source_of'//trname(n)
        lname_jls(k) = trname(n)//' industrial source'
        jls_ltop(k) = 1
        jls_power(k) = -2
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_source(2,n) = k
        sname_jls(k) = 'Biomass_Burning_source_of'//trname(n)
        lname_jls(k) = trname(n)//' biomass burning source'
        jls_ltop(k) = 1
        jls_power(k) = -2
        units_jls(k) = unit_string(jls_power(k),'kg/s')
!! #3 is (common) vegetation source. See above.
#endif 

      case ('ClOx','BrOx','HCl','HOCl','ClONO2',
     &      'HBr','HOBr','BrONO2','CFC')
        k = k + 1
        jls_3Dsource(nChemistry,n) = k
        sname_jls(k) = 'chemistry_source_of'//trname(n)
        lname_jls(k) = 'CHANGE OF '//trname(n)//' BY CHEMISTRY'
        jls_ltop(k) = LM
        jls_power(k) = -1.
        units_jls(k) = unit_string(jls_power(k),'kg/s')

      case ('DMS')
        k = k + 1
        jls_isrc(1,n) = k
        sname_jls(k) = 'Ocean_source_of'//trname(n)
        lname_jls(k) = 'DMS ocean source'
        jls_ltop(k) = 1
        jls_power(k) =0
        units_jls(k) = unit_string(jls_power(k),'kg/s')
C
        k = k + 1
        jls_source(2,n) = k
        sname_jls(k) = 'TKE_Contribution'//trname(n)
        lname_jls(k) = 'SGSWSP TKE'
        jwt_jls(k) = 2
        jls_ltop(k) = 1
        jls_power(k) =0
        scale_jls(k) = 1.
        units_jls(k) = unit_string(jls_power(k),'%')

        k = k + 1
        jls_source(3,n) = k
        sname_jls(k) = 'Wet_Conv_Contr'//trname(n)
        lname_jls(k) = 'SGSWSP Wet Conv'
        jwt_jls(k) = 2
        jls_ltop(k) = 1
        jls_power(k) =0
        scale_jls(k) = 1.
        units_jls(k) = unit_string(jls_power(k),'%')

        k = k + 1
        jls_source(4,n) = k
        sname_jls(k) = 'Dry_Conv_Contr'//trname(n)
        lname_jls(k) = 'SGSWSP Dry Conv'
        jwt_jls(k) = 2
        jls_ltop(k) = 1
        jls_power(k) =0
        scale_jls(k) = 1.
        units_jls(k) = unit_string(jls_power(k),'%')

        k = k + 1
        jls_source(5,n) = k
        sname_jls(k) = 'SGSWSP-old'//trname(n)
        lname_jls(k) = 'DMS SGSWP-old/old'
        jwt_jls(k) = 2
        jls_ltop(k) = 1
        jls_power(k) =0
        scale_jls(k) = 1.
        units_jls(k) = unit_string(jls_power(k),'%')

        k = k + 1
        jls_3Dsource(1,n) = k
        sname_jls(k) = 'Chemical_sink_of'//trname(n)
        lname_jls(k) = 'DMS chemical loss'
        jls_ltop(k) =LM
        jls_power(k) =0
        units_jls(k) = unit_string(jls_power(k),'kg/s')

       case ('MSA')
c put in chemical production of MSA
        k = k + 1
        jls_3Dsource(1,n) = k
        sname_jls(k) = 'chemistry_source_of'//trname(n)
        lname_jls(k) = 'Chemical production of MSA'
        jls_ltop(k) = LM
        jls_power(k) = -1
        units_jls(k) = unit_string(jls_power(k),'kg/s')
c gravitational settling of MSA
        k = k + 1
        jls_grav(n) = k
        sname_jls(k) = 'grav_sett_of'//trname(n)
        lname_jls(k) = 'Gravitational Settling of MSA'
        jls_ltop(k) = LM
        jls_power(k) = -3
        units_jls(k) = unit_string(jls_power(k),'kg/s')

       case ('SO2')
c industrial source
        k = k + 1
        jls_source(1,n) = k
        sname_jls(k) = 'Industrial_src_of'//trname(n)
        lname_jls(k) = 'SO2 industrial source'
        jls_ltop(k) = 1
        jls_power(k) =0
        units_jls(k) = unit_string(jls_power(k),'kg/s')
c volcanic production of SO2
        k = k + 1
        jls_3Dsource(1,n) = k
        sname_jls(k) = 'volcanic_source_of'//trname(n)
        lname_jls(k) = 'production of SO2 from volcanos'
        jls_ltop(k) = LM
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),'kg/s')
c aircraft production of SO2
        k = k + 1
        jls_3Dsource(2,n) = k
        sname_jls(k) = 'aircraft_source_of'//trname(n)
        lname_jls(k) = 'production of SO2 from aircraft'
        jls_ltop(k) = LM
        jls_power(k) = -2
        units_jls(k) = unit_string(jls_power(k),'kg/s')
c biomass burning source
        k = k + 1
        jls_3Dsource(3,n) = k
        sname_jls(k) = 'Biomass_src_of'//trname(n)
        lname_jls(k) = 'SO2 biomass source'
        jls_ltop(k) = LM 
        jls_power(k) =0
        units_jls(k) = unit_string(jls_power(k),'kg/s')
c put in chemical production of SO2
        k = k + 1
        jls_3Dsource(4,n) = k
        sname_jls(k) = 'dms_source_of'//trname(n)
        lname_jls(k) = 'production of SO2 from DMS'
        jls_ltop(k) = LM
        jls_power(k) =  1
        units_jls(k) = unit_string(jls_power(k),'kg/s')
c put in chemical sink of SO2
        k = k + 1
        jls_3Dsource(5,n) = k
        sname_jls(k) = 'chem_sink_of'//trname(n)
        lname_jls(k) = 'chemical sink of SO2'
        jls_ltop(k) = LM
        jls_power(k) =  1
        units_jls(k) = unit_string(jls_power(k),'kg/s')
c heterogenous sink SO2
        k = k + 1
        jls_3Dsource(6,n) = k
        sname_jls(k) = 'het_src_of'//trname(n)
        lname_jls(k) = 'SO2 Heterogenous sink'
        jls_ltop(k) = LM
        jls_power(k) = -1
        units_jls(k) = unit_string(jls_power(k),'kg/s')
#ifdef TRACERS_AEROSOLS_Koch
c convective chem cloud phase sink of SO2
        k = k + 1
        jls_incloud(1,n) = k
        sname_jls(k) = 'mc_cloud_chem_sink_of'//trname(n)
        lname_jls(k) = 'SO2 used in convective cloud chem'
        jls_ltop(k) = LM
        jls_power(k) = -1
        units_jls(k) = unit_string(jls_power(k),'kg/s')
c stratiform chem cloud phase sink of SO2
        k = k + 1
        jls_incloud(2,n) = k
        sname_jls(k) = 'ss_cloud_chem_sink_of'//trname(n)
        lname_jls(k) = 'SO2 used in stratiform cloud chem'
        jls_ltop(k) = LM
        jls_power(k) = -1
        units_jls(k) = unit_string(jls_power(k),'kg/s')
#endif
        case ('SO4')
c gas phase source of SO4
        k = k + 1
        jls_3Dsource(1,n) = k
        sname_jls(k) = 'gas_phase_source_of'//trname(n)
        lname_jls(k) = 'SO4 gas phase source'
        jls_ltop(k) = LM
        jls_power(k) = 1
        units_jls(k) = unit_string(jls_power(k),'kg/s')
c heterogeneous source of SO4
        k = k + 1
        jls_3Dsource(2,n) = k
        sname_jls(k) = 'heterogeneous_source_of'//trname(n)
        lname_jls(k) = 'Heterogenous Chemical source of SO4'
        jls_ltop(k) = LM
        jls_power(k) = -3
        units_jls(k) = unit_string(jls_power(k),'kg/s')
c volcanic source of SO4
        k = k + 1
        jls_3Dsource(3,n) = k
        sname_jls(k) = 'volcanic_source_of'//trname(n)
        lname_jls(k) = 'SO4 volcanic source'
        jls_ltop(k) = LM
        jls_power(k) = 1
        units_jls(k) = unit_string(jls_power(k),'kg/s')
c biomass source of SO4
        k = k + 1
        jls_3Dsource(4,n) = k
        sname_jls(k) = 'biomass_source_of'//trname(n)
        lname_jls(k) = 'SO4 biomass source'
        jls_ltop(k) = LM
        jls_power(k) = 1
        units_jls(k) = unit_string(jls_power(k),'kg/s')
#ifdef TRACERS_AEROSOLS_Koch
c convective cloud phase source of SO4
        k = k + 1
        jls_incloud(1,n) = k
        sname_jls(k) = 'mc_cloud_source_of'//trname(n)
        lname_jls(k) = 'SO4 made in convective clouds'
        jls_ltop(k) = LM
        jls_power(k) = -1
        units_jls(k) = unit_string(jls_power(k),'kg/s')
c stratiform cloud phase source of SO4
        k = k + 1
        jls_incloud(2,n) = k
        sname_jls(k) = 'ss_cloud_source_of'//trname(n)
        lname_jls(k) = 'SO4 made in stratiform clouds'
        jls_ltop(k) = LM
        jls_power(k) = -1
        units_jls(k) = unit_string(jls_power(k),'kg/s')
#endif
c industrial source
        k = k + 1
        jls_source(1,n) = k
        sname_jls(k) = 'Industrial_src_of'//trname(n)
        lname_jls(k) = 'SO4 industrial source'
        jls_ltop(k) = 1
        jls_power(k) =0
        units_jls(k) = unit_string(jls_power(k),'kg/s')
c gravitational settling of SO4
        k = k + 1
        jls_grav(n) = k
        sname_jls(k) = 'grav_sett_of'//trname(n)
        lname_jls(k) = 'Gravitational Settling of SO4'
        jls_ltop(k) = LM
        jls_power(k) = -3
        units_jls(k) = unit_string(jls_power(k),'kg/s')

#ifdef TRACERS_HETCHEM
        case ('SO4_d1')
c gas phase source of SO4_d1
        k = k + 1
        jls_3Dsource(1,n) = k
        sname_jls(k) = 'gas_phase_source_of'//trname(n)
        lname_jls(k) = 'SO4_d1 gas phase source'
        jls_ltop(k) = LM
        jls_power(k) = 1
        units_jls(k) = unit_string(jls_power(k),'kg/s')
c gravitational settling of SO4_d1
        k = k + 1
        jls_grav(n) = k
        sname_jls(k) = 'grav_sett_of'//trname(n)
        lname_jls(k) = 'Gravitational Settling of SO4_d1'
        jls_ltop(k) = LM
        jls_power(k) = -3
        units_jls(k) = unit_string(jls_power(k),'kg/s')

        case ('SO4_d2')
c gas phase source of SO4_d2
        k = k + 1
        jls_3Dsource(1,n) = k
        sname_jls(k) = 'gas_phase_source_of'//trname(n)
        lname_jls(k) = 'SO4_d2 gas phase source'
        jls_ltop(k) = LM
        jls_power(k) = 1
        units_jls(k) = unit_string(jls_power(k),'kg/s')
c gravitational settling of SO4_d2
        k = k + 1
        jls_grav(n) = k
        sname_jls(k) = 'grav_sett_of'//trname(n)
        lname_jls(k) = 'Gravitational Settling of SO4_d2'
        jls_ltop(k) = LM
        jls_power(k) = -3
        units_jls(k) = unit_string(jls_power(k),'kg/s')

        case ('SO4_d3')
c gas phase source of SO4_d3
        k = k + 1
        jls_3Dsource(1,n) = k
        sname_jls(k) = 'gas_phase_source_of'//trname(n)
        lname_jls(k) = 'SO4_d3 gas phase source'
        jls_ltop(k) = LM
        jls_power(k) = 1
        units_jls(k) = unit_string(jls_power(k),'kg/s')
c gravitational settling of SO4_d3
        k = k + 1
        jls_grav(n) = k
        sname_jls(k) = 'grav_sett_of'//trname(n)
        lname_jls(k) = 'Gravitational Settling of SO4_d3'
        jls_ltop(k) = LM
        jls_power(k) = -3
        units_jls(k) = unit_string(jls_power(k),'kg/s')

        case ('SO4_d4')
c gas phase source of SO4_d4
        k = k + 1
        jls_3Dsource(1,n) = k
        sname_jls(k) = 'gas_phase_source_of'//trname(n)
        lname_jls(k) = 'SO4_d4 gas phase source'
        jls_ltop(k) = LM
        jls_power(k) = 1
        units_jls(k) = unit_string(jls_power(k),'kg/s')
c gravitational settling of SO4_d4
        k = k + 1
        jls_grav(n) = k
        sname_jls(k) = 'grav_sett_of'//trname(n)
        lname_jls(k) = 'Gravitational Settling of SO4_d4'
        jls_ltop(k) = LM
        jls_power(k) = -3
        units_jls(k) = unit_string(jls_power(k),'kg/s')
#endif

        case ('Be7')
c cosmogenic source from file
        k = k + 1
        jls_3Dsource(1,n) = k
        sname_jls(k) = 'Cosmogenic_src_of'//trname(n)
        lname_jls(k) = 'Be7 cosmogenic src'
        jls_ltop(k) = lm
        jls_power(k) = -28
        units_jls(k) = unit_string(jls_power(k),'kg/s/mb/m^2')
        jwt_jls(k) = 3
c radioactive decay
        k = k + 1
        jls_decay(n) = k   ! special array for all radioactive sinks
        sname_jls(k) = 'Decay_of_'//trname(n)
        lname_jls(k) = 'Loss of Be7 by decay'
        jls_ltop(k) = lm
        jls_power(k) = -28
        units_jls(k) = unit_string(jls_power(k),'kg/s/mb/m^2')
        jwt_jls(k) = 3
c gravitational settling
        k = k + 1
        jls_grav(n) = k   ! special array grav. settling sinks
        sname_jls(k) = 'Grav_Settle_of_'//trname(n)
        lname_jls(k) = 'Loss of Be7 by grav settling'
        jls_ltop(k) = lm
        jls_power(k) = -28
        units_jls(k) = unit_string(jls_power(k),'kg/s/mb/m^2')
        jwt_jls(k) = 3

        case ('Be10')
c cosmogenic source from file/same as Be7
        k = k + 1
        jls_3Dsource(1,n) = k
        sname_jls(k) = 'Cosmogenic_src_of'//trname(n)
        lname_jls(k) = 'Be10 cosmogenic src'
        jls_ltop(k) = lm
        jls_power(k) = -28  !may need changing around
        units_jls(k) = unit_string(jls_power(k),'kg/s/mb/m^2')
        jwt_jls(k) = 3
c gravitational settling
        k = k + 1
        jls_grav(n) = k   ! special array grav. settling sinks
        sname_jls(k) = 'Grav_Settle_of_'//trname(n)
        lname_jls(k) = 'Loss of Be10 by grav settling'
        jls_ltop(k) = lm
        jls_power(k) = -28  !may need changing around
        units_jls(k) = unit_string(jls_power(k),'kg/s/mb/m^2')
        jwt_jls(k) = 3

        case ('Pb210')
c source of Pb210 from Rn222 decay
        k = k + 1
        jls_3Dsource(1,n) = k
        sname_jls(k) = 'Radioactive_src_of'//trname(n)
        lname_jls(k) = 'Pb210 radioactive src'
        jls_ltop(k) = lm
        jls_power(k) =-26   ! -10  !may need to be changed
        units_jls(k) = unit_string(jls_power(k),'kg/s/mb/m^2')
        jwt_jls(k) = 3
c radioactive decay
        k = k + 1
        jls_decay(n) = k   ! special array for all radioactive sinks
        sname_jls(k) = 'Decay_of_'//trname(n)
        lname_jls(k) = 'Loss of Pb210 by decay'
        jls_ltop(k) = lm
        jls_power(k) =-26   ! -10  !may need to be changed
        units_jls(k) = unit_string(jls_power(k),'kg/s/mb/m^2')
        jwt_jls(k) = 3
c gravitational settling
        k = k + 1
        jls_grav(n) = k   ! special array grav. settling sinks
        sname_jls(k) = 'Grav_Settle_of_'//trname(n)
        lname_jls(k) = 'Loss of Pb210 by grav settling'
        jls_ltop(k) = lm
        jls_power(k) = -28
        units_jls(k) = unit_string(jls_power(k),'kg/s/mb/m^2')
        jwt_jls(k) = 3

        case ('H2O2_s')
c gas phase source and sink of H2O2
        k = k + 1
        jls_3Dsource(1,n) = k
        sname_jls(k) = 'gas_phase_source_of'//trname(n)
        lname_jls(k) = 'H2O2 gas phase source'
        jls_ltop(k) = LM
        jls_power(k) = 2
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_3Dsource(2,n) = k
        sname_jls(k) = 'gas_phase_sink_of'//trname(n)
        lname_jls(k) = 'H2O2 gas phase sink'
        jls_ltop(k) = LM
        jls_power(k) = 2
        units_jls(k) = unit_string(jls_power(k),'kg/s')
c heterogenous sink H2O2_s
        k = k + 1
        jls_3Dsource(3,n) = k
        sname_jls(k) = 'het_src_of'//trname(n)
        lname_jls(k) = 'H2O2 Heterogenous sink'
        jls_ltop(k) = LM
        jls_power(k) = -1
        units_jls(k) = unit_string(jls_power(k),'kg/s')
#ifdef TRACERS_AEROSOLS_Koch
c convective chem cloud phase sink of H2O2
        k = k + 1
        jls_incloud(1,n) = k
        sname_jls(k) = 'mc_cloud_chem_sink_of'//trname(n)
        lname_jls(k) = 'H2O2 used in convective cloud chem'
        jls_ltop(k) = LM
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),'kg/s')
c stratiform chem cloud phase sink of H2O2
        k = k + 1
        jls_incloud(2,n) = k
        sname_jls(k) = 'ss_cloud_chem_sink_of'//trname(n)
        lname_jls(k) = 'H2O2 used in stratiform cloud chem'
        jls_ltop(k) = LM
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),'kg/s')
c photolysis rate
        k = k + 1
        jls_phot = k
        sname_jls(k) = 'photolysis_rate_of_H2O2'//trname(n)
        lname_jls(k) = 'photolysis rate of H2O2'
        jls_ltop(k) =LM
        jls_power(k) =-9
        units_jls(k) = unit_string(jls_power(k),'/s')
#endif
      case ('BCII')
        k = k + 1
        jls_3Dsource(1,n) = k   ! defined but not output
        sname_jls(k) = 'unused'   ! 'Aging_sink_of'//trname(n)
        lname_jls(k) = 'unused'   ! 'BCII aging sink'
        jls_ltop(k) = LM
        jls_power(k) = -1
        units_jls(k) = unit_string(jls_power(k),'kg/s')
c gravitational settling of BCII 
        k = k + 1
        jls_grav(n) = k
        sname_jls(k) = 'grav_sett_of'//trname(n)
        lname_jls(k) = 'Gravitational Settling of BCII'
        jls_ltop(k) = LM
        jls_power(k) = -2
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_source(1,n) = k
        sname_jls(k) = 'Indust_src_'//trname(n)
        lname_jls(k) = 'BCII Industrial source'
        jls_ltop(k) = 1
        jls_power(k) = -1
        units_jls(k) = unit_string(jls_power(k),'kg/s')

      case ('BCIA')
        k = k + 1
        jls_3Dsource(1,n) = k
        sname_jls(k) = 'Aging_source_of'//trname(n)
        lname_jls(k) = 'BCIA aging source'
        jls_ltop(k) = LM
        jls_power(k) = -1
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_3Dsource(2,n) = k
        sname_jls(k) = 'Aircraft_source_of'//trname(n)
        lname_jls(k) = 'BCIA aircraft source'
        jls_ltop(k) = LM
        jls_power(k) = -1
        units_jls(k) = unit_string(jls_power(k),'kg/s')
c gravitational settling of BCIA 
        k = k + 1
        jls_grav(n) = k
        sname_jls(k) = 'grav_sett_of'//trname(n)
        lname_jls(k) = 'Gravitational Settling of BCIA'
        jls_ltop(k) = LM
        jls_power(k) = -2
        units_jls(k) = unit_string(jls_power(k),'kg/s')

      case ('BCB')
c gravitational settling of BCB 
        k = k + 1
        jls_grav(n) = k
        sname_jls(k) = 'grav_sett_of'//trname(n)
        lname_jls(k) = 'Gravitational Settling of BCB'
        jls_ltop(k) = LM
        jls_power(k) = -2
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_3Dsource(1,n) = k
        sname_jls(k) = 'biomass_src_'//trname(n)
        lname_jls(k) = 'BCB Biomass source'
        jls_ltop(k) = LM 
        jls_power(k) = -1
        units_jls(k) = unit_string(jls_power(k),'kg/s')

      case ('OCII')
c gravitational settling of OCII 
        k = k + 1
        jls_grav(n) = k
        sname_jls(k) = 'grav_sett_of'//trname(n)
        lname_jls(k) = 'Gravitational Settling of OCII'
        jls_ltop(k) = LM
        jls_power(k) = -2
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_3Dsource(1,n) = k   ! defined but not output
        sname_jls(k) = 'unused'    ! 'Aging_sink_'//trname(n)
        lname_jls(k) = 'unused'    ! 'OCIA aging sink'
        jls_ltop(k) = LM
        jls_power(k) = 1
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_source(1,n) = k
        sname_jls(k) = 'indust_src_'//trname(n)
        lname_jls(k) = 'OCIA Industrial source'
        jls_ltop(k) = 1
        jls_power(k) = -1
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_source(2,n) = k
        sname_jls(k) = 'terpene_src_'//trname(n)
        lname_jls(k) = 'OCIA Terpene source'
        jls_ltop(k) = 1
        jls_power(k) = -1
        units_jls(k) = unit_string(jls_power(k),'kg/s')

      case ('OCIA')
        k = k + 1
        jls_3Dsource(1,n) = k
        sname_jls(k) = 'Aging_source_'//trname(n)
        lname_jls(k) = 'OCIA aging source'
        jls_ltop(k) = LM
        jls_power(k) = 1
        units_jls(k) = unit_string(jls_power(k),'kg/s')
c gravitational settling of OCIA 
        k = k + 1
        jls_grav(n) = k
        sname_jls(k) = 'grav_sett_of'//trname(n)
        lname_jls(k) = 'Gravitational Settling of OCIA'
        jls_ltop(k) = LM
        jls_power(k) = -2
        units_jls(k) = unit_string(jls_power(k),'kg/s')
      case ('OCB')
c gravitational settling of OCB 
        k = k + 1
        jls_grav(n) = k
        sname_jls(k) = 'grav_sett_of'//trname(n)
        lname_jls(k) = 'Gravitational Settling of OCB'
        jls_ltop(k) = LM
        jls_power(k) = -2
        units_jls(k) = unit_string(jls_power(k),'kg/s')
c emission
        k = k + 1
        jls_3Dsource(1,n) = k
        sname_jls(k) = 'biomass_src_'//trname(n)
        lname_jls(k) = 'OCB Biomass source'
        jls_ltop(k) = LM 
        jls_power(k) = -1
        units_jls(k) = unit_string(jls_power(k),'kg/s')

      case ('seasalt1')
        k = k + 1
        jls_isrc(1,n) = k
        sname_jls(k) = 'Ocean_source_of'//trname(n)
        lname_jls(k) = 'seasalt1 ocean source'
        jls_ltop(k) = 1
        jls_power(k) = 1
        units_jls(k) = unit_string(jls_power(k),'kg/s')

c gravitational settling of ss1
        k = k + 1
        jls_grav(n) = k
        sname_jls(k) = 'grav_sett_of'//trname(n)
        lname_jls(k) = 'Gravitational Settling of seasalt1'
        jls_ltop(k) = LM
        jls_power(k) = -2
        units_jls(k) = unit_string(jls_power(k),'kg/s')

       case ('seasalt2')
        k = k + 1
        jls_isrc(1,n) = k
        sname_jls(k) = 'Ocean_source_of'//trname(n)
        lname_jls(k) = 'seasalt2 ocean source'
        jls_ltop(k) = 1
        jls_power(k) =1
        units_jls(k) = unit_string(jls_power(k),'kg/s')

c gravitational settling of ss2
        k = k + 1
        jls_grav(n) = k
        sname_jls(k) = 'grav_sett_of'//trname(n)
        lname_jls(k) = 'Gravitational Settling of seasalt2'
        jls_ltop(k) = LM
        jls_power(k) =0
        units_jls(k) = unit_string(jls_power(k),'kg/s')

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
        CASE('Clay','Silt1','Silt2','Silt3','Silt4',
     &       'ClayIlli','ClayKaol','ClaySmec','ClayCalc','ClayQuar',
     &       'Sil1Quar','Sil1Feld','Sil1Calc','Sil1Hema','Sil1Gyps',
     &       'Sil2Quar','Sil2Feld','Sil2Calc','Sil2Hema','Sil2Gyps',
     &       'Sil3Quar','Sil3Feld','Sil3Calc','Sil3Hema','Sil3Gyps',
     &       'Sil1QuHe','Sil2QuHe','Sil3QuHe')

        k=k+1
          jls_source(nDustEmjl,n)=k
          lname_jls(k)='Emission of '//TRIM(trname(n))
          sname_jls(k)=TRIM(trname(n))//'_emission'
          jls_ltop(k)=1
          jls_power(k)=1
          units_jls(k)=unit_string(jls_power(k),'kg/s')
#ifdef TRACERS_DUST
        IF (imDust == 0) THEN
          k=k+1
          jls_source(nDustEm2jl,n)=k
          lname_jls(k)='Cubic emission of '//TRIM(trname(n))
          sname_jls(k)=TRIM(trname(n))//'_emission2'
          jls_ltop(k)=1
          jls_power(k)=1
          units_jls(k)=unit_string(jls_power(k),'kg/s')
        END IF
#endif
#ifndef TRACERS_DRYDEP
        k=k+1
          jls_source(nDustTurbjl,n)=k
          lname_jls(k)='Turbulent deposition of '//TRIM(trname(n))
          sname_jls(k)=TRIM(trname(n))//'_turb_depo'
          jls_ltop(k)=1
          jls_power(k)=1
          units_jls(k)=unit_string(jls_power(k),'kg/s')
#endif
        k=k+1
          jls_grav(n)=k
          lname_jls(k)='Gain by gravitational settling of '
     &         //TRIM(trname(n))
          sname_jls(k)=TRIM(trname(n))//'_grav_sett'
          jls_ltop(k)=Lm
          jls_power(k)=1
          units_jls(k)=unit_string(jls_power(k),'kg/s')
#ifndef TRACERS_WATER
        k=k+1
          jls_wet(n)=k
          lname_jls(k)='Loss by wet deposition of '//TRIM(trname(n))
          sname_jls(k)=TRIM(trname(n))//'_wet_depo'
          jls_ltop(k)=Lm
          jls_power(k)=1
          units_jls(k)=unit_string(jls_power(k),'kg/s')
#endif
#endif

C**** Here are some more examples of generalised diag. configuration
c      n = n_dust
c        k = k + 1
c        jls_grav(n) = k   ! special array grav. settling sinks
c        sname_jls(k) = 'Grav_Settle_of_'//trname(n)
c        lname_jls(k) = 'LOSS OF DUST BY SETTLING'
c        jls_ltop(k) = lm
c        jls_power(k) = -11.
c        units_jls(k) = unit_string(jls_power(k),'kg/s')

      end select
c
C**** Checks
      if (ntsurfsrc(n).gt.ntsurfsrcmax) then
!       write(6,*) ' ntsurfsrc too large for ',trname(n)
        write(6,*) ' Increase ntsurfsrcmax to at least',ntsurfsrc(n)
        call stop_model(
     &       ' Ntsurfsrc too large.  Increase ntsurfsrcmax',255)
      end if

      end do

C**** Additional Special JL diagnostics
C**** (not necessary associated with a particular tracer)
#ifdef TRACERS_SPECIAL_Shindell
        k = k + 1
        jls_OHcon=k
        sname_jls(k) = 'OH_conc'
        lname_jls(k) = 'OH concentration'
        jls_ltop(k)  = LTOP
        jls_power(k) = 5.
        scale_jls(k) = byim
        units_jls(k) = unit_string(jls_power(k),'molecules/cm3')
c
        k = k + 1
        jls_H2Omr=k
        sname_jls(k) = 'H2O_mr'
        lname_jls(k) = 'H2O mixing ratio (weighted by daylight)'
        jls_ltop(k)  = LTOP
        jls_power(k) = -4.
        scale_jls(k) = 1.
        units_jls(k) = unit_string(jls_power(k),'parts/vol')
c
        k = k + 1
        jls_day=k
        sname_jls(k) = 'daylight'   ! not output
        lname_jls(k) = 'Daylight weighting'
        jls_ltop(k)  = 1
        jls_power(k) = 0.
        scale_jls(k) = 100.*byim
        units_jls(k) = unit_string(jls_power(k),'%')
c
        k = k + 1
        jls_N2O5sulf=k
        sname_jls(k) = 'N2O5_sulf'
        lname_jls(k) = 'N2O5 sulfate sink'
        jls_ltop(k)  = LTOP
        jls_power(k) = -2.
        units_jls(k) = unit_string(jls_power(k),'kg/s')
#ifdef regional_Ox_tracers
c
        k = k + 1
        jls_Oxloss=k
        sname_jls(k) = 'Ox_loss'
        lname_jls(k) = 'Ox chemical loss'
        jls_ltop(k)  = LTOP
        jls_power(k) = -2.
        units_jls(k) = unit_string(jls_power(k),'molecules/cm3/s')
c
        k = k + 1
        jls_Oxprod=k
        sname_jls(k) = 'Ox_prod'
        lname_jls(k) = 'Ox chemical production'
        jls_ltop(k)  = LTOP
        jls_power(k) = -2.
        units_jls(k) = unit_string(jls_power(k),'molecules/cm3/s')
#endif
#endif

#ifdef TRACERS_AEROSOLS_Koch
c Oxidants
        k = k + 1
        jls_OHconk = k
        sname_jls(k) = 'OH_conc'
        lname_jls(k) = 'OH Concentration'
        jls_ltop(k) = LM
        jls_power(k) =5
        scale_jls(k) =byim
        units_jls(k) = unit_string(jls_power(k),'molec/cm3')

        k = k + 1
        jls_HO2con = k
        sname_jls(k) = 'HO2_conc'
        lname_jls(k) = 'HO2 Concentration'
        jls_ltop(k) =LM
        jls_power(k) =7
        scale_jls(k) =byim
        units_jls(k) = unit_string(jls_power(k),'molec/cm3')

        k = k + 1
        jls_NO3 = k
        sname_jls(k) = 'NO3_conc'
        lname_jls(k) = 'NO3 Concentration'
        jls_ltop(k) =LM
        jls_power(k) =5
        scale_jls(k) =byim
        units_jls(k) = unit_string(jls_power(k),'molec/cm3')
#endif

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
      k = k + 1
      jls_spec(nDustEv1jl)=k
      lname_jls(k)='No. dust events'
      sname_jls(k)='no_dust_ev1'
      jls_ltop(k)=1
      scale_jls(k)=Sday*byim/Dtsrc
      units_jls(k)='1/d'
      k = k + 1
      jls_spec(nDustEv2jl)=k
      lname_jls(k)='No. dust events above threshold wind'
      sname_jls(k)='no_dust_ev2'
      jls_ltop(k)=1
      scale_jls(k)=Sday*byim/Dtsrc
      units_jls(k)='1/d'
      k = k + 1
      jls_spec(nDustWthjl)=k
      lname_jls(k)='Threshold velocity for dust emission'
      sname_jls(k)='wtrsh'
      jls_ltop(k)=1
      scale_jls(k)=byim
      units_jls(k)='m/s'
#endif

      if (k.gt. ktajls) then
        write (6,*)
     &   'tjl_defs: Increase ktajls=',ktajls,' to at least ',k
        call stop_model('ktajls too small',255)
      end if

C**** Tracer sources, sinks and specials
C**** Defaults for ijts (sources, sinks, etc.)
      ijts_fc(:,:)=0
      ijts_3Dsource(:,:)=0
      ijts_aq(:)=0

C**** This needs to be 'hand coded' depending on circumstances
      k = 0
      do n=1,ntm
      select case (trname(n))

      case ('SF6')
      k = k+1
        ijts_source(1,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'SF6 Layer 1 SOURCE'
        sname_ijts(k) = 'SF6_CFC-GRID_SOURCE_LAYER_1'
        ijts_power(k) = -15
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc

      case ('Rn222')
      k = k+1
        ijts_source(1,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'Rn222 L 1 SOURCE'
        sname_ijts(k) = 'Rn222_SOURCE_Layer_1'
        ijts_power(k) = -21.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc

      case ('CO2')
! keep AIJ and AJL CO2 sources in same order !!
      k = k + 1
        ijts_source(1,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        sname_ijts(k) = 'Fossil_fuel_source_'//trname(n)
        lname_ijts(k) = 'CO2 Fossil fuel src'
        ijts_power(k) = -11.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(2,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        sname_ijts(k) = 'fertilization_sink_'//trname(n)
        lname_ijts(k) = 'CO2 fertilization'
        ijts_power(k) = -11.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(3,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        sname_ijts(k) = 'Northern_forest_regrowth_'//trname(n)
        lname_ijts(k) = 'CO2 North forest regrowth'
        ijts_power(k) = -11.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(4,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        sname_ijts(k) = 'Land_Use_Modification_'//trname(n)
        lname_ijts(k) = 'CO2 from Land use mods'
        ijts_power(k) = -11.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(5,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        sname_ijts(k) = 'Ecosystem_exchange_'//trname(n)
        lname_ijts(k) = 'CO2 Ecosystem exch'
        ijts_power(k) = -11.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(6,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        sname_ijts(k) = 'Ocean_exchange_'//trname(n)
        lname_ijts(k) = 'CO2 Ocean exchange'
        ijts_power(k) = -11.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc

      case ('N2O')
#ifdef TRACERS_SPECIAL_Lerner
      k = k + 1
        ijts_source(1,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'N2O CHANGE IN L 1'
        sname_ijts(k) = 'N2O_CHANGE_IN_L_1'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
#endif
        k = k + 1
        ijts_3Dsource(1,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'N2O Chemistry'
        sname_ijts(k) = 'N2O_chem'
        ijts_power(k) = -12.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc

      case ('CFC11')
      k = k + 1
        ijts_source(1,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CFC_11 L 1 SOURCE'
        sname_ijts(k) = 'CFC_11_SOURCE_LAYER_1'
        ijts_power(k) = -15
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_3Dsource(1,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CFC_11 Stratospheric Chemistry Sink'
        sname_ijts(k) = 'CFC_11_strat_sink'
        ijts_power(k) = -18
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc

      case ('14CO2')
      k = k + 1
        ijts_source(1,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = '14CO2 L 1 Sink'
        sname_ijts(k) = '14CO2_L1_Sink'
        ijts_power(k) = -21
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc

      case ('NOx')
#ifdef EDGAR_HYDE_SOURCES
      k = k+1
        ijts_source(1,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'NOx EH fossil fuel combustion source'
        sname_ijts(k) = 'NOx_EH_ffuel_source'
        ijts_power(k) = -14
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k+1
        ijts_source(2,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'NOx EH industrial processes source'
        sname_ijts(k) = 'NOx_EH_indust_source'
        ijts_power(k) = -14
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k+1
        ijts_source(3,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'NOx EH agricultural land activities source'
        sname_ijts(k) = 'NOx_EH_Ag_land_source'
        ijts_power(k) = -14
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc  
      k = k+1
        ijts_source(4,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'NOx EH agricultural waste burning source'
        sname_ijts(k) = 'NOx_EH_Ag_waste_source'
        ijts_power(k) = -14
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc  
      k = k+1
        ijts_source(5,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'NOx EH Biofuel prod, trans, comb source'
        sname_ijts(k) = 'NOx_EH_biofuel_source'
        ijts_power(k) = -14
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k+1
        ijts_source(6,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'NOx EH deforestation source'
        sname_ijts(k) = 'NOx_EH_deforest_source'
        ijts_power(k) = -14
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k+1
        ijts_source(7,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'NOx EH savannah burning source'
        sname_ijts(k) = 'NOx_EH_savannah_source'
        ijts_power(k) = -14
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
#else
      k = k+1
        ijts_source(1,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'NOx fossil fuel source'
        sname_ijts(k) = 'NOx_fossil_fuel_source'
        ijts_power(k) = -14
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k+1
        ijts_source(2,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'NOx biomass burning source'
        sname_ijts(k) = 'NOx_biomass_burning_source'
        ijts_power(k) = -14
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k+1
        ijts_source(3,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'NOx Soil (sink)'
        sname_ijts(k) = 'NOx_Soil_(sink)'
        ijts_power(k) = -14
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
#endif
      k = k + 1
        ijts_3Dsource(nChemistry,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'NOx Chemistry'
        sname_ijts(k) = 'NOx_chem'
        ijts_power(k) = -12.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
#ifndef SHINDELL_STRAT_CHEM
      k = k + 1
        ijts_3Dsource(nStratwrite,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'NOx Stratospheric Overwrite'
        sname_ijts(k) = 'NOx_strat'
        ijts_power(k) = -12.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
#endif
      k = k + 1
        ijts_3Dsource(nLightning,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'NOx Lightning Source'
        sname_ijts(k) = 'NOx_lightning'
        ijts_power(k) = -12.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_3Dsource(nAircraft,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'NOx Aircraft Source'
        sname_ijts(k) = 'NOx_aircraft'
        ijts_power(k) = -12.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc

      case ('CO')
#ifdef EDGAR_HYDE_SOURCES
      k = k+1
        ijts_source(1,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CO EH fossil fuel combustion source'
        sname_ijts(k) = 'CO_EH_ffuel_source'
        ijts_power(k) = -13.
        units_ijts(k) = unit_string(ijts_power(k),'kg')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k+1
        ijts_source(2,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CO EH industrial processes source'
        sname_ijts(k) = 'CO_EH_indust_source'
        ijts_power(k) = -13.
        units_ijts(k) = unit_string(ijts_power(k),'kg')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc  
      k = k+1
        ijts_source(3,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CO EH agricultural waste burning source'
        sname_ijts(k) = 'CO_EH_Ag_waste_source'
        ijts_power(k) = -13.
        units_ijts(k) = unit_string(ijts_power(k),'kg')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc     
      k = k+1
        ijts_source(4,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CO EH savannah burning source'
        sname_ijts(k) = 'CO_EH_savannah_source'
        ijts_power(k) = -13.
        units_ijts(k) = unit_string(ijts_power(k),'kg')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k+1
        ijts_source(5,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CO EH biofuel prod, trans, comb source'
        sname_ijts(k) = 'CO_EH_biofuel_source'
        ijts_power(k) = -13.
        units_ijts(k) = unit_string(ijts_power(k),'kg')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k+1
        ijts_source(6,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CO EH deforestation source'
        sname_ijts(k) = 'CO_EH_deforest_source'
        ijts_power(k) = -13.
        units_ijts(k) = unit_string(ijts_power(k),'kg')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
#else
      k = k+1
        ijts_source(1,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CO industrial source'
        sname_ijts(k) = 'CO_industrial_source'
        ijts_power(k) = -13.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k+1
        ijts_source(2,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CO biomass burning source'
        sname_ijts(k) = 'CO_biomass_burning_source'
        ijts_power(k) = -13.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
#endif
      k = k + 1
        ijts_3Dsource(nChemistry,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CO Chemistry'
        sname_ijts(k) = 'CO_chem'
        ijts_power(k) = -12.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
#ifndef SHINDELL_STRAT_CHEM
      k = k + 1
        ijts_3Dsource(nStratwrite,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CO Stratospheric Overwrite'
        sname_ijts(k) = 'CO_strat'
        ijts_power(k) = -12.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
#endif

      case ('CH4')
      k = k + 1
        ijts_source(6,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 sink due to soil absorp.'
        sname_ijts(k) = 'CH4_soil_sink.'
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(7,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 Termite source'
        sname_ijts(k) = 'CH4_Termite_source'
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(9,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 Ocean source'
        sname_ijts(k) = 'CH4_Ocean_source'
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(10,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 Fresh Water lake source'
        sname_ijts(k) = 'CH4_lake_source'
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(11,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 Misc. Ground source'
        sname_ijts(k) = 'CH4_Misc._Ground_source'
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(14,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 Wetlands+Tundra source'
        sname_ijts(k) = 'CH4_Wetlands+Tundra_source'
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
#ifdef EDGAR_HYDE_SOURCES
      k = k + 1
        ijts_source(1,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 EH Animal source'
        sname_ijts(k) = 'CH4_EH_Animal_source'
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(2,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 EH Fossil Fuel Combustion source'
        sname_ijts(k) = 'CH4_EH_ffuel_comb_source'
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(3,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 EH Landfills source'
        sname_ijts(k) = 'CH4_EH_Landfill_source'
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(4,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 EH Fossil Fuel Production source'
        sname_ijts(k) = 'CH4_ffuel_prod_source'
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(5,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 EH Agricultural Waste Burning source'
        sname_ijts(k) = 'CH4_EH_Ag_Waste_source'
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(8,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 EH Deforestation source'
        sname_ijts(k) = 'CH4_EH_Deforest_source'
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(12,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 EH Savannah Burning source'
        sname_ijts(k) = 'CH4_EH_Savannah_source'
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(13,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 EH Biofuel P,T,C source'
        sname_ijts(k) = 'CH4_EH_Biofuel_source'
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(15,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 EH Agricultural Land Activities source'
        sname_ijts(k) = 'CH4_EH_Ag_Land_source'
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
#else
      k = k + 1
        ijts_source(1,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 Animal source'
        sname_ijts(k) = 'CH4_Animal_source'
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(2,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 Coal Mine source'
        sname_ijts(k) = 'CH4_Coal_Mine_source'
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(3,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 Gas Leak source'
        sname_ijts(k) = 'CH4_Gas_Leak_source'
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(4,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 Gas Venting source'
        sname_ijts(k) = 'CH4_Gas_Venting_source'
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(5,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 Municipal solid waste src'
        sname_ijts(k) = 'CH4_MSW_src'
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(8,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 Coal combustion source'
        sname_ijts(k) = 'CH4_Coal_source'
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(12,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 Biomass burning source'
        sname_ijts(k) = 'CH4_Biomass_source'
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(13,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 Rice cultivation source'
        sname_ijts(k) = 'CH4_Rice_source'
        ijts_power(k) = -13
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
#endif        
#ifdef TRACERS_SPECIAL_Shindell
      k = k + 1
        ijts_3Dsource(nChemistry,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 Chemistry'
        sname_ijts(k) = 'CH4_chem'
        ijts_power(k) = -12.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
#ifndef SHINDELL_STRAT_CHEM
      k = k + 1
        ijts_3Dsource(nStratwrite,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 Stratospheric Overwrite'
        sname_ijts(k) = 'CH4_strat'
        ijts_power(k) = -12.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
#endif
#else
      k = k + 1
        ijts_3Dsource(1,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 Tropospheric Chemistry'
        sname_ijts(k) = 'CH4_trop_chem'
        ijts_power(k) = -12.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_3Dsource(2,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 Stratospheric Chemistry'
        sname_ijts(k) = 'CH4_strat_chem'
        ijts_power(k) = -12.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
#endif

      case ('O3')
      k = k + 1
        ijts_3Dsource(1,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'O3 Tropospheric Chemistry'
        sname_ijts(k) = 'O3_trop_chem'
        ijts_power(k) = -10.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_3Dsource(2,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'O3 Stratospheric Chemistry'
        sname_ijts(k) = 'O3_strat_chem'
        ijts_power(k) = -10.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc


#ifdef TRACERS_WATER
        case ('Water')
          ! nothing I can think of....
        case ('H2O18')
        case ('HDO')
        case ('HTO')
#endif

      case ('Isoprene')
      k = k+1
        ijts_source(1,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'Isoprene vegetation source'
        sname_ijts(k) = 'Isoprene_vegetation_source'
        ijts_power(k) = -13.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_3Dsource(nChemistry,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'Isoprene Chemistry'
        sname_ijts(k) = 'Isoprene_chem'
        ijts_power(k) = -12.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
#ifndef SHINDELL_STRAT_CHEM
      k = k + 1
        ijts_3Dsource(nStratwrite,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'Isoprene Stratospheric Overwrite'
        sname_ijts(k) = 'Isoprene_strat'
        ijts_power(k) = -12.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
#endif

      case ('Alkenes','Paraffin')
#ifdef EDGAR_HYDE_SOURCES
      k = k+1
        ijts_source(1,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = trname(n)//' EH fossil fuel combustion source'
        sname_ijts(k) = trname(n)//'_EH_ffuel_comb_source'
        ijts_power(k) = -13.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k+1
        ijts_source(2,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = trname(n)//' EH industrial processes source'
        sname_ijts(k) = trname(n)//'_EH_indust_source'
        ijts_power(k) = -13.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
!! #3 is the (common) vegetation source. See below.
      k = k+1
        ijts_source(4,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = trname(n)//' EH fossil fuel production source'
        sname_ijts(k) = trname(n)//'_EH_ffuel_prod_source'
        ijts_power(k) = -13.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k+1
        ijts_source(5,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 
     &  trname(n)//' EH agricultural waste burning source'
        sname_ijts(k) = trname(n)//'_EH_Ag_waste_source'
        ijts_power(k) = -13.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k+1
        ijts_source(6,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = trname(n)//' EH savannah burning source'
        sname_ijts(k) = trname(n)//'_EH_savannah_source'
        ijts_power(k) = -13.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k+1
        ijts_source(7,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = trname(n)//' EH biofuel P,T,C source'
        sname_ijts(k) = trname(n)//'_EH_biofuel_source'
        ijts_power(k) = -13.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k+1
        ijts_source(8,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = trname(n)//' EH deforestation source'
        sname_ijts(k) = trname(n)//'_EH_deforest_source'
        ijts_power(k) = -13.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
#else
      k = k+1
        ijts_source(1,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = trname(n)//' industrial source'
        sname_ijts(k) = trname(n)//'_industrial_source'
        ijts_power(k) = -13.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k+1
        ijts_source(2,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = trname(n)//' biomass burning source'
        sname_ijts(k) = trname(n)//'_biomass_burning_source'
        ijts_power(k) = -13.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
!! #3 is the (common) vegetation source. See below.
#endif
      k = k+1
        ijts_source(3,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = trname(n)//' vegetation source'
        sname_ijts(k) = trname(n)//'_vegetation_source'
        ijts_power(k) = -13.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_3Dsource(nChemistry,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = trname(n)//' Chemistry'
        sname_ijts(k) = trname(n)//'_chem'
        ijts_power(k) = -12.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
#ifndef SHINDELL_STRAT_CHEM
      k = k + 1
        ijts_3Dsource(nStratwrite,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = trname(n)//' Stratospheric Overwrite'
        sname_ijts(k) = trname(n)//'_strat'
        ijts_power(k) = -12.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
#endif

      case ('ClOx','BrOx','HCl','HOCl','ClONO2',
     &      'HBr','HOBr','BrONO2','CFC')
      k = k + 1
        ijts_3Dsource(nChemistry,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = trname(n)//' Chemistry'
        sname_ijts(k) = trname(n)//'_chem'
        ijts_power(k) = -12.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc

      case ('H2O2','CH3OOH','Ox','N2O5','HNO3','HCHO',
     &      'HO2NO2','PAN','AlkylNit','OxREG1','OxREG2',
     &      'OxREG3','OxREG4','OxREG5','OxREG6')
      k = k + 1
        ijts_3Dsource(nChemistry,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = trname(n)//' Chemistry'
        sname_ijts(k) = trname(n)//'_chem'
        ijts_power(k) = -12.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
#ifndef SHINDELL_STRAT_CHEM
      k = k + 1
        ijts_3Dsource(nStratwrite,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = trname(n)//' Stratospheric Overwrite'
        sname_ijts(k) = trname(n)//'_strat'
        ijts_power(k) = -12.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
#endif
        select case(trname(n))
        case('Ox')
        k = k + 1
          ijts_fc(1,n) = k
          ijts_index(k) = n
          ia_ijts(k) = ia_rad
          lname_ijts(k) = trname(n)//' SW radiative forcing'
          sname_ijts(k) = 'swf_'//trim(trname(n))
          ijts_power(k) = -2.
          units_ijts(k) = unit_string(ijts_power(k),'W/m2')
          scale_ijts(k) = 10.**(-ijts_power(k))
        k = k + 1
          ijts_fc(2,n) = k
          ijts_index(k) = n
          ia_ijts(k) = ia_rad
          lname_ijts(k) = trname(n)//' LW radiative forcing'
          sname_ijts(k) = 'lwf_'//trim(trname(n))
          ijts_power(k) = -2.
          units_ijts(k) = unit_string(ijts_power(k),'W/m2')
          scale_ijts(k) = 10.**(-ijts_power(k))
        end select

      case ('BCII')
        k = k + 1
        ijts_source(1,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'BC Industrial source'
        sname_ijts(k) = 'BC_Industrial_source'
        ijts_power(k) = -12.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
        k = k + 1
        ijts_3Dsource(1,n) = k   ! defined but not output
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'unused'    ! 'BCII Aging sink'
        sname_ijts(k) = 'unused'    ! 'BCII_Aging_sink'
        ijts_power(k) = -12.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc

      case ('BCIA')
        k = k + 1
        ijts_3Dsource(1,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'BC Aging source'
        sname_ijts(k) = 'BC_Aging_Source'
        ijts_power(k) = -12.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
        k = k + 1
        ijts_3Dsource(2,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'BCIA Aircraft source'
        sname_ijts(k) = 'BCIA_Aircraft_src'
        ijts_power(k) = -12.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
#ifdef TRACERS_AEROSOLS_Koch
c BCI optical thickness
        k = k + 1
        ijts_tau(1,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_rad
        lname_ijts(k) = 'BCI optical thickness'
        sname_ijts(k) = 'tau_BCI'
        ijts_power(k) = -2.
        units_ijts(k) = unit_string(ijts_power(k),' ')
        scale_ijts(k) = 10.**(-ijts_power(k))
c BCI clear sky optical thickness
        k = k + 1
        ijts_tau(2,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_rad
        lname_ijts(k) = 'BCI clr sky optical thickness'
        sname_ijts(k) = 'tau_CS_BCI'
        ijts_power(k) = -2.
        units_ijts(k) = unit_string(ijts_power(k),' ')
        scale_ijts(k) = 10.**(-ijts_power(k))
c BCI shortwave radiative forcing
        k = k + 1
        ijts_fc(1,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_rad
        lname_ijts(k) = 'BCI SW radiative forcing'
        sname_ijts(k) = 'swf_BCI'
        ijts_power(k) = -2.
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
c BCI longwave radiative forcing
        k = k + 1
        ijts_fc(2,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_rad
        lname_ijts(k) = 'BCI LW radiative forcing'
        sname_ijts(k) = 'lwf_BCI'
        ijts_power(k) = -2.
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
c BCI clear sky shortwave radiative forcing
        k = k + 1
        ijts_fc(5,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_rad
        lname_ijts(k) = 'BCI clr sky SW rad forcing'
        sname_ijts(k) = 'swf_CS_BCI'
        ijts_power(k) = -2.
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
c BCI longwave radiative forcing
        k = k + 1
        ijts_fc(6,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_rad
        lname_ijts(k) = 'BCI clr sky LW rad forcing'
        sname_ijts(k) = 'lwf_CS_BCI'
        ijts_power(k) = -2.
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
c BCI clear sky shortwave radiative forcing
        k = k + 1
        ijts_fc(5,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_rad
        lname_ijts(k) = 'BCI clr sky SW rad forcing'
        sname_ijts(k) = 'swrf_CS_BCI'
        ijts_power(k) = -2.
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
c BCI longwave radiative forcing
        k = k + 1
        ijts_fc(6,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_rad
        lname_ijts(k) = 'BCI clr sky LW rad forcing'
        sname_ijts(k) = 'lwrf_CS_BCI'
        ijts_power(k) = -2.
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
#endif
      case ('BCB')
        k = k + 1
        ijts_3Dsource(1,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'BC Biomass source'
        sname_ijts(k) = 'BC_Biomass_source'
        ijts_power(k) = -12.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
#ifdef TRACERS_AEROSOLS_Koch
c BCB optical thickness
        k = k + 1
        ijts_tau(1,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_rad
        lname_ijts(k) = 'BCB optical thickness'
        sname_ijts(k) = 'tau_BCB'
        ijts_power(k) = -2.
        units_ijts(k) = unit_string(ijts_power(k),' ')
        scale_ijts(k) = 10.**(-ijts_power(k))
c BCB clear sky optical thickness
        k = k + 1
        ijts_tau(2,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_rad
        lname_ijts(k) = 'BCB clr sky optical thickness'
        sname_ijts(k) = 'tau_CS_BCB'
        ijts_power(k) = -2.
        units_ijts(k) = unit_string(ijts_power(k),' ')
        scale_ijts(k) = 10.**(-ijts_power(k))
c BCB shortwave radiative forcing
        k = k + 1
        ijts_fc(1,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_rad
        lname_ijts(k) = 'BCB SW radiative forcing'
        sname_ijts(k) = 'swf_BCB'
        ijts_power(k) = -2.
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
c BCB longwave radiative forcing
        k = k + 1
        ijts_fc(2,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_rad
        lname_ijts(k) = 'BCB LW radiative forcing'
        sname_ijts(k) = 'lwf_BCB'
        ijts_power(k) = -2.
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
c BCB clear sky shortwave radiative forcing
        k = k + 1
        ijts_fc(5,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_rad
        lname_ijts(k) = 'BCB clr sky SW rad forcing'
        sname_ijts(k) = 'swf_CS_BCB'
        ijts_power(k) = -2.
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
c BCB clear sky longwave radiative forcing
        k = k + 1
        ijts_fc(6,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_rad
        lname_ijts(k) = 'BCB clr sky LW rad forcing'
        sname_ijts(k) = 'lwf_CS_BCB'
        ijts_power(k) = -2.
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
c BCB clear sky shortwave radiative forcing
        k = k + 1
        ijts_fc(5,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_rad
        lname_ijts(k) = 'BCB clr sky SW rad forcing'
        sname_ijts(k) = 'swrf_CS_BCB'
        ijts_power(k) = -2.
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
c BCB clear sky longwave radiative forcing
        k = k + 1
        ijts_fc(6,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_rad
        lname_ijts(k) = 'BCB clr sky LW rad forcing'
        sname_ijts(k) = 'lwrf_CS_BCB'
        ijts_power(k) = -2.
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
#endif

      case ('OCII')
        k = k + 1
        ijts_source(1,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'OC Industrial source'
        sname_ijts(k) = 'OC_Industrial_source'
        ijts_power(k) = -12.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
        k = k + 1
        ijts_source(2,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'Terpene source'
        sname_ijts(k) = 'OC_Terpene_source'
        ijts_power(k) = -12.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
        k = k + 1
        ijts_3Dsource(1,n) = k   ! defined but not output
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'unused'       ! 'OCII Aging sink'
        sname_ijts(k) = 'unused'       ! 'OCII_Aging_Sink'
        ijts_power(k) = -12.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc

      case ('OCIA')
        k = k + 1
        ijts_3Dsource(1,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'OC Aging source'
        sname_ijts(k) = 'OC_Aging_Source'
        ijts_power(k) = -12.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
#ifdef TRACERS_AEROSOLS_Koch
c OC optical thickness
        k = k + 1
        ijts_tau(1,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_rad
        lname_ijts(k) = 'OC optical thickness'
        sname_ijts(k) = 'tau_OC'
        ijts_power(k) = -2.
        units_ijts(k) = unit_string(ijts_power(k),' ')
        scale_ijts(k) = 10.**(-ijts_power(k))
c OC clear sky optical thickness
        k = k + 1
        ijts_tau(2,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_rad
        lname_ijts(k) = 'OC clr sky optical thickness'
        sname_ijts(k) = 'tau_CS_OC'
        ijts_power(k) = -2.
        units_ijts(k) = unit_string(ijts_power(k),' ')
        scale_ijts(k) = 10.**(-ijts_power(k))
c OC shortwave radiative forcing
        k = k + 1
        ijts_fc(1,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_rad
        lname_ijts(k) = 'OC SW radiative forcing'
        sname_ijts(k) = 'swf_OC'
        ijts_power(k) = -2.
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
c OC longwave radiative forcing
        k = k + 1
        ijts_fc(2,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_rad
        lname_ijts(k) = 'OC LW radiative forcing'
        sname_ijts(k) = 'lwf_OC'
        ijts_power(k) = -2.
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
c OC clear sky shortwave radiative forcing
        k = k + 1
        ijts_fc(5,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_rad
        lname_ijts(k) = 'OC clr sky SW rad forcing'
        sname_ijts(k) = 'swf_CS_OC'
        ijts_power(k) = -2.
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
c OC clear sky longwave radiative forcing
        k = k + 1
        ijts_fc(6,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_rad
        lname_ijts(k) = 'OC clr sky LW rad forcing'
        sname_ijts(k) = 'lwf_CS_OC'
        ijts_power(k) = -2.
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
c OC clear sky shortwave radiative forcing
        k = k + 1
        ijts_fc(5,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_rad
        lname_ijts(k) = 'OC clr sky SW rad forcing'
        sname_ijts(k) = 'swrf_CS_OC'
        ijts_power(k) = -2.
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
c OC clear sky longwave radiative forcing
        k = k + 1
        ijts_fc(6,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_rad
        lname_ijts(k) = 'OC clr sky LW rad forcing'
        sname_ijts(k) = 'lwrf_CS_OC'
        ijts_power(k) = -2.
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
#endif
      case ('OCB')
        k = k + 1
        ijts_3Dsource(1,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'OC Biomass source'
        sname_ijts(k) = 'OC_Biomass_source'
        ijts_power(k) = -12.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc

      case ('DMS')
        k = k + 1
        ijts_isrc(1,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'DMS Ocean source'
        sname_ijts(k) = 'DMS_Ocean_source'
        ijts_power(k) = -12.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc

        k = k + 1
        ijts_3Dsource(1,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'DMS Chem sink'
        sname_ijts(k) = 'DMS_Chem_sink'
        ijts_power(k) = -12.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc

      case ('MSA')
c put in chemical production of MSA
        k = k + 1
        ijts_3Dsource(1,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'MSA Chemical source'
        sname_ijts(k) = 'MSA_Chemical_source'
        ijts_power(k) = -17.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc

      case ('SO2')
c production of SO2 from volcanic emissions
        k = k + 1
        ijts_3Dsource(1,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'SO2 source from volcanos'
        sname_ijts(k) = 'SO2_source_from_volcanos'
        ijts_power(k) = -15.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
c production of SO2 from aircraft
        k = k + 1
        ijts_3Dsource(2,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'SO2 source from aircraft'
        sname_ijts(k) = 'SO2_source_from_aricraft'
        ijts_power(k) = -15.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
c emissions of biomass SO2
        k = k + 1
        ijts_3Dsource(3,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'Biomass SO2 source'
        sname_ijts(k) = 'SO2_source_from_biomass'
        ijts_power(k) = -15.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
c put in production of SO2 from DMS
        k = k + 1
        ijts_3Dsource(4,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'SO2 source from DMS'
        sname_ijts(k) = 'SO2_source_from_DMS'
        ijts_power(k) = -15.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
c put in chemical loss of SO2
        k = k + 1
        ijts_3Dsource(5,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'SO2 Chemical sink'
        sname_ijts(k) = 'SO2_chem_sink'
        ijts_power(k) = -15.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
c emissions of industrial SO2
        k = k + 1
        ijts_source(1,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'Industrial SO2 source'
        sname_ijts(k) = 'SO2_source_from_industry'
        ijts_power(k) = -15.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
c put in production of SO2 from heter chem
        k = k + 1
        ijts_3Dsource(6,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'SO2 heter chem sink'
        sname_ijts(k) = 'SO2_het_chem_sink'
        ijts_power(k) = -15.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc

      case ('SO4')
c put in production of SO4 from gas phase
        k = k + 1
        ijts_3Dsource(1,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'SO4 gas phase source'
        sname_ijts(k) = 'SO4_gas_phase_source'
        ijts_power(k) = -15.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
c SO4 from industrial emissions
        k = k + 1
        ijts_source(1,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'Industrial SO4 source'
        sname_ijts(k) = 'SO4_source_from_industry'
        ijts_power(k) = -15.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
c put in loss of SO4 from heter chem
        k = k + 1
        ijts_3Dsource(2,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'SO4 heter chem sink'
        sname_ijts(k) = 'SO4_het_chem_sink'
        ijts_power(k) = -10.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
#ifdef TRACERS_AEROSOLS_Koch
c put in source of SO4 from aqueous chem
        k = k + 1
        ijts_aq(n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'SO4 aqueous chem source'
        sname_ijts(k) = 'SO4_aq_chem_source'
        ijts_power(k) = -15.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
c SO4 optical thickness
        k = k + 1
        ijts_tau(1,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_rad
        lname_ijts(k) = 'SO4 optical thickness'
        sname_ijts(k) = 'tau_'//trim(trname(n))
        ijts_power(k) = -2.
        units_ijts(k) = unit_string(ijts_power(k),' ')
        scale_ijts(k) = 10.**(-ijts_power(k))
c SO4 clear sky optical thickness
        k = k + 1
        ijts_tau(2,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_rad
        lname_ijts(k) = 'SO4 clr sky optical thickness'
        sname_ijts(k) = 'tau_CS_'//trim(trname(n))
        ijts_power(k) = -2.
        units_ijts(k) = unit_string(ijts_power(k),' ')
        scale_ijts(k) = 10.**(-ijts_power(k))
c SO4 shortwave radiative forcing
        k = k + 1
        ijts_fc(1,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_rad
        lname_ijts(k) = 'SO4 SW radiative forcing'
        sname_ijts(k) = 'swf_'//trim(trname(n))
        ijts_power(k) = -2.
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
c SO4 longwave radiative forcing
        k = k + 1
        ijts_fc(2,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_rad
        lname_ijts(k) = 'SO4 LW radiative forcing'
        sname_ijts(k) = 'lwf_'//trim(trname(n))
        ijts_power(k) = -2.
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
c SO4 clear sky shortwave radiative forcing
        k = k + 1
        ijts_fc(5,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_rad
        lname_ijts(k) = 'SO4 clr sky SW rad forcing'
        sname_ijts(k) = 'swf_CS'//trim(trname(n))
        ijts_power(k) = -2.
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
c SO4 clear sky longwave radiative forcing
        k = k + 1
        ijts_fc(6,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_rad
        lname_ijts(k) = 'SO4 clr sky LW rad forcing'
        sname_ijts(k) = 'lwf_CS'//trim(trname(n))
        ijts_power(k) = -2.
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
c SO4 clear sky shortwave radiative forcing
        k = k + 1
        ijts_fc(5,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_rad
        lname_ijts(k) = 'SO4 clr sky SW rad forcing'
        sname_ijts(k) = 'swrf_CS'//trim(trname(n))
        ijts_power(k) = -2.
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
c SO4 clear sky longwave radiative forcing
        k = k + 1
        ijts_fc(6,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_rad
        lname_ijts(k) = 'SO4 clr sky LW rad forcing'
        sname_ijts(k) = 'lwrf_CS'//trim(trname(n))
        ijts_power(k) = -2.
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
#endif
#ifdef TRACERS_HETCHEM
      case ('SO4_d1')
c chemical production of SO4 from SO2 on dust
        k = k + 1
        ijts_source(1,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'SO4d1 Chemical source'
        sname_ijts(k) = 'SO4d1_Chemical_source'
        ijts_power(k) = -10.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
c SO4 dry dep
        k = k + 1
        ijts_source(2,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'SO4d1 Dry Dep'
        sname_ijts(k) = 'SO4d1_Dry_Dep'
        ijts_power(k) = -10.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
c SO4 wet dep
        k = k + 1
        ijts_source(3,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'SO4d1 Wet Dep'
        sname_ijts(k) = 'SO4d1_Wet_Dep'
        ijts_power(k) = -10.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
c SO4 settling
        k = k + 1
        ijts_source(4,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'Gravitational settling of '//trname(n)
        sname_ijts(k) = 'Grav_Settle_of_'//trname(n)
        ijts_power(k) = -10.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      case ('SO4_d2')
c chemical production of SO4 from SO2 on dust
        k = k + 1
        ijts_source(1,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'SO4d2 Chemical source'
        sname_ijts(k) = 'SO4d2_Chemical_source'
        ijts_power(k) = -10.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
c SO4 dry dep
        k = k + 1
        ijts_source(2,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'SO4d2 Dry Dep'
        sname_ijts(k) = 'SO4d2_Dry_Dep'
        ijts_power(k) = -10.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      case ('SO4_d3')
c chemical production of SO4 from SO2 on dust
        k = k + 1
        ijts_source(1,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'SO4d3 Chemical source'
        sname_ijts(k) = 'SO4d3_Chemical_source'
        ijts_power(k) = -10.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
c SO4 dry dep
        k = k + 1
        ijts_source(2,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'SO4d3 Dry Dep'
        sname_ijts(k) = 'SO4d3_Dry_Dep'
        ijts_power(k) = -10.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      case ('SO4_d4')
c chemical production of SO4 from SO2 on dust
        k = k + 1
        ijts_source(1,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'SO4d4 Chemical source'
        sname_ijts(k) = 'SO4d4_Chemical_source'
        ijts_power(k) = -10.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
c SO4 dry dep
        k = k + 1
        ijts_source(2,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'SO4d4 Dry Dep'
        sname_ijts(k) = 'SO4d4_Dry_Dep'
        ijts_power(k) = -10.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
#endif
      case ('H2O2_s')
c put in production of H2O2 from gas phase
        k = k + 1
        ijts_3Dsource(1,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'H2O2 gas phase source'
        sname_ijts(k) = 'H2O2_gas_phase_source'
        ijts_power(k) = -10.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
c put in production of H2O2 from gas phase
        k = k + 1
        ijts_3Dsource(2,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'H2O2 gas phase sink'
        sname_ijts(k) = 'H2O2_gas_phase_sink'
        ijts_power(k) = -10.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
c put in production of H2O2 from heter chem
        k = k + 1
        ijts_3Dsource(3,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'H2O2 heter chem sink'
        sname_ijts(k) = 'H2O2_het_chem_sink'
        ijts_power(k) = -10.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc

        case ('Be7')
c cosmogenic source from file
        k = k + 1
        ijts_3Dsource(1,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'Cosmogenic source of '//trname(n)
        sname_ijts(k) = 'Be7_cosmo_src'
        ijts_power(k) = -25
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc

        case ('Be10')
c cosmogenic source from file
        k = k + 1
        ijts_3Dsource(1,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'Cosmogenic source of '//trname(n)
        sname_ijts(k) = 'Be10_cosmo_src'
        ijts_power(k) = -25
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc

        case ('Pb210')
c source of Pb210 from Rn222 decay
        k = k + 1
        ijts_3Dsource(1,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'Radioactive source of '//trname(n)
        sname_ijts(k) = 'Pb210_radio_src'
        ijts_power(k) = -24
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc

      case ('seasalt1')
        k = k + 1
        ijts_isrc(1,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'seasalt1 Ocean source'
        sname_ijts(k) = 'seasalt1_Ocean_source'
        ijts_power(k) = -12.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
#ifdef TRACERS_AEROSOLS_Koch
c ss1 optical thickness
        k = k + 1
        ijts_tau(1,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_rad
        lname_ijts(k) = 'ss1 optical thickness'
        sname_ijts(k) = 'tau_'//trim(trname(n))
        ijts_power(k) = -2.
        units_ijts(k) = unit_string(ijts_power(k),' ')
        scale_ijts(k) = 10.**(-ijts_power(k))
c ss1 clear sky optical thickness
        k = k + 1
        ijts_tau(2,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_rad
        lname_ijts(k) = 'ss1 clr sky optical thickness'
        sname_ijts(k) = 'tau_CS_'//trim(trname(n))
        ijts_power(k) = -2.
        units_ijts(k) = unit_string(ijts_power(k),' ')
        scale_ijts(k) = 10.**(-ijts_power(k))
c SS shortwave radiative forcing
        k = k + 1
        ijts_fc(1,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_rad
        lname_ijts(k) = 'SS SW radiative forcing'
        sname_ijts(k) = 'swf_SS'
        ijts_power(k) = -2.
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
c SS longwave radiative forcing
        k = k + 1
        ijts_fc(2,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_rad
        lname_ijts(k) = 'SS LW radiative forcing'
        sname_ijts(k) = 'lwf_SS'
        ijts_power(k) = -2.
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
c SS clear sky shortwave radiative forcing
        k = k + 1
        ijts_fc(5,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_rad
        lname_ijts(k) = 'SS clr sky SW rad forcing'
        sname_ijts(k) = 'swf_CS_SS'
        ijts_power(k) = -2.
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
c SS clear sky longwave radiative forcing
        k = k + 1
        ijts_fc(6,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_rad
        lname_ijts(k) = 'SS clr sky LW rad forcing'
        sname_ijts(k) = 'lwf_CS_SS'
        ijts_power(k) = -2.
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
c SS clear sky shortwave radiative forcing
        k = k + 1
        ijts_fc(5,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_rad
        lname_ijts(k) = 'SS clr sky SW rad forcing'
        sname_ijts(k) = 'swrf_CS_SS'
        ijts_power(k) = -2.
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
c SS clear sky longwave radiative forcing
        k = k + 1
        ijts_fc(6,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_rad
        lname_ijts(k) = 'SS clr sky LW rad forcing'
        sname_ijts(k) = 'lwrf_CS_SS'
        ijts_power(k) = -2.
        units_ijts(k) = unit_string(ijts_power(k),'W/m2')
        scale_ijts(k) = 10.**(-ijts_power(k))
#endif
       case ('seasalt2')
        k = k + 1
        ijts_isrc(1,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'seasalt2 Ocean source'
        sname_ijts(k) = 'seasalt2_Ocean_source'
        ijts_power(k) = -12.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
#ifdef TRACERS_AEROSOLS_Koch
c ss2 optical thickness
        k = k + 1
        ijts_tau(1,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_rad
        lname_ijts(k) = 'ss2 optical thickness'
        sname_ijts(k) = 'tau_'//trim(trname(n))
        ijts_power(k) = -2.
        units_ijts(k) = unit_string(ijts_power(k),' ')
        scale_ijts(k) = 10.**(-ijts_power(k))
c ss2 clear sky optical thickness
        k = k + 1
        ijts_tau(2,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_rad
        lname_ijts(k) = 'ss2 clr sky optical thickness'
        sname_ijts(k) = 'tau_CS_'//trim(trname(n))
        ijts_power(k) = -2.
        units_ijts(k) = unit_string(ijts_power(k),' ')
        scale_ijts(k) = 10.**(-ijts_power(k))
C**** Additional Special IJ diagnostics 
      do L=1,LTOP
        k = k + 1
          ijs_dms_dens(L)=k
          ijts_index(k) = ntm
          ia_ijts(k) = ia_src
          write(lname_ijts(k),'(a20,i2.2)') 'DMS concentration L=',L
          write(sname_ijts(k),'(a8,i2.2)') 'DMS_con_',L
          ijts_power(k) = 6.
          units_ijts(k) = unit_string(ijts_power(k),'molecules/cm3')
          scale_ijts(k) = 10.**(-ijts_power(k))
      end do
      do L=1,LTOP
       k = k + 1
          ijs_so2_dens(L)=k
          ijts_index(k) = ntm
          ia_ijts(k) = ia_src
          write(lname_ijts(k),'(a20,i2.2)') 'SO2 concentration L=',L
          write(sname_ijts(k),'(a8,i2.2)') 'SO2_con_',L
          ijts_power(k) = 6.
          units_ijts(k) = unit_string(ijts_power(k),'molecules/cm3')
          scale_ijts(k) = 10.**(-ijts_power(k))
      end do
      do L=1,LTOP
        k = k + 1
          ijs_so4_dens(L)=k
          ijts_index(k) = ntm
          ia_ijts(k) = ia_src
          write(lname_ijts(k),'(a20,i2.2)') 'SO4 concentration L=',L
          write(sname_ijts(k),'(a8,i2.2)') 'SO4_con_',L
          ijts_power(k) = 6.
          units_ijts(k) = unit_string(ijts_power(k),'molecules/cm3')
          scale_ijts(k) = 10.**(-ijts_power(k))
      end do
#endif

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
      CASE('Clay','Silt1','Silt2','Silt3','Silt4',
     &   'ClayIlli','ClayKaol','ClaySmec','ClayCalc','ClayQuar',
     &   'Sil1Quar','Sil1Feld','Sil1Calc','Sil1Hema','Sil1Gyps',
     &   'Sil2Quar','Sil2Feld','Sil2Calc','Sil2Hema','Sil2Gyps',
     &   'Sil3Quar','Sil3Feld','Sil3Calc','Sil3Hema','Sil3Gyps',
     &   'Sil1QuHe','Sil2QuHe','Sil3QuHe')
        k=k+1
        ijts_source(nDustEmij,n)=k
        lname_ijts(k)='Emission of '//TRIM(trname(n))
        sname_ijts(k)=TRIM(trname(n))//'_emission'
        ijts_index(k)=n
        ia_ijts(k)=ia_src
        ijts_power(k) = -13.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
#ifdef TRACERS_DUST
        IF (imDust == 0) THEN
        k=k+1
        ijts_source(nDustEm2ij,n)=k
        lname_ijts(k)='Cubic emission of '//TRIM(trname(n))
        sname_ijts(k)=TRIM(trname(n))//'_emission2'
        ijts_index(k)=n
        ia_ijts(k)=ia_src
        ijts_power(k) = -13.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
        END IF
#endif
#ifndef TRACERS_DRYDEP
      k=k+1
        ijts_source(nDustTurbij,n)=k
        lname_ijts(k)='Turbulent Deposition of '//TRIM(trname(n))
        sname_ijts(k)=TRIM(trname(n))//'_turb_depo'
        ijts_index(k)=n
        ia_ijts(k)=ia_src
        ijts_power(k) = -13.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
#endif
#ifndef TRACERS_WATER
      k=k+1
        ijts_wet(n)=k
        lname_ijts(k)='Wet deposition of '//TRIM(trname(n))
        sname_ijts(k)=TRIM(trname(n))//'_wet_depo'
        ijts_index(k)=n
        ia_ijts(k)=ia_src
        ijts_power(k) = -13.
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
#endif
        SELECT CASE (trname(n))
        CASE ('Clay')
c dust optical thickness of four clay sub size classes
          k = k + 1
          ijts_tausub(1,n,1) = k
          ijts_index(k) = n
          ia_ijts(k) = ia_rad
          lname_ijts(k) = trim(trname(n))//'1 optical thickness'
          sname_ijts(k) = 'tau_'//trim(trname(n))//'1'
          ijts_power(k) = -2.
          units_ijts(k) = unit_string(ijts_power(k),' ')
          scale_ijts(k) = 10.**(-ijts_power(k))
          k = k + 1
          ijts_tausub(1,n,2) = k
          ijts_index(k) = n
          ia_ijts(k) = ia_rad
          lname_ijts(k) = trim(trname(n))//'2 optical thickness'
          sname_ijts(k) = 'tau_'//trim(trname(n))//'2'
          ijts_power(k) = -2.
          units_ijts(k) = unit_string(ijts_power(k),' ')
          scale_ijts(k) = 10.**(-ijts_power(k))
          k = k + 1
          ijts_tausub(1,n,3) = k
          ijts_index(k) = n
          ia_ijts(k) = ia_rad
          lname_ijts(k) = trim(trname(n))//'3 optical thickness'
          sname_ijts(k) = 'tau_'//trim(trname(n))//'3'
          ijts_power(k) = -2.
          units_ijts(k) = unit_string(ijts_power(k),' ')
          scale_ijts(k) = 10.**(-ijts_power(k))
          k = k + 1
          ijts_tausub(1,n,4) = k
          ijts_index(k) = n
          ia_ijts(k) = ia_rad
          lname_ijts(k) = trim(trname(n))//'4 optical thickness'
          sname_ijts(k) = 'tau_'//trim(trname(n))//'4'
          ijts_power(k) = -2.
          units_ijts(k) = unit_string(ijts_power(k),' ')
          scale_ijts(k) = 10.**(-ijts_power(k))
c dust clear sky optical thickness of four clay sub size classes
          k = k + 1
          ijts_tausub(2,n,1) = k
          ijts_index(k) = n
          ia_ijts(k) = ia_rad
          lname_ijts(k) = trim(trname(n))//'1 CS optical thickness'
          sname_ijts(k) = 'tau_CS_'//trim(trname(n))//'1'
          ijts_power(k) = -2.
          units_ijts(k) = unit_string(ijts_power(k),' ')
          scale_ijts(k) = 10.**(-ijts_power(k))
          k = k + 1
          ijts_tausub(2,n,2) = k
          ijts_index(k) = n
          ia_ijts(k) = ia_rad
          lname_ijts(k) = trim(trname(n))//'2 CS optical thickness'
          sname_ijts(k) = 'tau_CS_'//trim(trname(n))//'2'
          ijts_power(k) = -2.
          units_ijts(k) = unit_string(ijts_power(k),' ')
          scale_ijts(k) = 10.**(-ijts_power(k))
          k = k + 1
          ijts_tausub(2,n,3) = k
          ijts_index(k) = n
          ia_ijts(k) = ia_rad
          lname_ijts(k) = trim(trname(n))//'3 CS optical thickness'
          sname_ijts(k) = 'tau_CS_'//trim(trname(n))//'3'
          ijts_power(k) = -2.
          units_ijts(k) = unit_string(ijts_power(k),' ')
          scale_ijts(k) = 10.**(-ijts_power(k))
          k = k + 1
          ijts_tausub(2,n,4) = k
          ijts_index(k) = n
          ia_ijts(k) = ia_rad
          lname_ijts(k) = trim(trname(n))//'4 CS optical thickness'
          sname_ijts(k) = 'tau_CS_'//trim(trname(n))//'4'
          ijts_power(k) = -2.
          units_ijts(k) = unit_string(ijts_power(k),' ')
          scale_ijts(k) = 10.**(-ijts_power(k))
c dust shortwave radiative forcing of four clay sub size classes
          k = k + 1
          ijts_fcsub(1,n,1) = k
          ijts_index(k) = n
          ia_ijts(k) = ia_rad
          lname_ijts(k) = trim(trname(n))//'1 SW radiative forcing'
          sname_ijts(k) = 'swf_'//trim(trname(n))//'1'
          ijts_power(k) = -2.
          units_ijts(k) = unit_string(ijts_power(k),'W/m2')
          scale_ijts(k) = 10.**(-ijts_power(k))
          k = k + 1
          ijts_fcsub(1,n,2) = k
          ijts_index(k) = n
          ia_ijts(k) = ia_rad
          lname_ijts(k) = trim(trname(n))//'2 SW radiative forcing'
          sname_ijts(k) = 'swf_'//trim(trname(n))//'2'
          ijts_power(k) = -2.
          units_ijts(k) = unit_string(ijts_power(k),'W/m2')
          scale_ijts(k) = 10.**(-ijts_power(k))
          k = k + 1
          ijts_fcsub(1,n,3) = k
          ijts_index(k) = n
          ia_ijts(k) = ia_rad
          lname_ijts(k) = trim(trname(n))//'3 SW radiative forcing'
          sname_ijts(k) = 'swf_'//trim(trname(n))//'3'
          ijts_power(k) = -2.
          units_ijts(k) = unit_string(ijts_power(k),'W/m2')
          scale_ijts(k) = 10.**(-ijts_power(k))
          k = k + 1
          ijts_fcsub(1,n,4) = k
          ijts_index(k) = n
          ia_ijts(k) = ia_rad
          lname_ijts(k) = trim(trname(n))//'4 SW radiative forcing'
          sname_ijts(k) = 'swf_'//trim(trname(n))//'4'
          ijts_power(k) = -2.
          units_ijts(k) = unit_string(ijts_power(k),'W/m2')
          scale_ijts(k) = 10.**(-ijts_power(k))
c dust longwave radiative forcing of four clay sub size classes
          k = k + 1
          ijts_fcsub(2,n,1) = k
          ijts_index(k) = n
          ia_ijts(k) = ia_rad
          lname_ijts(k) = trim(trname(n))//'1 LW radiative forcing'
          sname_ijts(k) = 'lwf_'//trim(trname(n))//'1'
          ijts_power(k) = -2.
          units_ijts(k) = unit_string(ijts_power(k),'W/m2')
          scale_ijts(k) = 10.**(-ijts_power(k))
          k = k + 1
          ijts_fcsub(2,n,2) = k
          ijts_index(k) = n
          ia_ijts(k) = ia_rad
          lname_ijts(k) = trim(trname(n))//'2 LW radiative forcing'
          sname_ijts(k) = 'lwf_'//trim(trname(n))//'2'
          ijts_power(k) = -2.
          units_ijts(k) = unit_string(ijts_power(k),'W/m2')
          scale_ijts(k) = 10.**(-ijts_power(k))
          k = k + 1
          ijts_fcsub(2,n,3) = k
          ijts_index(k) = n
          ia_ijts(k) = ia_rad
          lname_ijts(k) = trim(trname(n))//'3 LW radiative forcing'
          sname_ijts(k) = 'lwf_'//trim(trname(n))//'3'
          ijts_power(k) = -2.
          units_ijts(k) = unit_string(ijts_power(k),'W/m2')
          scale_ijts(k) = 10.**(-ijts_power(k))
          k = k + 1
          ijts_fcsub(2,n,4) = k
          ijts_index(k) = n
          ia_ijts(k) = ia_rad
          lname_ijts(k) = trim(trname(n))//'4 LW radiative forcing'
          sname_ijts(k) = 'lwf_'//trim(trname(n))//'4'
          ijts_power(k) = -2.
          units_ijts(k) = unit_string(ijts_power(k),'W/m2')
          scale_ijts(k) = 10.**(-ijts_power(k))
c dust shortwave radiative forcing at surface of four clay sub size classes
          k = k + 1
          ijts_fcsub(3,n,1) = k
          ijts_index(k) = n
          ia_ijts(k) = ia_rad
          lname_ijts(k) = trim(trname(n))//'1 SW Surf radiative forcing'
          sname_ijts(k) = 'swf_surf_'//trim(trname(n))//'1'
          ijts_power(k) = -2.
          units_ijts(k) = unit_string(ijts_power(k),'W/m2')
          scale_ijts(k) = 10.**(-ijts_power(k))
          k = k + 1
          ijts_fcsub(3,n,2) = k
          ijts_index(k) = n
          ia_ijts(k) = ia_rad
          lname_ijts(k) = trim(trname(n))//'2 SW Surf radiative forcing'
          sname_ijts(k) = 'swf_surf_'//trim(trname(n))//'2'
          ijts_power(k) = -2.
          units_ijts(k) = unit_string(ijts_power(k),'W/m2')
          scale_ijts(k) = 10.**(-ijts_power(k))
          k = k + 1
          ijts_fcsub(3,n,3) = k
          ijts_index(k) = n
          ia_ijts(k) = ia_rad
          lname_ijts(k) = trim(trname(n))//'3 SW Surf radiative forcing'
          sname_ijts(k) = 'swf_surf_'//trim(trname(n))//'3'
          ijts_power(k) = -2.
          units_ijts(k) = unit_string(ijts_power(k),'W/m2')
          scale_ijts(k) = 10.**(-ijts_power(k))
          k = k + 1
          ijts_fcsub(3,n,4) = k
          ijts_index(k) = n
          ia_ijts(k) = ia_rad
          lname_ijts(k) = trim(trname(n))//'4 SW Surf radiative forcing'
          sname_ijts(k) = 'swf_surf_'//trim(trname(n))//'4'
          ijts_power(k) = -2.
          units_ijts(k) = unit_string(ijts_power(k),'W/m2')
          scale_ijts(k) = 10.**(-ijts_power(k))
c dust longwave radiative forcing at surface of four sub size classes
          k = k + 1
          ijts_fcsub(4,n,1) = k
          ijts_index(k) = n
          ia_ijts(k) = ia_rad
          lname_ijts(k) = trim(trname(n))//'1 LW Surf radiative forcing'
          sname_ijts(k) = 'lwf_surf_'//trim(trname(n))//'1'
          ijts_power(k) = -2.
          units_ijts(k) = unit_string(ijts_power(k),'W/m2')
          scale_ijts(k) = 10.**(-ijts_power(k))
          k = k + 1
          ijts_fcsub(4,n,2) = k
          ijts_index(k) = n
          ia_ijts(k) = ia_rad
          lname_ijts(k) = trim(trname(n))//'2 LW Surf radiative forcing'
          sname_ijts(k) = 'lwf_surf_'//trim(trname(n))//'2'
          ijts_power(k) = -2.
          units_ijts(k) = unit_string(ijts_power(k),'W/m2')
          scale_ijts(k) = 10.**(-ijts_power(k))
          k = k + 1
          ijts_fcsub(4,n,3) = k
          ijts_index(k) = n
          ia_ijts(k) = ia_rad
          lname_ijts(k) = trim(trname(n))//'3 LW Surf radiative forcing'
          sname_ijts(k) = 'lwf_surf_'//trim(trname(n))//'3'
          ijts_power(k) = -2.
          units_ijts(k) = unit_string(ijts_power(k),'W/m2')
          scale_ijts(k) = 10.**(-ijts_power(k))
          k = k + 1
          ijts_fcsub(4,n,4) = k
          ijts_index(k) = n
          ia_ijts(k) = ia_rad
          lname_ijts(k) = trim(trname(n))//'4 LW Surf radiative forcing'
          sname_ijts(k) = 'lwf_surf_'//trim(trname(n))//'4'
          ijts_power(k) = -2.
          units_ijts(k) = unit_string(ijts_power(k),'W/m2')
          scale_ijts(k) = 10.**(-ijts_power(k))
        CASE('Silt1','Silt2','Silt3','Silt4',
     &     'ClayIlli','ClayKaol','ClaySmec','ClayCalc','ClayQuar',
     &     'Sil1Quar','Sil1Feld','Sil1Calc','Sil1Hema','Sil1Gyps',
     &     'Sil2Quar','Sil2Feld','Sil2Calc','Sil2Hema','Sil2Gyps',
     &     'Sil3Quar','Sil3Feld','Sil3Calc','Sil3Hema','Sil3Gyps',
     &     'Sil1QuHe','Sil2QuHe','Sil3QuHe')
c dust optical thickness
          k = k + 1
          ijts_tau(1,n) = k
          ijts_index(k) = n
          ia_ijts(k) = ia_rad   
          lname_ijts(k) = trim(trname(n))//' optical thickness'
          sname_ijts(k) = 'tau_'//trim(trname(n))
          ijts_power(k) = -2.
          units_ijts(k) = unit_string(ijts_power(k),' ')
          scale_ijts(k) = 10.**(-ijts_power(k))
c dust clear sky optical thickness
          k = k + 1
          ijts_tau(2,n) = k
          ijts_index(k) = n
          ia_ijts(k) = ia_rad   
          lname_ijts(k) = trim(trname(n))//' CS optical thickness'
          sname_ijts(k) = 'tau_CS_'//trim(trname(n))
          ijts_power(k) = -2.
          units_ijts(k) = unit_string(ijts_power(k),' ')
          scale_ijts(k) = 10.**(-ijts_power(k))
c dust shortwave radiative forcing
          k = k + 1
          ijts_fc(1,n) = k
          ijts_index(k) = n
          ia_ijts(k) = ia_rad  
          lname_ijts(k) = trim(trname(n))//' SW radiative forcing'
          sname_ijts(k) = 'swf_'//trim(trname(n))
          ijts_power(k) = -2.
          units_ijts(k) = unit_string(ijts_power(k),'W/m2')
          scale_ijts(k) = 10.**(-ijts_power(k))
c dust longwave radiative forcing
          k = k + 1
          ijts_fc(2,n) = k
          ijts_index(k) = n
          ia_ijts(k) = ia_rad  
          lname_ijts(k) = trim(trname(n))//' LW radiative forcing'
          sname_ijts(k) = 'lwf_'//trim(trname(n))
          ijts_power(k) = -2.
          units_ijts(k) = unit_string(ijts_power(k),'W/m2')
          scale_ijts(k) = 10.**(-ijts_power(k))
c dust shortwave radiative forcing at surface
          k = k + 1
          ijts_fc(3,n) = k
          ijts_index(k) = n
          ia_ijts(k) = ia_rad  
          lname_ijts(k) = trim(trname(n))//' SW Surf radiative forcing'
          sname_ijts(k) = 'swf_surf_'//trim(trname(n))
          ijts_power(k) = -2.
          units_ijts(k) = unit_string(ijts_power(k),'W/m2')
          scale_ijts(k) = 10.**(-ijts_power(k))
c dust longwave radiative forcing at surface
          k = k + 1
          ijts_fc(4,n) = k
          ijts_index(k) = n
          ia_ijts(k) = ia_rad  
          lname_ijts(k) = trim(trname(n))//' LW Surf radiative forcing'
          sname_ijts(k) = 'lwf_surf_'//trim(trname(n))
          ijts_power(k) = -2.
          units_ijts(k) = unit_string(ijts_power(k),'W/m2')
          scale_ijts(k) = 10.**(-ijts_power(k))
        END SELECT
#endif

      end select
      end do

C**** Additional Special IJ diagnostics
C**** (not necessary associated with a particular tracer)
#ifdef TRACERS_SPECIAL_Shindell
      k = k+1
        ijs_flash=k
        ijts_index(k) = ntm
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'Lightning Flash Rate'
        sname_ijts(k) = 'lightning_flash'
        ijts_power(k) = -10.
        units_ijts(k) = unit_string(ijts_power(k),'flash/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k+1
        ijs_CtoG=k
        ijts_index(k) = ntm
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'Cloud to Ground Lightning Flash Rate'
        sname_ijts(k) = 'CtoG_flash'
        ijts_power(k) = -10.
        units_ijts(k) = unit_string(ijts_power(k),'flash/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
#ifdef regional_Ox_tracers
      k = k+1
        ijs_Oxloss=k
        ijts_index(k) = ntm
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'Ox chemical loss rate'
        sname_ijts(k) = 'Ox_loss'
        ijts_power(k) = -2.
        units_ijts(k) = unit_string(ijts_power(k),'molecules/cm3/s')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k+1
        ijs_Oxprod=k
        ijts_index(k) = ntm
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'Ox chemical production rate'
        sname_ijts(k) = 'Ox_prod'
        ijts_power(k) = -2.
        units_ijts(k) = unit_string(ijts_power(k),'molecules/cm3/s')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
#endif
      do L=1,LTOP
        k = k + 1
          ijs_OH(L)=k
          ijts_index(k) = ntm
          ia_ijts(k) = ia_src
          write(lname_ijts(k),'(a19,i2.2)') 'OH concentration L=',L
          write(sname_ijts(k),'(a7,i2.2)') 'OH_con_',L
          ijts_power(k) = 5.
          units_ijts(k) = unit_string(ijts_power(k),'molecules/cm3')
          scale_ijts(k) = 10.**(-ijts_power(k))
      end do
      do L=1,LTOP
        k = k + 1
          ijs_NO3(L)=k
          ijts_index(k) = ntm
          ia_ijts(k) = ia_src
          write(lname_ijts(k),'(a20,i2.2)') 'NO3 concentration L=',L
          write(sname_ijts(k),'(a8,i2.2)') 'NO3_con_',L
          ijts_power(k) = 5.
          units_ijts(k) = unit_string(ijts_power(k),'molecules/cm3')
          scale_ijts(k) = 10.**(-ijts_power(k))
      end do
      do L=1,LTOP
        k = k + 1
          ijs_HO2(L)=k
          ijts_index(k) = ntm
          ia_ijts(k) = ia_src
          write(lname_ijts(k),'(a20,i2.2)') 'HO2 concentration L=',L
          write(sname_ijts(k),'(a8,i2.2)') 'HO2_con_',L
          ijts_power(k) = 7.
          units_ijts(k) = unit_string(ijts_power(k),'molecules/cm3')
          scale_ijts(k) = 10.**(-ijts_power(k))
      end do
      do L=1,LTOP
        k = k + 1
          ijs_JH2O2(L)=k
          ijts_index(k) = ntm
          ia_ijts(k) = ia_src
          write(lname_ijts(k),'(a23,i2.2)') 'H2O2 photolysis rate L=',L
          write(sname_ijts(k),'(a8,i2.2)') 'J(H2O2)_',L
          ijts_power(k) = 2.
          units_ijts(k) = unit_string(ijts_power(k),'s-1')
          scale_ijts(k) = 10.**(-ijts_power(k))
      end do
#endif
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
      k = k + 1
      ijts_spec(nDustEv1ij)=k
      lname_ijts(k)='No. dust events'
      sname_ijts(k)='no_dust_ev1'
      ijts_index(k)=Ntm_dust
      ia_ijts(k)=ia_src
      scale_ijts(k)=Sday/Dtsrc
      units_ijts(k)='1/d'
      k = k + 1
      ijts_spec(nDustEv2ij)=k
      lname_ijts(k)='No. dust events above threshold wind'
      sname_ijts(k)='no_dust_ev2'
      ijts_index(k)=Ntm_dust
      ia_ijts(k)=ia_src
      scale_ijts(k)=Sday/Dtsrc
      units_ijts(k)='1/d'
      k = k + 1
      ijts_spec(nDustWthij)=k
      lname_ijts(k)='Threshold velocity for dust emission'
      sname_ijts(k)='wtrsh'
      ijts_index(k)=Ntm_dust
      ia_ijts(k)=ia_src
      scale_ijts(k)=1.
      units_ijts(k)='m/s'
#endif
#ifdef TRACERS_DUST
c       do n=1,ntm
c       select case(trname(n))
c     CASE('Clay','Silt1','Silt2','Silt3','Silt4')
c     do L=1,LTOP
c       k = k + 1
c         ijts_3Dtau(l,n)=k
c         ijts_index(k) = n
c         ia_ijts(k) = ia_rad
c         write(lname_ijts(k),'(a12,i2.2)') trname(n)' tau L=    ',L
c         write(sname_ijts(k),'(a12,i2.2)') 'tau_3D_'//trname(n) ,L
c         ijts_power(k) = -2.
c         units_ijts(k) = unit_string(ijts_power(k),' ')
c         scale_ijts(k) = 10.**(-ijts_power(k))
c     end do
c     end select
c     end do
#endif

      if (k .gt. ktaijs) then
        write (6,*)'ijt_defs: Increase ktaijs=',ktaijs,' to at least ',k
        call stop_model('ktaijs too small',255)
      end if

C**** Initialize conservation diagnostics
C**** To add a new conservation diagnostic:
C****       Set up a QCON, and call SET_TCON to allocate array numbers,
C****       set up scales, titles, etc.
C**** QCON denotes when the conservation diags should be accumulated
C**** QSUM says whether that diag is to be used in summation (if the
C****      routine DIAGTCB is used, this must be false).
C**** 1:NPTS+1 ==> INST,  DYN,   COND,   RAD,   PREC,   LAND,  SURF,
C****            FILTER,STRDG/OCEAN, DAILY, OCEAN1, OCEAN2,
C**** First 12 are standard for all tracers and GCM
      QCON=(/ t,                                           !instant.
     *        T,  T,  F,  F,  T,  T,  T,  T,  F,  F,  F,   !2-12 (npts)
     *        F,  F,  F,  F,  F,  F,  F,  F,  F,  F,       !13-22
     *        F,  F,  F,  F,  F,  F,  F,  F,  F,  F  /)    !21-ktcon-1
      QSUM=(/ f,                                           !instant.
     *        T,  T,  F,  F,  T,  T,  T,  T,  F,  F,  F,   !2-12 (npts)
     *        F,  F,  F,  F,  F,  F,  F,  F,  F,  F,       !13-22
     *        F,  F,  F,  F,  F,  F,  F,  F,  F,  F  /)    !21-ktcon-1
      do n=1,ntm
        kt_power_inst(n)   = ntm_power(n)+2
        kt_power_change(n) = ntm_power(n)-4
        scale_inst(n)   = 10d0**(-kt_power_inst(n))
        scale_change(n) = 10d0**(-kt_power_change(n))
        inst_unit(n) = unit_string(kt_power_inst(n),  'kg/m^2)')
         sum_unit(n) = unit_string(kt_power_change(n),'kg/m^2 s)')
      end do

C**** set some defaults
      itcon_mc(:)=0
      itcon_ss(:)=0
      itcon_surf(:,:)=0
      itcon_3Dsrc(:,:)=0
      itcon_decay(:)=0
#ifdef TRACERS_DRYDEP
      itcon_dd(:,:)=0
#endif

      k = 0
      do n=1,ntm
      select case (trname(n))

      case ('Air')
      CALL SET_TCON(QCON,TRNAME(N),QSUM,inst_unit(n),
     *     sum_unit(n),scale_inst(n),scale_change(n), N,CONPTs)

      case ('SF6')
      CALL SET_TCON(QCON,TRNAME(N),QSUM,inst_unit(n),
     *     sum_unit(n),scale_inst(n),scale_change(n), N,CONPTs)

      case ('Rn222')
      itcon_decay(n) = 13
      qcon(itcon_decay(n)) = .true.; conpts(1) = 'DECAY'
      qsum(itcon_decay(n)) = .true.
      CALL SET_TCON(QCON,TRNAME(N),QSUM,inst_unit(n),
     *     sum_unit(n),scale_inst(n),scale_change(n), N,CONPTs)

      case ('CO2')
      itcon_surf(1,N) = 13
      qcon(itcon_surf(1,N)) = .true.; conpts(1) = 'FossilFuel'
      itcon_surf(2,N) = 14
      qcon(itcon_surf(2,N)) = .true.; conpts(2) = 'Fertilization'
      itcon_surf(3,N) = 15
      qcon(itcon_surf(3,N)) = .true.; conpts(3) = 'Forest Regrowth'
      itcon_surf(4,N) = 16
      qcon(itcon_surf(4,N)) = .true.; conpts(4) = 'Land Use'
      itcon_surf(5,N) = 17
      qcon(itcon_surf(5,N)) = .true.; conpts(5) = 'Ecosystem Exch'
      itcon_surf(6,N) = 18
      qcon(itcon_surf(6,N)) = .true.; conpts(6) = 'Ocean Exch'
      qsum(itcon_surf(1:6,N)) = .false.  ! prevent summing twice
      CALL SET_TCON(QCON,TRNAME(N),QSUM,inst_unit(n),
     *     sum_unit(n),scale_inst(n),scale_change(n), N,CONPTs)
      qcon(13:) = .false.  ! reset to defaults for next tracer

      case ('N2O')
#ifdef TRACERS_SPECIAL_Lerner
      itcon_surf(1,N) = 13
      qcon(itcon_surf(1,N)) = .true.; conpts(1) = 'Reset in L1'
      itcon_3Dsrc(1,N) = 14
      qcon(itcon_3Dsrc(1,N)) = .true.; conpts(2) = 'Strat. Chem'
      qsum(itcon_3Dsrc(1,N)) = .true.
      CALL SET_TCON(QCON,TRNAME(N),QSUM,inst_unit(n),
     *     sum_unit(n),scale_inst(n),scale_change(n), N,CONPTs)
      qsum(13:) = .false.  ! reset to defaults for next tracer
#endif
#ifdef TRACERS_SPECIAL_Shindell
      kt_power_change(n) = -14
      scale_change(n) = 10d0**(-kt_power_change(n))
      sum_unit(n) = unit_string(kt_power_change(n),'kg/m^2 s)')
      g=13; itcon_3Dsrc(1,N) = g
      qcon(itcon_3Dsrc(1,N)) = .true.; conpts(g-12)= 'Chemistry'
      qsum(itcon_3Dsrc(1,N)) = .true.
#ifdef TRACERS_DRYDEP
      if(dodrydep(n)) then
        g=g+1; itcon_dd(n,1)=g
        qcon(itcon_dd(n,1)) = .true. ; conpts(g-12) = 'DRY DEP'
        qsum(itcon_dd(n,1)) = .false.
      end if
#endif
      CALL SET_TCON(QCON,TRNAME(N),QSUM,inst_unit(n),
     *     sum_unit(n),scale_inst(n),scale_change(n), N,CONPTs)
      qcon(13:) = .false.  ! reset to defaults for next tracer
      qsum(13:) = .false.  ! reset to defaults for next tracer
#endif

      case ('CFC11')
      itcon_surf(1,N) = 13
      qcon(itcon_surf(1,N)) = .true.; conpts(1) = 'L1 Source'
      itcon_3Dsrc(1,N) = 14
      qcon(itcon_3Dsrc(1,N)) = .true.; conpts(2) = 'Strat. Chem'
      qsum(itcon_3Dsrc(1,N)) = .true.
      CALL SET_TCON(QCON,TRNAME(N),QSUM,inst_unit(n),
     *     sum_unit(n),scale_inst(n),scale_change(n), N,CONPTs)
      qcon(13:) = .false.  ! reset to defaults for next tracer
      qsum(13:) = .false.  ! reset to defaults for next tracer

      case ('14CO2')
      itcon_surf(1,N) = 13
      qcon(itcon_surf(1,N)) = .true.; conpts(1) = 'Bombs and drift'
      CALL SET_TCON(QCON,TRNAME(N),QSUM,inst_unit(n),
     *     sum_unit(n),scale_inst(n),scale_change(n), N,CONPTs)
      qcon(13:) = .false.  ! reset to defaults for next tracer
      qsum(13:) = .false.  ! reset to defaults for next tracer

      case ('CH4')
#ifdef TRACERS_SPECIAL_Shindell
      kt_power_change(n) = -13
      scale_change(n) = 10d0**(-kt_power_change(n))
      sum_unit(n) = unit_string(kt_power_change(n),'kg/m^2 s)')
      g=13; itcon_3Dsrc(nChemistry,N) = g
      qcon(itcon_3Dsrc(nChemistry,N)) =.true.; conpts(g-12)='Chemistry'
      qsum(itcon_3Dsrc(nChemistry,N)) = .true.
#ifndef SHINDELL_STRAT_CHEM
      g=g+1; itcon_3Dsrc(nStratwrite,N) = g
      qcon(itcon_3Dsrc(nStratwrite,N))=.true.
      conpts(g-12)='Strat Overwrite'
      qsum(itcon_3Dsrc(nStratwrite,N)) = .true.
#endif
      g=g+1; itcon_surf(6,N) = g; qcon(itcon_surf(6,N)) = .true.
      conpts(g-12)= 'Soil sink'
      g=g+1; itcon_surf(7,N) = g; qcon(itcon_surf(7,N)) = .true.
      conpts(g-12)= 'Termite Source'
      g=g+1; itcon_surf(9,N) = g; qcon(itcon_surf(9,N)) = .true.
      conpts(g-12)= 'Ocean source'
      g=g+1; itcon_surf(10,N) = g; qcon(itcon_surf(10,N)) = .true.
      conpts(g-12)= 'Lake source'
      g=g+1; itcon_surf(11,N) = g; qcon(itcon_surf(11,N)) = .true.
      conpts(g-12)='Misc. Ground source'      
      g=g+1; itcon_surf(14,N) = g; qcon(itcon_surf(14,N)) = .true.
      conpts(g-12)= 'Wetlands+Tundra'
#ifdef EDGAR_HYDE_SOURCES
      g=g+1; itcon_surf(1,N) = g; qcon(itcon_surf(1,N)) = .true.
      conpts(g-12)= 'EH Animal source'
      g=g+1; itcon_surf(2,N) = g; qcon(itcon_surf(2,N)) = .true.
      conpts(g-12)= 'EH ffuel comb source'
      g=g+1; itcon_surf(3,N) = g; qcon(itcon_surf(3,N)) = .true.
      conpts(g-12)= 'EH landfills source'
      g=g+1; itcon_surf(4,N) = g; qcon(itcon_surf(4,N)) = .true.
      conpts(g-12)= 'EH ffuel prod source'
      g=g+1; itcon_surf(5,N) = g; qcon(itcon_surf(5,N)) = .true.
      conpts(g-12)= 'EH Ag Waste Burn source'
      g=g+1; itcon_surf(8,N) = g; qcon(itcon_surf(8,N)) = .true.
      conpts(g-12)= 'EH deforestation source'
      g=g+1; itcon_surf(12,N) = g; qcon(itcon_surf(12,N)) = .true.
      conpts(g-12)= 'EH savann. burn source'
      g=g+1; itcon_surf(13,N) = g; qcon(itcon_surf(13,N)) = .true.
      conpts(g-12)= 'EH biofuel PTC source'
      g=g+1; itcon_surf(15,N) = g; qcon(itcon_surf(15,N)) = .true.
      conpts(g-12)= 'EH Ag Land Act source'
      qsum(itcon_surf(1:15,N)) = .false.  ! prevent summing twice
#else
      g=g+1; itcon_surf(1,N) = g
      qcon(itcon_surf(1,N)) = .true.; conpts(g-12)= 'Animal source'
      g=g+1; itcon_surf(2,N) = g
      qcon(itcon_surf(2,N)) = .true.; conpts(g-12)= 'Coal Mine source'
      g=g+1; itcon_surf(3,N) = g
      qcon(itcon_surf(3,N)) = .true.; conpts(g-12)= 'Gas Leak source'
      g=g+1; itcon_surf(4,N) = g
      qcon(itcon_surf(4,N)) = .true.; conpts(g-12)= 'Gas Vent source'
      g=g+1; itcon_surf(5,N) = g
      qcon(itcon_surf(5,N)) = .true.; conpts(g-12)= 'City Dump source'
      g=g+1; itcon_surf(8,N) = g
      qcon(itcon_surf(8,N)) = .true.; conpts(g-12)= 'Coal Combustion'
      g=g+1; itcon_surf(12,N) = g
      qcon(itcon_surf(12,N)) = .true.; conpts(g-12)= 'Biomass Burning'
      g=g+1; itcon_surf(13,N) = g
      qcon(itcon_surf(13,N)) = .true.; conpts(g-12)= 'Rice source'
      qsum(itcon_surf(1:14,N)) = .false.  ! prevent summing twice
#endif
#ifdef TRACERS_DRYDEP
      if(dodrydep(n)) then
        g=g+1; itcon_dd(n,1)= g
        qcon(itcon_dd(n,1)) = .true. ; conpts(g-12) = 'DRY DEP'
        qsum(itcon_dd(n,1)) = .false.
      end if
#endif
      CALL SET_TCON(QCON,TRNAME(N),QSUM,inst_unit(n),
     *     sum_unit(n),scale_inst(n),scale_change(n), N,CONPTs)
      qcon(13:) = .false.  ! reset to defaults for next tracer
      qsum(13:) = .false.  ! reset to defaults for next tracer
#else
      itcon_surf(1,N) = 13
      qcon(itcon_surf(1,N)) = .true.; conpts(1) = 'Animal source'
      itcon_surf(2,N) = 14
      qcon(itcon_surf(2,N)) = .true.; conpts(2) = 'Coal Mine source'
      itcon_surf(3,N) = 15
      qcon(itcon_surf(3,N)) = .true.; conpts(3) = 'Gas Leak source'
      itcon_surf(4,N) = 16
      qcon(itcon_surf(4,N)) = .true.; conpts(4) = 'Gas Vent source'
      itcon_surf(5,N) = 17
      qcon(itcon_surf(5,N)) = .true.; conpts(5) = 'City Dump source'
      itcon_surf(6,N) = 18
      qcon(itcon_surf(6,N)) = .true.; conpts(6) = 'Soil sink'
      itcon_surf(7,N) = 19
      qcon(itcon_surf(7,N)) = .true.; conpts(7) = 'Termite Source'
      itcon_surf(8,N) = 20
      qcon(itcon_surf(8,N)) = .true.; conpts(8) = 'Coal Combustion'
      itcon_surf(9,N) = 21
      qcon(itcon_surf(9,N)) = .true.; conpts(9) = 'Ocean source'
      itcon_surf(10,N) = 22
      qcon(itcon_surf(10,N)) = .true.; conpts(10) = 'Lake source'
      itcon_surf(11,N) = 23
      qcon(itcon_surf(11,N)) = .true.; conpts(11) ='Misc. Ground source'
      itcon_surf(12,N) = 24
      qcon(itcon_surf(12,N)) = .true.; conpts(12) = 'Biomass Burning'
      itcon_surf(13,N) = 25
      qcon(itcon_surf(13,N)) = .true.; conpts(13) = 'Rice source'
      itcon_surf(14,N) = 26
      qcon(itcon_surf(14,N)) = .true.; conpts(14) = 'Wetlands+Tundra'
      qsum(itcon_surf(1:14,N)) = .false.  ! prevent summing twice
      itcon_3Dsrc(1,N) = 27
      qcon(itcon_3Dsrc(1,N)) = .true.; conpts(15) = 'Tropos. Chem'
      qsum(itcon_3Dsrc(1,N)) = .true.
      itcon_3Dsrc(2,N) = 28
      qcon(itcon_3Dsrc(2,N)) = .true.; conpts(16) = 'Stratos. Chem'
      qsum(itcon_3Dsrc(2,N)) = .true.
      CALL SET_TCON(QCON,TRNAME(N),QSUM,inst_unit(n),
     *     sum_unit(n),scale_inst(n),scale_change(n), N,CONPTs)
      qcon(13:) = .false.  ! reset to defaults for next tracer
      qsum(13:) = .false.  ! reset to defaults for next tracer
#endif

      case ('O3')
      itcon_3Dsrc(1,N) = 13
      qcon(itcon_3Dsrc(1,N)) = .true.; conpts(1) = 'Tropos. Chem'
      qsum(itcon_3Dsrc(1,N)) = .true.
      itcon_3Dsrc(2,N) = 14
      qcon(itcon_3Dsrc(2,N)) = .true.; conpts(2) = 'Stratos. Chem'
      qsum(itcon_3Dsrc(2,N)) = .true.
      CALL SET_TCON(QCON,TRNAME(N),QSUM,inst_unit(n),
     *     sum_unit(n),scale_inst(n),scale_change(n), N,CONPTs)
      qcon(13:) = .false.  ! reset to defaults for next tracer
      qsum(13:) = .false.  ! reset to defaults for next tracer


      case ('Ox','OxREG1','OxREG2','OxREG3','OxREG4','OxREG5','OxREG6',
     &      'N2O5','HNO3','H2O2','CH3OOH','HCHO','HO2NO2','PAN',
     &      'AlkylNit','ClOx','BrOx','HCl','HOCl','ClONO2','HBr',
     &      'HOBr','BrONO2','CFC') ! N2O done above
       select case (trname(n)) 
       case ('N2O5','CH3OOH','HCHO','HO2NO2','PAN','AlkylNit',
     &       'ClOx','BrOx','HCl','HOCl','ClONO2','HBr',
     &       'HOBr','BrONO2','CFC')
         kt_power_change(n) = -14
       case ('HNO3','H2O2')
         kt_power_change(n) = -13
       case default
         kt_power_change(n) = -12 
       end select
      scale_change(n) = 10d0**(-kt_power_change(n))
      sum_unit(n) = unit_string(kt_power_change(n),'kg/m^2 s)')
      g=13; itcon_3Dsrc(nChemistry,N) = g
      qcon(itcon_3Dsrc(nChemistry,N)) = .true.
      conpts(g-12) = 'Chemistry'
      qsum(itcon_3Dsrc(nChemistry,N)) = .true.
#ifndef SHINDELL_STRAT_CHEM
      g=g+1; itcon_3Dsrc(nStratwrite,N) = g
      qcon(itcon_3Dsrc(nStratwrite,N)) = .true.
      conpts(g-12)='Strat Overwrite'
      qsum(itcon_3Dsrc(nStratwrite,N)) = .true.
#endif
#ifdef TRACERS_WATER
      if(dowetdep(n)) then
        g=g+1; itcon_mc(n) = g
        qcon(itcon_mc(n)) = .true.  ; conpts(g-12) = 'MOIST CONV'
        g=g+1; itcon_ss(n) = g
        qcon(itcon_ss(n)) = .true.  ; conpts(g-12) = 'LS COND'
      end if
#endif
#ifdef TRACERS_DRYDEP
      if(dodrydep(n)) then
        g=g+1; itcon_dd(n,1)=g
        qcon(itcon_dd(n,1)) = .true. ; conpts(g-12) = 'DRY DEP'
      end if
#endif
      CALL SET_TCON(QCON,TRNAME(N),QSUM,inst_unit(n),
     *     sum_unit(n),scale_inst(n),scale_change(n), N,CONPTs)
      qcon(13:) = .false.  ! reset to defaults for next tracer
      qsum(13:) = .false.  ! reset to defaults for next tracer


      case ('NOx')
      kt_power_change(n) = -14
      scale_change(n) = 10d0**(-kt_power_change(n))
      sum_unit(n) = unit_string(kt_power_change(n),'kg/m^2 s)')
      g=13; itcon_3Dsrc(nChemistry,N) = g
      qcon(itcon_3Dsrc(nChemistry,N)) = .true.
      conpts(g-12) = 'Chemistry'
      qsum(itcon_3Dsrc(nChemistry,N)) = .true.
#ifndef SHINDELL_STRAT_CHEM
      g=g+1; itcon_3Dsrc(nStratwrite,N) = g
      qcon(itcon_3Dsrc(nStratwrite,N)) = .true.
      conpts(g-12)='Strat Overwrite'
      qsum(itcon_3Dsrc(nStratwrite,N)) = .true.
#endif
      g=g+1; itcon_3Dsrc(nLightning,N) = g
      qcon(itcon_3Dsrc(nLightning,N)) = .true.
      conpts(g-12) = 'Lightning'
      qsum(itcon_3Dsrc(nLightning,N)) = .true.
      g=g+1; itcon_3Dsrc(nAircraft,N) = g
      qcon(itcon_3Dsrc(nAircraft,N)) = .true.
      conpts(g-12) = 'Aircraft'
      qsum(itcon_3Dsrc(nAircraft,N)) = .true.
#ifdef EDGAR_HYDE_SOURCES
      g=g+1; itcon_surf(1,N) = g; qcon(itcon_surf(1,N)) = .true.
      conpts(g-12) = 'EH Fossil Fuel Comb.'
      g=g+1; itcon_surf(2,N) = g; qcon(itcon_surf(2,N)) = .true.
      conpts(g-12) = 'EH Indust. Proc. Comb.'
      g=g+1; itcon_surf(3,N) = g; qcon(itcon_surf(3,N)) = .true.
      conpts(g-12) = 'EH Ag Land Activities'
      g=g+1; itcon_surf(4,N) = g; qcon(itcon_surf(4,N)) = .true.
      conpts(g-12) = 'EH Ag Waste Burning'
      g=g+1; itcon_surf(5,N) = g; qcon(itcon_surf(5,N)) = .true.
      conpts(g-12) = 'EH Biofuel PTC'
      g=g+1; itcon_surf(6,N) = g; qcon(itcon_surf(6,N)) = .true.
      conpts(g-12) = 'EH Deforestation'
      g=g+1; itcon_surf(7,N) = g; qcon(itcon_surf(7,N)) = .true.
      conpts(g-12) = 'EH Savannah Burning'
#else
      g=g+1; itcon_surf(1,N) = g
      qcon(itcon_surf(1,N)) = .true.; conpts(g-12) = 'Fossil Fuels'
      g=g+1; itcon_surf(2,N) = g
      qcon(itcon_surf(2,N)) = .true.; conpts(g-12) = 'Biomass Burning'
      g=g+1; itcon_surf(3,N) = g
      qcon(itcon_surf(3,N)) = .true.; conpts(g-12) = 'Soil'
#endif
#ifdef TRACERS_DRYDEP
      if(dodrydep(n)) then
        g=g+1; itcon_dd(n,1)=g
        qcon(itcon_dd(n,1)) = .true. ; conpts(g-12) = 'DRY DEP'
      end if
#endif
      CALL SET_TCON(QCON,TRNAME(N),QSUM,inst_unit(n),
     *     sum_unit(n),scale_inst(n),scale_change(n), N,CONPTs)
      qcon(13:) = .false.  ! reset to defaults for next tracer
      qsum(13:) = .false.  ! reset to defaults for next tracer


      case ('CO')
      kt_power_change(n) = -13
      scale_change(n) = 10d0**(-kt_power_change(n))
      sum_unit(n) = unit_string(kt_power_change(n),'kg/m^2 s)')
      g=13; itcon_3Dsrc(nChemistry,N) = g
      qcon(itcon_3Dsrc(nChemistry,N)) = .true.
      conpts(g-12) = 'Chemistry'
      qsum(itcon_3Dsrc(nChemistry,N)) = .true.
#ifndef SHINDELL_STRAT_CHEM
      g=g+1; itcon_3Dsrc(nStratwrite,N) = g
      qcon(itcon_3Dsrc(nStratwrite,N)) = .true.
      conpts(g-12)='Strat Overwrite'
      qsum(itcon_3Dsrc(nStratwrite,N)) = .true.
#endif
#ifdef EDGAR_HYDE_SOURCES
      g=g+1; itcon_surf(1,N) = g; qcon(itcon_surf(1,N)) = .true.
      conpts(g-12) = 'EH fossil fuel comb.'
      g=g+1; itcon_surf(2,N) = g; qcon(itcon_surf(2,N)) = .true.
      conpts(g-12) = 'EH industrial proc.'
      g=g+1; itcon_surf(3,N) = g; qcon(itcon_surf(3,N)) = .true.
      conpts(g-12) = 'EH ag waste burning'
      g=g+1; itcon_surf(4,N) = g; qcon(itcon_surf(4,N)) = .true.
      conpts(g-12) = 'EH savannah burning'
      g=g+1; itcon_surf(5,N) = g; qcon(itcon_surf(5,N)) = .true.
      conpts(g-12) = 'EH biofuel PTC'
      g=g+1; itcon_surf(6,N) = g; qcon(itcon_surf(6,N)) = .true.
      conpts(g-12) = 'EH deforestation'
#else
      g=g+1; itcon_surf(1,N) = g
      qcon(itcon_surf(1,N)) = .true.; conpts(g-12) = 'Industrial'
      g=g+1; itcon_surf(2,N) = g
      qcon(itcon_surf(2,N)) = .true.; conpts(g-12) = 'Biomass Burning'
#endif
#ifdef TRACERS_DRYDEP
      if(dodrydep(n)) then
        g=g+1; itcon_dd(n,1)=g
        qcon(itcon_dd(n,1)) = .true. ; conpts(g-12) = 'DRY DEP'
        qsum(itcon_dd(n,1)) = .false.
      end if
#endif
      CALL SET_TCON(QCON,TRNAME(N),QSUM,inst_unit(n),
     *     sum_unit(n),scale_inst(n),scale_change(n), N,CONPTs)
      qcon(13:) = .false.  ! reset to defaults for next tracer
      qsum(13:) = .false.  ! reset to defaults for next tracer

      case ('Isoprene')
      kt_power_change(n) = -13
      scale_change(n) = 10d0**(-kt_power_change(n))
      sum_unit(n) = unit_string(kt_power_change(n),'kg/m^2 s)')
      g=13; itcon_3Dsrc(nChemistry,N) = g
      qcon(itcon_3Dsrc(nChemistry,N)) = .true.
      conpts(g-12) = 'Chemistry'
      qsum(itcon_3Dsrc(nChemistry,N)) = .true.
#ifndef SHINDELL_STRAT_CHEM
      g=g+1; itcon_3Dsrc(nStratwrite,N) = g
      qcon(itcon_3Dsrc(nStratwrite,N)) = .true.
      conpts(g-12)='Strat Overwrite'
      qsum(itcon_3Dsrc(nStratwrite,N)) = .true.
#endif
      g=g+1; itcon_surf(1,N) = g
      qcon(itcon_surf(1,N)) = .true.; conpts(g-12) = 'Vegetation'
      qsum(itcon_surf(1,N)) = .false.
#ifdef TRACERS_DRYDEP
      if(dodrydep(n)) then
        g=g+1; itcon_dd(n,1)=g
        qcon(itcon_dd(n,1)) = .true. ; conpts(g-12) = 'DRY DEP'
        qsum(itcon_dd(n,1)) = .false.
      end if
#endif
      CALL SET_TCON(QCON,TRNAME(N),QSUM,inst_unit(n),
     *     sum_unit(n),scale_inst(n),scale_change(n), N,CONPTs)
      qcon(13:) = .false.  ! reset to defaults for next tracer
      qsum(13:) = .false.  ! reset to defaults for next tracer

      case ('Alkenes','Paraffin')
      kt_power_change(n) = -13
      scale_change(n) = 10d0**(-kt_power_change(n))
      sum_unit(n) = unit_string(kt_power_change(n),'kg/m^2 s)')
      g=13; itcon_3Dsrc(nChemistry,N) = g
      qcon(itcon_3Dsrc(nChemistry,N)) = .true.
      conpts(g-12) = 'Chemistry'
      qsum(itcon_3Dsrc(nChemistry,N)) = .true.
#ifndef SHINDELL_STRAT_CHEM
      g=g+1; itcon_3Dsrc(nStratwrite,N) = g
      qcon(itcon_3Dsrc(nStratwrite,N)) = .true.
      conpts(g-12)='Strat Overwrite'
      qsum(itcon_3Dsrc(nStratwrite,N)) = .true.
#endif
      g=g+1; itcon_surf(3,N) = g
      qcon(itcon_surf(3,N)) = .true.; conpts(g-12) = 'Vegetation'
#ifdef EDGAR_HYDE_SOURCES
      g=g+1; itcon_surf(1,N) = g; qcon(itcon_surf(1,N)) = .true.
      conpts(g-12) = 'EH Fossil Fuel Comb'
      g=g+1; itcon_surf(2,N) = g; qcon(itcon_surf(2,N)) = .true.
      conpts(g-12) = 'EH Idustrial Proc'
!! #3 is (common) vegetation source.  See above.
      g=g+1; itcon_surf(4,N) = g; qcon(itcon_surf(4,N)) = .true.
      conpts(g-12) = 'EH Fossil Fuel Prod'
      g=g+1; itcon_surf(5,N) = g; qcon(itcon_surf(5,N)) = .true.
      conpts(g-12) = 'EH Ag Waste Burning'
      g=g+1; itcon_surf(6,N) = g; qcon(itcon_surf(6,N)) = .true.
      conpts(g-12) = 'EH Savanna Burning'
      g=g+1; itcon_surf(7,N) = g; qcon(itcon_surf(7,N)) = .true.
      conpts(g-12) = 'EH Biofuel PTC'
      g=g+1; itcon_surf(8,N) = g; qcon(itcon_surf(8,N)) = .true.
      conpts(g-12) = 'EH Deforestation'
#else
      g=g+1; itcon_surf(1,N) = g
      qcon(itcon_surf(1,N)) = .true.; conpts(g-12) = 'Industrial'
      g=g+1; itcon_surf(2,N) = g
      qcon(itcon_surf(2,N)) = .true.; conpts(g-12) = 'Biomass Burning'
!! #3 is (common) vegetation source.  See above.
#endif      
#ifdef TRACERS_DRYDEP
      if(dodrydep(n)) then
        g=g+1; itcon_dd(n,1)=g
        qcon(itcon_dd(n,1)) = .true. ; conpts(g-12) = 'DRY DEP'
      end if
#endif
      CALL SET_TCON(QCON,TRNAME(N),QSUM,inst_unit(n),
     *     sum_unit(n),scale_inst(n),scale_change(n), N,CONPTs)
      qcon(13:) = .false.  ! reset to defaults for next tracer
      qsum(13:) = .false.  ! reset to defaults for next tracer

      case ('Water')
      itcon_mc(n) = 13
      qcon(itcon_mc(n)) = .true.  ; conpts(1) = 'MOIST CONV'
      qsum(itcon_mc(n)) = .false.
      itcon_ss(n) = 14
      qcon(itcon_ss(n)) = .true.  ; conpts(2) = 'LS COND'
      qsum(itcon_ss(n)) = .false.
      CALL SET_TCON(QCON,TRNAME(N),QSUM,inst_unit(n),
     *     sum_unit(n),scale_inst(n),scale_change(n), N,CONPTs)

      case ('H2O18')
      itcon_mc(n) = 13
      qcon(itcon_mc(n)) = .true.  ; conpts(1) = 'MOIST CONV'
      qsum(itcon_mc(n)) = .false.
      itcon_ss(n) = 14
      qcon(itcon_ss(n)) = .true.  ; conpts(2) = 'LS COND'
      qsum(itcon_ss(n)) = .false.
      CALL SET_TCON(QCON,TRNAME(N),QSUM,inst_unit(n),
     *     sum_unit(n),scale_inst(n),scale_change(n), N,CONPTs)

      case ('HDO')
      itcon_mc(n) = 13
      qcon(itcon_mc(n)) = .true.  ; conpts(1) = 'MOIST CONV'
      qsum(itcon_mc(n)) = .false.
      itcon_ss(n) = 14
      qcon(itcon_ss(n)) = .true.  ; conpts(2) = 'LS COND'
      qsum(itcon_ss(n)) = .false.
      CALL SET_TCON(QCON,TRNAME(N),QSUM,inst_unit(n),
     *     sum_unit(n),scale_inst(n),scale_change(n), N,CONPTs)

      case ('HTO')
      itcon_mc(n) = 13
      qcon(itcon_mc(n)) = .true.  ; conpts(1) = 'MOIST CONV'
      qsum(itcon_mc(n)) = .false.
      itcon_ss(n) = 14
      qcon(itcon_ss(n)) = .true.  ; conpts(2) = 'LS COND'
      qsum(itcon_ss(n)) = .false.
      itcon_decay(n) = 15
      qcon(itcon_decay(n)) = .true.; conpts(3) = 'DECAY'
      qsum(itcon_decay(n)) = .true.
      CALL SET_TCON(QCON,TRNAME(N),QSUM,inst_unit(n),
     *     sum_unit(n),scale_inst(n),scale_change(n), N,CONPTs)

      case ('DMS')
      itcon_surf(1,N) = 13
      qcon(itcon_surf(1,N)) = .true.; conpts(1) = 'Ocean Src'
      qsum(itcon_surf(1,N))=.false.
      itcon_3Dsrc(1,N) = 14
      qcon(itcon_3Dsrc(1,N)) = .true.; conpts(2) = 'Chem'
      qsum(itcon_3Dsrc(1,N))= .true.
      CALL SET_TCON(QCON,TRNAME(N),QSUM,inst_unit(n),
     *     sum_unit(n),scale_inst(n),scale_change(n), N,CONPTs)
      qcon(13:) = .false.  ! reset to defaults for next tracer
      qsum(13:) = .false.  ! reset to defaults for next tracer

      case ('MSA')
      itcon_3Dsrc(1,N) = 13
      qcon(itcon_3Dsrc(1,N)) = .true.; conpts(1) = 'Chem'
      qsum(itcon_3Dsrc(1,N)) = .true.
      itcon_mc(n) = 14
      qcon(itcon_mc(n)) = .true.  ; conpts(2) = 'MOIST CONV'
      qsum(itcon_mc(n)) = .false.
      itcon_ss(n) = 15
      qcon(itcon_ss(n)) = .true.  ; conpts(3) = 'LS COND'
      qsum(itcon_ss(n)) = .false.
#ifdef TRACERS_DRYDEP
      if(dodrydep(n)) then
        itcon_dd(n,1)=16 
        qcon(itcon_dd(n,1)) = .true. ; conpts(4) = 'TURB DEP'
        qsum(itcon_dd(n,1)) = .false.
        itcon_dd(n,2)=17
        qcon(itcon_dd(n,2)) = .true. ; conpts(5) = 'GRAV SET'
        qsum(itcon_dd(n,2)) = .false.
      end if
#endif
      CALL SET_TCON(QCON,TRNAME(N),QSUM,inst_unit(n),
     *     sum_unit(n),scale_inst(n),scale_change(n), N,CONPTs)
      qcon(13:) = .false.  ! reset to defaults for next tracer
      qsum(13:) = .false.  ! reset to defaults for next tracer

      case ('SO2')
      itcon_3Dsrc(1,N) = 13
      qcon(itcon_3Dsrc(1,N)) = .true.; conpts(1) = 'Volcanic src'
      qsum(itcon_3Dsrc(1,N)) = .true.
      itcon_3Dsrc(2,N) = 14
      qcon(itcon_3Dsrc(2,N)) = .true.; conpts(2) = 'Aircraft src'
      qsum(itcon_3Dsrc(2,N))=.true.
      itcon_3Dsrc(3,N) = 15
      qcon(itcon_3Dsrc(3,N)) = .true.; conpts(3) = 'Biomass src'
      qsum(itcon_3Dsrc(3,N))=.true.
      itcon_3Dsrc(4,N) = 16
      qcon(itcon_3Dsrc(4,N)) = .true.; conpts(4) = 'Chem src'
      qsum(itcon_3Dsrc(4,N)) = .true.
      itcon_3Dsrc(5,N) = 17
      qcon(itcon_3Dsrc(5,N)) = .true.; conpts(5) = 'Chem sink'
      qsum(itcon_3Dsrc(5,N)) = .true.
      itcon_3Dsrc(6,N) = 18
      qcon(itcon_3Dsrc(6,N)) = .true.; conpts(6) = 'Heter sink'
      qsum(itcon_3Dsrc(6,N)) = .true.
      itcon_surf(1,N) = 19
      qcon(itcon_surf(1,N)) = .true.; conpts(7) = 'Industrial src'
      qsum(itcon_surf(1,N))=.false.
      itcon_mc(n) =20
      qcon(itcon_mc(n)) = .true.  ; conpts(8) = 'MOIST CONV'
      qsum(itcon_mc(n)) = .false.
      itcon_ss(n) =21
      qcon(itcon_ss(n)) = .true.  ; conpts(9) = 'LS COND'
      qsum(itcon_ss(n)) = .false.
#ifdef TRACERS_DRYDEP
      if(dodrydep(n)) then
        itcon_dd(n,1)=22
        qcon(itcon_dd(n,1)) = .true. ; conpts(10) = 'TURB DEP'
        qsum(itcon_dd(n,1)) = .false.
      end if
#endif

      CALL SET_TCON(QCON,TRNAME(N),QSUM,inst_unit(n),
     *     sum_unit(n),scale_inst(n),scale_change(n), N,CONPTs)
      qcon(13:) = .false.  ! reset to defaults for next tracer
      qsum(13:) = .false.  ! reset to defaults for next tracer

      case ('SO4')
      itcon_3Dsrc(1,N) = 13
      qcon(itcon_3Dsrc(1,N)) = .true.; conpts(1) = 'Gas phase src'
      qsum(itcon_3Dsrc(1,N)) = .true.
      itcon_3Dsrc(2,N) = 14
      qcon(itcon_3Dsrc(2,N)) = .true.; conpts(2) = 'Heter src'
      qsum(itcon_3Dsrc(2,N)) = .true.
      itcon_3Dsrc(3,N) = 15
      qcon(itcon_3Dsrc(3,N)) = .true.; conpts(3) = 'volcanic src'
      qsum(itcon_3Dsrc(3,N)) = .true.
      itcon_3Dsrc(4,N) = 16
      qcon(itcon_3Dsrc(4,N)) = .true.; conpts(4) = 'biomass src'
      qsum(itcon_3Dsrc(4,N)) = .true.
      itcon_surf(1,N) = 17
      qcon(itcon_surf(1,N)) = .true.; conpts(5) = 'Industrial src'
      qsum(itcon_surf(1,N)) = .false.
      itcon_mc(n) =18
      qcon(itcon_mc(n)) = .true.  ; conpts(6) = 'MOIST CONV'
      qsum(itcon_mc(n)) = .false.
      itcon_ss(n) =19
      qcon(itcon_ss(n)) = .true.  ; conpts(7) = 'LS COND'
      qsum(itcon_ss(n)) = .false.
#ifdef TRACERS_DRYDEP
      if(dodrydep(n)) then
        itcon_dd(n,1)=20
        qcon(itcon_dd(n,1)) = .true. ; conpts(8) = 'TURB DEP'
        qsum(itcon_dd(n,1)) = .false.
        itcon_dd(n,2)=21
        qcon(itcon_dd(n,2)) = .true. ; conpts(9) = 'GRAV SET'
        qsum(itcon_dd(n,2)) = .false.
      end if
#endif
      CALL SET_TCON(QCON,TRNAME(N),QSUM,inst_unit(n),
     *     sum_unit(n),scale_inst(n),scale_change(n), N,CONPTs)
      qcon(13:) = .false.  ! reset to defaults for next tracer
      qsum(13:) = .false.  ! reset to defaults for next tracer

      case ('BCII')
      itcon_3Dsrc(1,N) = 13
      qcon(itcon_3Dsrc(1,N)) = .true.; conpts(1) = 'Aging loss'
      qsum(itcon_3Dsrc(1,N)) = .true.
      itcon_3Dsrc(2,N) = 14
      qcon(itcon_3Dsrc(2,N)) = .true.; conpts(2) = 'Aircraft source'
      qsum(itcon_3Dsrc(2,N)) = .true.
      itcon_surf(1,N) = 15
      qcon(itcon_surf(1,N)) = .true.; conpts(3) = 'Industrial src'
      qsum(itcon_surf(1,N)) = .false.
      itcon_mc(n) = 16
      qcon(itcon_mc(n)) = .true.  ; conpts(4) = 'MOIST CONV'
      qsum(itcon_mc(n)) = .false.
      itcon_ss(n) = 17
      qcon(itcon_ss(n)) = .true.  ; conpts(5) = 'LS COND'
      qsum(itcon_ss(n)) = .false.
#ifdef TRACERS_DRYDEP
      if(dodrydep(n)) then
        itcon_dd(n,1)=18
        qcon(itcon_dd(n,1)) = .true. ; conpts(6) = 'TURB DEP'
        qsum(itcon_dd(n,1)) = .false.
        itcon_dd(n,2)=19
        qcon(itcon_dd(n,2)) = .true. ; conpts(7) = 'GRAV SET'
        qsum(itcon_dd(n,2)) = .false.
      end if
#endif
      CALL SET_TCON(QCON,TRNAME(N),QSUM,inst_unit(n),
     *     sum_unit(n),scale_inst(n),scale_change(n), N,CONPTs)
      qcon(13:) = .false.  ! reset to defaults for next tracer
      qsum(13:) = .false.  ! reset to defaults for next tracer

      case ('BCIA')
      itcon_3Dsrc(1,N) = 13
      qcon(itcon_3Dsrc(1,N)) = .true.; conpts(1) = 'Aging source'
      qsum(itcon_3Dsrc(1,N)) = .true.
      itcon_3Dsrc(2,N) = 14
      qcon(itcon_3Dsrc(1,N)) = .true.; conpts(2) = 'Aircraft source'
      qsum(itcon_3Dsrc(1,N)) = .true.
      itcon_mc(n) = 15
      qcon(itcon_mc(n)) = .true.  ; conpts(3) = 'MOIST CONV'
      qsum(itcon_mc(n)) = .false.
      itcon_ss(n) = 16
      qcon(itcon_ss(n)) = .true.  ; conpts(4) = 'LS COND'
      qsum(itcon_ss(n)) = .false.
#ifdef TRACERS_DRYDEP
      if(dodrydep(n)) then
        itcon_dd(n,1)=17
        qcon(itcon_dd(n,1)) = .true. ; conpts(5) = 'TURB DEP'
        qsum(itcon_dd(n,1)) = .false.
        itcon_dd(n,2)=18
        qcon(itcon_dd(n,2)) = .true. ; conpts(6) = 'GRAV SET'
        qsum(itcon_dd(n,2)) = .false.
      end if
#endif
      CALL SET_TCON(QCON,TRNAME(N),QSUM,inst_unit(n),
     *     sum_unit(n),scale_inst(n),scale_change(n), N,CONPTs)
      qcon(13:) = .false.  ! reset to defaults for next tracer
      qsum(13:) = .false.  ! reset to defaults for next tracer

      case ('OCIA')
      itcon_3Dsrc(1,N) = 13
      qcon(itcon_3Dsrc(1,N)) = .true.; conpts(1) = 'Aging source'
      qsum(itcon_3Dsrc(1,N)) = .true.
      itcon_mc(n) = 14
      qcon(itcon_mc(n)) = .true.  ; conpts(3) = 'MOIST CONV'
      qsum(itcon_mc(n)) = .false.
      itcon_ss(n) = 15
      qcon(itcon_ss(n)) = .true.  ; conpts(4) = 'LS COND'
      qsum(itcon_ss(n)) = .false.
#ifdef TRACERS_DRYDEP
      if(dodrydep(n)) then
        itcon_dd(n,1)=16
        qcon(itcon_dd(n,1)) = .true. ; conpts(5) = 'TURB DEP'
        qsum(itcon_dd(n,1)) = .false.
        itcon_dd(n,2)=17
        qcon(itcon_dd(n,2)) = .true. ; conpts(6) = 'GRAV SET'
        qsum(itcon_dd(n,2)) = .false.
      end if
#endif
      CALL SET_TCON(QCON,TRNAME(N),QSUM,inst_unit(n),
     *     sum_unit(n),scale_inst(n),scale_change(n), N,CONPTs)
      qcon(13:) = .false.  ! reset to defaults for next tracer
      qsum(13:) = .false.  ! reset to defaults for next tracer

      case ('BCB', 'OCB')
      itcon_3Dsrc(1,N) = 13
      qcon(itcon_3Dsrc(1,N)) = .true.; conpts(1) = 'Biomass src'
      qsum(itcon_3Dsrc(1,N)) = .true.
      itcon_mc(n) = 14
      qcon(itcon_mc(n)) = .true.  ; conpts(2) = 'MOIST CONV'
      qsum(itcon_mc(n)) = .false.
      itcon_ss(n) = 15
      qcon(itcon_ss(n)) = .true.  ; conpts(3) = 'LS COND'
      qsum(itcon_ss(n)) = .false.
#ifdef TRACERS_DRYDEP
      if(dodrydep(n)) then
        itcon_dd(n,1)=16
        qcon(itcon_dd(n,1)) = .true. ; conpts(4) = 'TURB DEP'
        qsum(itcon_dd(n,1)) = .false.
        itcon_dd(n,2)=17
        qcon(itcon_dd(n,2)) = .true. ; conpts(5) = 'GRAV SET'
        qsum(itcon_dd(n,2)) = .false.
      end if
#endif
      CALL SET_TCON(QCON,TRNAME(N),QSUM,inst_unit(n),
     *     sum_unit(n),scale_inst(n),scale_change(n), N,CONPTs)
      qcon(13:) = .false.  ! reset to defaults for next tracer
      qsum(13:) = .false.  ! reset to defaults for next tracer

      case ('OCII')
      itcon_3Dsrc(1,N) = 13
      qcon(itcon_3Dsrc(1,N)) = .true.; conpts(1) = 'Aging loss'
      qsum(itcon_3Dsrc(1,N)) = .true.
      itcon_surf(1,N) = 14
      qcon(itcon_surf(1,N)) = .true.; conpts(2) = 'Industrial src'
      qsum(itcon_surf(1,N)) = .false.
      itcon_surf(2,N) = 15
      qcon(itcon_surf(2,N)) = .true.; conpts(3) = 'Terpene src'
      qsum(itcon_surf(2,N)) = .false.
      itcon_mc(n) = 16
      qcon(itcon_mc(n)) = .true.  ; conpts(4) = 'MOIST CONV'
      qsum(itcon_mc(n)) = .false.
      itcon_ss(n) = 17
      qcon(itcon_ss(n)) = .true.  ; conpts(5) = 'LS COND'
      qsum(itcon_ss(n)) = .false.
#ifdef TRACERS_DRYDEP
      if(dodrydep(n)) then
        itcon_dd(n,1)=18
        qcon(itcon_dd(n,1)) = .true. ; conpts(6) = 'TURB DEP'
        qsum(itcon_dd(n,1)) = .false.
        itcon_dd(n,2)=19
        qcon(itcon_dd(n,2)) = .true. ; conpts(7) = 'GRAV SET'
        qsum(itcon_dd(n,2)) = .false.
      end if
#endif
      CALL SET_TCON(QCON,TRNAME(N),QSUM,inst_unit(n),
     *     sum_unit(n),scale_inst(n),scale_change(n), N,CONPTs)
      qcon(13:) = .false.  ! reset to defaults for next tracer
      qsum(13:) = .false.  ! reset to defaults for next tracer

#ifdef TRACERS_HETCHEM
      case ('SO4_s1','SO4_s2','SO4_d1', 'SO4_d2', 'SO4_d3', 'SO4_d4')
      itcon_3Dsrc(1,N) = 13
      qcon(itcon_3Dsrc(1,N)) = .true.; conpts(1) = 'Gas phase change'
      qsum(itcon_3Dsrc(1,N)) = .true.
      itcon_mc(n) =14
      qcon(itcon_mc(n)) = .true.  ; conpts(2) = 'MOIST CONV'
      qsum(itcon_mc(n)) = .false.
      itcon_ss(n) =15
      qcon(itcon_ss(n)) = .true.  ; conpts(3) = 'LS COND'
      qsum(itcon_ss(n)) = .false.
#ifdef TRACERS_DRYDEP
      if(dodrydep(n)) then
        itcon_dd(n,1)=16
        qcon(itcon_dd(n,1)) = .true. ; conpts(4) = 'TURB DEP'
        qsum(itcon_dd(n,1)) = .false.
        itcon_dd(n,2)=17
        qcon(itcon_dd(n,2)) = .true. ; conpts(5) = 'GRAV SET'
        qsum(itcon_dd(n,2)) = .false.
      end if
#endif
      CALL SET_TCON(QCON,TRNAME(N),QSUM,inst_unit(n),
     *     sum_unit(n),scale_inst(n),scale_change(n), N,CONPTs)
      qcon(13:) = .false.  ! reset to defaults for next tracer
      qsum(13:) = .false.  ! reset to defaults for next tracer
#endif

      case ('Be7')
      itcon_3Dsrc(1,N) =13
      qcon(itcon_3Dsrc(1,N)) = .true.  ; conpts(1) = 'COSMO SRC'
      qsum(itcon_3Dsrc(1,N)) = .true.
      itcon_mc(n) =14
      qcon(itcon_mc(n)) = .true.  ; conpts(2) = 'MOIST CONV'
      qsum(itcon_mc(n)) = .false.
      itcon_ss(n) =15
      qcon(itcon_ss(n)) = .true.  ; conpts(3) = 'LS COND'
      qsum(itcon_ss(n)) = .false.
#ifdef TRACERS_DRYDEP
      if(dodrydep(n)) then
        itcon_dd(n,1)=16
        qcon(itcon_dd(n,1)) = .true. ; conpts(4) = 'TURB DEP'
        qsum(itcon_dd(n,1)) = .false.
        itcon_dd(n,2)=17
        qcon(itcon_dd(n,2)) = .true. ; conpts(5) = 'GRAV SET'
        qsum(itcon_dd(n,2)) = .false.
      end if
#endif
      itcon_decay(n) = 18
      qcon(itcon_decay(n)) = .true.; conpts(6) = 'DECAY'
      qsum(itcon_decay(n)) = .true.

      CALL SET_TCON(QCON,TRNAME(N),QSUM,inst_unit(n),
     *     sum_unit(n),scale_inst(n),scale_change(n), N,CONPTs)
      qcon(13:) = .false.  ! reset to defaults for next tracer
      qsum(13:) = .false.  ! reset to defaults for next tracer

      case ('Be10')
      itcon_3Dsrc(1,N) =13
      qcon(itcon_3Dsrc(1,N)) = .true.  ; conpts(1) = 'COSMO SRC'
      qsum(itcon_3Dsrc(1,N)) = .true.
      itcon_mc(n) =14
      qcon(itcon_mc(n)) = .true.  ; conpts(2) = 'MOIST CONV'
      qsum(itcon_mc(n)) = .false.
      itcon_ss(n) =15
      qcon(itcon_ss(n)) = .true.  ; conpts(3) = 'LS COND'
      qsum(itcon_ss(n)) = .false.
#ifdef TRACERS_DRYDEP
      if(dodrydep(n)) then
        itcon_dd(n,1)=16
        qcon(itcon_dd(n,1)) = .true. ; conpts(4) = 'TURB DEP'
        qsum(itcon_dd(n,1)) = .false.
        itcon_dd(n,2)=17
        qcon(itcon_dd(n,2)) = .true. ; conpts(5) = 'GRAV SET'
        qsum(itcon_dd(n,2)) = .false.
      end if
#endif

      CALL SET_TCON(QCON,TRNAME(N),QSUM,inst_unit(n),
     *     sum_unit(n),scale_inst(n),scale_change(n), N,CONPTs)
      qcon(13:) = .false.  ! reset to defaults for next tracer
      qsum(13:) = .false.  ! reset to defaults for next tracer

      case ('Pb210')
      itcon_3Dsrc(1,N) =13
      qcon(itcon_3Dsrc(1,N)) = .true.  ; conpts(1) = 'RADIO SRC'
      qsum(itcon_3Dsrc(1,N)) = .true.
      itcon_mc(n) =14
      qcon(itcon_mc(n)) = .true.  ; conpts(2) = 'MOIST CONV'
      qsum(itcon_mc(n)) = .false.
      itcon_ss(n) =15
      qcon(itcon_ss(n)) = .true.  ; conpts(3) = 'LS COND'
      qsum(itcon_ss(n)) = .false.
#ifdef TRACERS_DRYDEP
      if(dodrydep(n)) then
        itcon_dd(n,1)=16
        qcon(itcon_dd(n,1)) = .true. ; conpts(4) = 'TURB DEP'
        qsum(itcon_dd(n,1)) = .false.
        itcon_dd(n,2)=17
        qcon(itcon_dd(n,2)) = .true. ; conpts(5) = 'GRAV SET'
        qsum(itcon_dd(n,2)) = .false.
      end if
#endif
      itcon_decay(n) = 18
      qcon(itcon_decay(n)) = .true.; conpts(6) = 'DECAY'
      qsum(itcon_decay(n)) = .true.

      CALL SET_TCON(QCON,TRNAME(N),QSUM,inst_unit(n),
     *     sum_unit(n),scale_inst(n),scale_change(n), N,CONPTs)
      qcon(13:) = .false.  ! reset to defaults for next tracer
      qsum(13:) = .false.  ! reset to defaults for next tracer

      case ('H2O2_s')
      itcon_3Dsrc(1,N) = 13
      qcon(itcon_3Dsrc(1,N)) = .true.; conpts(1) = 'Gas phase src'
      qsum(itcon_3Dsrc(1,N)) = .true.
      itcon_3Dsrc(2,N) = 14
      qcon(itcon_3Dsrc(2,N)) = .true.; conpts(2) = 'Gas phase sink'
      qsum(itcon_3Dsrc(2,N)) = .true.
      itcon_3Dsrc(3,N) = 15
      qcon(itcon_3Dsrc(3,N)) = .true.; conpts(3) = 'Heter sink'
      qsum(itcon_3Dsrc(3,N)) = .true.
      itcon_mc(n) =16
      qcon(itcon_mc(n)) = .true.  ; conpts(4) = 'MOIST CONV'
      qsum(itcon_mc(n)) = .false.
      itcon_ss(n) =17
      qcon(itcon_ss(n)) = .true.  ; conpts(5) = 'LS COND'
      qsum(itcon_ss(n)) = .false.
#ifdef TRACERS_DRYDEP
      if(dodrydep(n)) then
        itcon_dd(n,1)=18
        qcon(itcon_dd(n,1)) = .true. ; conpts(6) = 'TURB DEP'
        qsum(itcon_dd(n,1)) = .false.
      end if
#endif

      CALL SET_TCON(QCON,TRNAME(N),QSUM,inst_unit(n),
     *     sum_unit(n),scale_inst(n),scale_change(n), N,CONPTs)
      qcon(13:) = .false.  ! reset to defaults for next tracer
      qsum(13:) = .false.  ! reset to defaults for next tracer

      CASE('seasalt1','seasalt2',
     &     'Clay','Silt1','Silt2','Silt3','Silt4',
     &     'ClayIlli','ClayKaol','ClaySmec','ClayCalc','ClayQuar',
     &     'Sil1Quar','Sil1Feld','Sil1Calc','Sil1Hema','Sil1Gyps',
     &     'Sil2Quar','Sil2Feld','Sil2Calc','Sil2Hema','Sil2Gyps',
     &     'Sil3Quar','Sil3Feld','Sil3Calc','Sil3Hema','Sil3Gyps',
     &     'Sil1QuHe','Sil2QuHe','Sil3QuHe')
      itcon_mc(n) =13
      qcon(itcon_mc(n)) = .true.  ; conpts(1) = 'MOIST CONV'
      qsum(itcon_mc(n)) = .false.
      itcon_ss(n) =14
      qcon(itcon_ss(n)) = .true.  ; conpts(2) = 'LS COND'
      qsum(itcon_ss(n)) = .false.
#ifdef TRACERS_DRYDEP
      if(dodrydep(n)) then
        itcon_dd(n,1)=15
        qcon(itcon_dd(n,1)) = .true. ; conpts(3) = 'TURB DEP'
        qsum(itcon_dd(n,1)) = .false.
        itcon_dd(n,2)=16
        qcon(itcon_dd(n,2)) = .true. ; conpts(4) = 'GRAV SET'
        qsum(itcon_dd(n,2)) = .false.
      end if
#endif
#ifndef TRACERS_WATER
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
      itcon_wt(n)=17
      qcon(itcon_wt(n)) = .true. ; conpts(5) = 'WET DEP'
      qsum(itcon_wt(n)) = .false.
#endif
#endif
      CALL SET_TCON(QCON,TRNAME(N),QSUM,inst_unit(n),
     *     sum_unit(n),scale_inst(n),scale_change(n), N,CONPTs)
      qcon(13:) = .false.  ! reset to defaults for next tracer
      qsum(13:) = .false.  ! reset to defaults for next tracer

C**** Here are some more examples of conservation diag configuration
C**** Separate Moist convection/Large scale condensation
c      itcon_mc(n)=xx
c      qcon(itcon_mc(n))=.true.  ; conpts(yy) = 'MOIST CONV'
c      qsum(itcon_mc(n)) = .false.
c      itcon_ss(n)=xx
c      qcon(itcon_ss(n))=.true.  ; conpts(yy) = 'LS COND'
c      qsum(itcon_ss(n)) = .false.

       end select
       end do

C**** print out total tracer diagnostic array size
C
#ifdef TRACERS_DRYDEP
C Read landuse parameters and coefficients for tracer dry deposition:
      CALL RDLAND
      CALL RDDRYCF
#endif
#ifdef TRACERS_SPECIAL_Shindell
      call cheminit ! **** Initialize the chemistry ****
      call special_layers_init
#endif
#endif
#ifdef TRACERS_COSMO
      do n=1,ntm
        if (trname(n) .eq. "Be7" .OR. trname(n) .eq. "Be10") then
          call read_Be_source
          exit
        end if
      end do
#endif
      return
      end subroutine init_tracer


      SUBROUTINE tracer_IC
!@sum tracer_IC initializes tracers when they are first switched on
!@auth Jean Lerner
#ifdef TRACERS_ON
      USE CONSTANT, only: mair,rhow
#ifdef TRACERS_AEROSOLS_Koch
     *,sday
#endif
      USE MODEL_COM, only: itime,im,jm,lm,ls1,jday,JEQ,ptop,psf,sig
#ifdef TRACERS_AEROSOLS_Koch
     * ,dtsrc
#endif
#ifdef TRACERS_WATER
     *  ,q,wm,flice,fearth
#endif
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
     &     ,JMperY,JDperY
#endif
      USE DOMAIN_DECOMP, only : GRID, GET
      USE SOMTQ_COM, only : qmom,mz,mzz
      USE TRACER_COM, only: ntm,trm,trmom,itime_tr0,trname,needtrs,
     *   tr_mm
#ifdef TRACERS_AEROSOLS_Koch
     *   ,imAER,n_SO2,imPI
#endif
#ifdef regional_Ox_tracers
     *   ,NregOx,n_Ox
#endif
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
     &   ,imDUST,kjm,table1,x11,x21
#endif
#ifdef TRACERS_WATER
     *  ,trwm,trw0,tr_wd_TYPE,nWATER
      USE LANDICE, only : ace1li,ace2li
      USE LANDICE_COM, only : trli0,trsnowli,trlndi,snowli
      USE SEAICE, only : xsi,ace1i
      USE SEAICE_COM, only : rsi,msi,snowi,trsi,trsi0,ssi
      USE LAKES_COM, only : trlake,mwl,mldlk,flake
      USE GHY_COM, only : tr_wbare,tr_wvege,tr_wsn_ij,wbare,wvege
     &     ,wsn_ij,nsn_ij,fr_snow_ij
      USE FLUXES, only : gtracer
#endif
      USE GEOM, only: dxyp,bydxyp,lat_dg
      USE DYNAMICS, only: am,byam  ! Air mass of each box (kg/m^2)
      USE PBLCOM, only: npbl,trabl
#ifdef TRACERS_WATER
     *  ,qabl
#endif
#ifdef TRACERS_SPECIAL_Lerner
      USE LINOZ_CHEM_COM, only: tlt0m,tltzm, tltzzm
      USE PRATHER_CHEM_COM, only: nstrtc
#endif
      USE FILEMANAGER, only: openunit,closeunit
#ifdef TRACERS_SPECIAL_Shindell
      USE RAD_COM, ONLY : O3_tracer_save
      USE TRCHEM_Shindell_COM,only:O3MULT,COlat,MDOFM
     &  ,COalt,JCOlat,OxIC,byO3MULT,pfix_CH4_S,pfix_CH4_N,PI_run,
     &  fix_CH4_chemistry,PIratio_N,PIratio_CO_T,PIratio_CO_S,
     &  PIratio_other
#ifdef SHINDELL_STRAT_CHEM
     &  ,ClOxalt,BrOxalt,ClONO2alt,HClalt,N2OICX,CFCIC
#endif
#endif
#ifdef TRACERS_AEROSOLS_Koch
      USE AEROSOL_SOURCES, only: DMSinput,BCI_src,OCI_src,
     * BCB_src,OCB_src,OCT_src
     * ,DMS_AER,SS1_AER,SS2_AER
     * ,SO2_src_3D,SO2_biosrc_3D,SO2_src
#endif
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
      USE tracers_dust,ONLY : hbaij,ricntd,dryhr,frclay,frsilt,vtrsh,
     &     ers_data,gin_data,table,x1,x2,x3,lim,ljm,lkm,d_dust
#if (defined TRACERS_MINERALS) || (defined TRACERS_QUARZHEM)
     &     ,Mtrac,minfr
#endif
#endif
      IMPLICIT NONE
      INTEGER i,n,l,j,iu_data,ipbl,it,lr
      CHARACTER*80 title
      REAL*8 CFC11ic,conv
      REAL*8 :: trinit =1., tmominit=0.

      REAL*8, DIMENSION(im,GRID%J_STRT_HALO:GRID%J_STOP_HALO,lm) :: 
     *                                                    ic14CO2 
      REAL*4, DIMENSION(im,jm,lm) ::  CO2ic   !each proc. reads global array
      REAL*4, DIMENSION(jm,lm)    ::  N2Oic   !each proc. reads global array
      REAL*8, DIMENSION(GRID%J_STRT_HALO:GRID%J_STOP_HALO,lm) :: 
     *                                                      CH4ic

!@param bymair 1/molecular wt. of air = 1/mair
!@param byjm 1./JM
      REAL*8, PARAMETER :: bymair = 1.d0/mair, byjm =1.d0/JM
#ifdef TRACERS_SPECIAL_Shindell
!@var imonth dummy index for choosing the right month
!@var ICfactor varying factor for altering initial conditions
      INTEGER imonth, J2
      REAL*8 ICfactor
!@var PRES local nominal pressure for vertical interpolations
      REAL*8, DIMENSION(LM) :: PRES
#ifdef regional_Ox_tracers
!@var byNregOx reciprocal of the number of regional Ox tracers
      REAL*8 byNregOx
#endif
#endif
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_DUST) ||\
    (defined TRACERS_MINERALS) || (defined TRACERS_QUARZHEM)
      include 'netcdf.inc'
#endif
#ifdef TRACERS_AEROSOLS_Koch
      integer start(3),count(3),status,ncidu,id1
      REAL*8 dmsconc
      INTEGER mon_unit, mont,ii,jj,ir,mm,iuc,m,mmm,ll
      real*8 carbstuff,ccnv
#endif
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
      INTEGER :: io_data,k
      INTEGER startd(3),countd(3),statusd
      INTEGER idd1,idd2,idd3,idd4,ncidd1,ncidd2,ncidd3,ncidd4
      REAL*4 :: work(Im,Jm)
      REAL*8 :: sum
      LOGICAL,SAVE :: ifirst=.TRUE.
#endif

      INTEGER J_0, J_1
      INTEGER J_0H, J_1H
      LOGICAL HAVE_SOUTH_POLE, HAVE_NORTH_POLE

C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT=J_0,       J_STOP=J_1, 
     *               J_STRT_HALO=J_0H, J_STOP_HALO=J_1H,
     *               HAVE_SOUTH_POLE=HAVE_SOUTH_POLE,
     *               HAVE_NORTH_POLE=HAVE_NORTH_POLE)

#ifdef regional_Ox_tracers
      byNregOx=1.d0/float(NregOx)
#endif
#ifdef SHINDELL_STRAT_CHEM
      PRES(:)=SIG(:)*(PSF-PTOP)+PTOP
#endif
      do n=1,ntm
      if (itime.eq.itime_tr0(n)) then

C**** set some defaults for air mass tracers
      trm(:,J_0:J_1,:,n) = 0.
      trmom(:,:,J_0:J_1,:,n) = 0.

#ifdef TRACERS_WATER
C**** set some defaults for water tracers
      trwm(:,J_0:J_1,:,n)=0.
      trlake(n,:,:,J_0:J_1)=0.
      trsi(n,:,:,J_0:J_1)=0.
      trlndi(n,:,J_0:J_1)=0.
      trsnowli(n,:,J_0:J_1)=0.
      tr_wbare(n,:,:,J_0:J_1)=0.
      tr_wvege(n,:,:,J_0:J_1)=0.
      tr_wsn_ij(n,:,:,:,J_0:J_1)=0.
#endif
      select case (trname(n))

        case default
c          write(6,*) 'In TRACER_IC:',trname(n),' does not exist '
          call stop_model("TRACER_IC",255)

        case ('Air')
          do l=1,lm
          do j=J_0,J_1
            trm(:,j,l,n) = am(l,:,j)*dxyp(j)
          end do; enddo
          if (HAVE_SOUTH_POLE) then
             do i=2,im
                trm(i,1,:,n) =  trm(1,1,:,n) !poles
             enddo
          endif
          if (HAVE_NORTH_POLE) then
             do i=2,im
                trm(i,jm,:,n) = trm(1,jm,:,n) !poles
             enddo
          endif

        case ('SF6')
          ! defaults ok

        case ('Be7', 'Be10', 'Pb210', 'Rn222')
          ! defaults ok

        case ('CO2')
          call openunit('CO2_IC',iu_data,.true.,.true.)
C**** ESMF: Each processor reads the global array co2ic
          read (iu_data) title,co2ic
          call closeunit(iu_data)
          write(6,*) title,' read from CO2_IC'
          do l=1,lm         !ppmv==>ppmm
          do j=J_0,J_1
            trm(:,j,l,n) = co2ic(:,j,l)*am(l,:,j)*dxyp(j)*1.54d-6
          enddo; enddo

        case ('N2O')
#ifdef TRACERS_SPECIAL_Lerner
          call openunit('N2O_IC',iu_data,.true.,.true.)
C**** ESMF: Each processor reads the global array: N2Oic
          read (iu_data) title,N2Oic     ! unit is PPMM/(M*DXYP)
          call closeunit(iu_data)
          write(6,*) title,' read from N2O_IC'
          do l=1,lm         !ppmv==>ppmm
          do j=J_0,J_1
            trm(:,j,l,n) = am(l,:,j)*dxyp(j)*N2Oic(j,l)
          enddo; enddo
#endif
#ifdef SHINDELL_STRAT_CHEM
          do l=1,lm; do j=J_0,J_1; do i=1,im
            trm(i,j,l,n) = N2OICX(i,j,l)
          end do   ; end do   ; end do
#endif

        case ('CFC11')   !!! should start April 1
          CFC11ic = 268.D-12*136.5/29.029    !268 PPTV
          do l=1,lm
          do j=J_0,J_1
            trm(:,j,l,n) = am(l,:,j)*dxyp(j)*CFC11ic
          enddo; enddo

        case ('14CO2')   !!! this tracer is supposed to start 10/16
#ifdef TRACERS_SPECIAL_Lerner
          call get_14CO2_IC(ic14CO2)
          do l=1,lm         !ppmv==>ppmm
          do j=J_0,J_1
            trm(:,j,l,n) = am(l,:,j)*dxyp(j)*ic14CO2(:,j,l)*1.d-18
          enddo; enddo
#endif

        case ('CH4')
#ifdef TRACERS_SPECIAL_Shindell
          select case(fix_CH4_chemistry)
          case(1)
            do l=1,lm; do j=J_0,J_1/2
              trm(:,j,l,n)=am(l,:,j)*dxyp(j)*TR_MM(n)*bymair*pfix_CH4_S
            end do; end do
            do l=1,lm; do j=(J_1/2)+1,J_1
              trm(:,j,l,n)=am(l,:,j)*dxyp(j)*TR_MM(n)*bymair*pfix_CH4_N
            end do; end do
          case default
            call get_CH4_IC ! defines trm(:,:,:,n_CH4) within
          end select
#endif
#ifdef TRACERS_SPECIAL_Lerner
          call get_wofsy_gas_IC(trname(n),CH4ic)
          do l=1,lm         !ppbv==>ppbm
          do j=J_0,J_1
            trm(:,j,l,n) = am(l,:,j)*dxyp(j)*CH4ic(j,l)*0.552d-9
          enddo; enddo
#endif

        case ('O3')
          do l=1,lm
          do j=J_0,J_1
            trm(:,j,l,n) = am(l,:,j)*dxyp(j)*20.d-9*tr_mm(n)/mair
          enddo; enddo
#ifdef TRACERS_SPECIAL_Lerner
          do l=lm,lm+1-nstrtc,-1
          lr = lm+1-l
            do j=J_0,J_1
            if (tlt0m(j,lr,5) /= 0.) then
            trm(:,j,l,n) =
     *          tlt0m(j,lr,1)*am(l,:,j)*dxyp(j)*tr_mm(n)/mair
            trmom(mz,:,j,l,n)  =
     *          tltzm(j,lr,1)*am(l,:,j)*dxyp(j)*tr_mm(n)/mair
            trmom(mzz,:,j,l,n)  =
     *         tltzzm(j,lr,1)*am(l,:,j)*dxyp(j)*tr_mm(n)/mair
            end if
            end do
          end do
#endif

#ifdef TRACERS_WATER
      case ('Water', 'H2O18', 'HDO', 'HTO')

C**** initial atmospheric conc. needs to be defined for each tracer
        select case (trname(n))
        case ('Water')
          trinit=1.
C**** for gradients defined on air mass
          tmominit = 1.
C**** for gradients defined on water mass (should be an option?)
c     tmominit = 0.
        case ('H2O18')        ! d18O=-80
          trinit=0.92d0*trw0(n)
          tmominit = trinit
        case ('HDO')   ! dD=-630
          trinit=0.37d0*trw0(n)
          tmominit = trinit
        case ('HTO')
          trinit=0.
          tmominit = trinit
        end select

        do l=1,lm
        do j=J_0,J_1
          trm(:,j,l,n) =  q(:,j,l)*am(l,:,j)*dxyp(j)*trinit
          trwm(:,j,l,n)= wm(:,j,l)*am(l,:,j)*dxyp(j)*trinit
          do i=1,im
            trmom(:,i,j,l,n) = qmom(:,i,j,l)*am(l,i,j)*dxyp(j)*tmominit
          end do
        end do
        end do
        if (HAVE_SOUTH_POLE) then
           do i=2,im
              trm(i,1,:,n) =  trm(1,1,:,n) !poles
              trwm(i, 1,:,n)= trwm(1, 1,:,n) !poles
              trmom(:,i, 1,:,n)=0.
           enddo
        endif
        if (HAVE_NORTH_POLE) then
           do i=2,im
              trm(i,jm,:,n) = trm(1,jm,:,n) !poles
              trwm(i,jm,:,n)= trwm(1,jm,:,n) !poles
              trmom(:,i,jm,:,n)=0.
          enddo
        endif
        if (trname(n).eq."HTO") then ! initialise bomb source
          do l=ls1-1,ls1+1      ! strat. source
            do j=J_0,J_1   ! lat 44 N - 56 N
              if (nint(lat_dg(j,1)).ge.44 .and. nint(lat_dg(j,1)).le.56)
     *             trm(:,j,l,n)= q(:,j,l)*am(l,:,j)*dxyp(j)*1d10*1d-18
            end do
          end do
        end if

        do j=J_0,J_1
          do i=1,im
C**** lakes
            if (flake(i,j).gt.0) then
              trlake(n,1,i,j)=trw0(n)*mldlk(i,j)*rhow*flake(i,j)*dxyp(j)
              if (mwl(i,j)-mldlk(i,j)*rhow*flake(i,j)*dxyp(j).gt.1d-10
     *             *mwl(i,j)) then
                trlake(n,2,i,j)=trw0(n)*mwl(i,j)-trlake(n,1,i,j)
              else
                trlake(n,2,i,j)=0.
              end if
              gtracer(n,1,i,j)=trw0(n)
            elseif (fearth(i,j).gt.0) then
              trlake(n,1,i,j)=trw0(n)*mwl(i,j)
              trlake(n,2,i,j)=0.
            else
              trlake(n,1:2,i,j)=0.
            end if
c**** ice
            if (msi(i,j).gt.0) then
              trsi(n,1,i,j)=trsi0(n)*
     *             (xsi(1)*(snowi(i,j)+ace1i)-ssi(1,i,j))
              trsi(n,2,i,j)=trsi0(n)*
     *             (xsi(2)*(snowi(i,j)+ace1i)-ssi(2,i,j))
              trsi(n,3,i,j)=trsi0(n)*(xsi(3)*msi(i,j)-ssi(3,i,j))
              trsi(n,4,i,j)=trsi0(n)*(xsi(4)*msi(i,j)-ssi(4,i,j))
              gtracer(n,2,i,j)=trsi0(n)
            end if
c**** landice
            if (flice(i,j).gt.0) then
              trlndi(n,i,j)=trli0(n)*(ace1li+ace2li)
              trsnowli(n,i,j)=trli0(n)*snowli(i,j)
              gtracer(n,3,i,j)=trli0(n)
            else
              trlndi(n,i,j)=0.
              trsnowli(n,i,j)=0.
              gtracer(n,3,i,j)=0.
            end if
c**** earth
            if (fearth(i,j).gt.0) then
              conv=rhow         ! convert from m to kg/m^2
              tr_wbare  (n,:,i,j)=trw0(n)*wbare (:,i,j)*conv
              tr_wvege  (n,:,i,j)=trw0(n)*wvege (:,i,j)*conv
              tr_wsn_ij(n,1:nsn_ij(1,i,j),1,i,j)=
     &             trw0(n)*wsn_ij(1:nsn_ij(1,i,j),1,i,j)
     &             *fr_snow_ij(1,i,j)*conv
              tr_wsn_ij(n,1:nsn_ij(2,i,j),2,i,j)=
     &             trw0(n)*wsn_ij(1:nsn_ij(2,i,j),2,i,j)
     &             *fr_snow_ij(2,i,j)*conv
              !trsnowbv(n,2,i,j)=trw0(n)*snowbv(2,i,j)*conv
              gtracer (n,4,i,j)=trw0(n)
            else
              tr_wbare  (n,:,i,j)=0.
              tr_wvege  (n,:,i,j)=0.
              tr_wsn_ij(n,:,:,i,j)=0.
              !trsnowbv(n,1,i,j)=0.
              !trsnowbv(n,2,i,j)=0.
              gtracer(n,4,i,j)=0.
            end if
          end do
          end do
#endif

#ifdef TRACERS_SPECIAL_Shindell
        case ('Ox')
          do l=1,lm; do j=J_0,J_1; do i=1,im
            trm(I,J,L,n) = OxIC(I,J,L)
            O3_tracer_save(L,I,J)=OxIC(I,J,L)*byO3MULT*bydxyp(j)
          end do   ; end do   ; end do
#endif
#ifdef regional_Ox_tracers
        case ('OxREG1','OxREG2','OxREG3','OxREG4','OxREG5','OxREG6')
          trm(:,J_0:J_1,:,n)=trm(:,J_0:J_1,:,n_Ox)*byNregOx
#endif

        case ('NOx')
#ifdef TRACERS_SPECIAL_Shindell
          select case(PI_run)
          case(1)     ; ICfactor=PIratio_N
          case default; ICfactor=1.d0
          end select
          do l=1,lm; do j=J_0,J_1; do i=1,im
            trm(i,j,l,n) = am(l,i,j)*dxyp(j)*1.d-11*ICfactor
#ifdef SHINDELL_STRAT_CHEM
            if(PRES(L).lt.10.)trm(i,j,l,n)=trm(i,j,l,n)*3.d2
#endif
          end do; end do; end do
#endif

#if (defined TRACERS_SPECIAL_Shindell) && (defined SHINDELL_STRAT_CHEM)
        case ('ClOx')
          do l=1,lm; do j=J_0,J_1; do i=1,im
            trm(i,j,l,n) =
     &      am(l,i,j)*dxyp(j)*TR_MM(n)*bymair*1.d-11*ClOxalt(l)
          end do; end do; end do

        case ('BrOx')
          do l=1,lm; do j=J_0,J_1; do i=1,im
            trm(i,j,l,n) = 
     &      am(l,i,j)*dxyp(j)*TR_MM(n)*bymair*1.d-11*BrOxalt(l)
          end do; end do; end do

        case ('HCl')
          do l=1,lm; do j=J_0,J_1; do i=1,im
            trm(i,j,l,n) = 
     &      am(l,i,j)*dxyp(j)*TR_MM(n)*bymair*1.d-11*HClalt(l)
          end do; end do; end do

        case ('ClONO2')
          do l=1,lm; do j=J_0,J_1; do i=1,im
            trm(i,j,l,n) = 
     &      am(l,i,j)*dxyp(j)*TR_MM(n)*bymair*1.d-11*ClONO2alt(l)
          end do; end do; end do
#endif
        case ('N2O5')
#ifdef TRACERS_SPECIAL_Shindell
          select case(PI_run)
          case(1)     ; ICfactor=PIratio_N
          case default; ICfactor=1.d0
          end select
          do l=1,lm; do j=J_0,J_1; do i=1,im
            trm(i,j,l,n) = am(l,i,j)*dxyp(j)*1.d-12*ICfactor
          end do; end do; end do
#endif

        case ('HNO3')
#ifdef TRACERS_SPECIAL_Shindell
          select case(PI_run)
          case(1)     ; ICfactor=PIratio_N
          case default; ICfactor=1.d0
          end select
          do l=1,lm; do j=J_0,J_1; do i=1,im
            trm(i,j,l,n) = am(l,i,j)*dxyp(j)*1.d-10*ICfactor
#ifdef SHINDELL_STRAT_CHEM
            if(PRES(L).lt.50.and.PRES(L).gt.10.)
     &      trm(i,j,l,n)=trm(i,j,l,n)*1.d2
#endif
          end do; end do; end do
#endif

        case ('H2O2')
          do l=1,lm; do j=J_0,J_1; do i=1,im
            trm(i,j,l,n) = am(l,i,j)*dxyp(j)*5.d-10
          end do; end do; end do

        case ('CH3OOH')
          do l=1,lm; do j=J_0,J_1; do i=1,im
            trm(i,j,l,n) = am(l,i,j)*dxyp(j)*1.d-11
          end do; end do; end do

        case ('HCHO')
          do l=1,lm; do j=J_0,J_1; do i=1,im
            trm(i,j,l,n) = am(l,i,j)*dxyp(j)*1.d-11
          end do; end do; end do

        case ('HO2NO2')
#ifdef TRACERS_SPECIAL_Shindell
          select case(PI_run)
          case(1)     ; ICfactor=PIratio_N
          case default; ICfactor=1.d0
          end select
          do l=1,lm; do j=J_0,J_1; do i=1,im
            trm(i,j,l,n) = am(l,i,j)*dxyp(j)*1.d-12*ICfactor
          end do; end do; end do

        case('CO')
C         COlat=ppbv, COalt=no unit, TR_MM(n)*bymair=ratio of mol.wt.,
C         AM=kg/m2, and DXYP=m2:
          DO L=1,LM
            select case(PI_run)
            case(1) ! pre-industrial
              if(L.le.LS1-1) then
                ICfactor=PIratio_CO_T ! troposphere
              else
                ICfactor=PIratio_CO_S ! stratosphere
              end if
            case default; ICfactor=1.d0
            end select
            DO J=J_0,J_1
              J2=MAX(1,NINT(float(J)*float(JCOlat)*BYJM))
              DO I=1,IM
                trm(i,j,l,n)=COlat(J2)*COalt(L)*1.D-9*TR_MM(n)*bymair*
     &          am(L,I,J)*DXYP(J)*ICfactor
              END DO
            END DO
          END DO
          J2=0

        case ('PAN')
          select case(PI_run)
          case(1)     ; ICfactor=PIratio_other
          case default; ICfactor=1.d0
          end select
          do l=1,lm; do j=J_0,J_1; do i=1,im
            trm(i,j,l,n) = 
     &      am(l,i,j)*dxyp(j)*TR_MM(n)*bymair*4.d-11*ICfactor
          end do; end do; end do

        case ('Isoprene')
          select case(PI_run)
          case(1)     ; ICfactor=PIratio_other
          case default; ICfactor=1.d0
          end select
          do l=1,lm; do j=J_0,J_1; do i=1,im
            trm(i,j,l,n) = 
     &      am(l,i,j)*dxyp(j)*TR_MM(n)*bymair*1.d-11*ICfactor
          end do; end do; end do

        case ('AlkylNit')
          select case(PI_run)
          case(1)     ; ICfactor=PIratio_other
          case default; ICfactor=1.d0
          end select
          do l=1,lm; do j=J_0,J_1; do i=1,im
            trm(i,j,l,n) = 
     &      am(l,i,j)*dxyp(j)*TR_MM(n)*bymair*2.d-10*ICfactor
          end do; end do; end do

        case('Alkenes')
          select case(PI_run)
          case(1)     ; ICfactor=PIratio_other
          case default; ICfactor=1.d0
          end select
          do l=1,lm; do j=J_0,J_1; do i=1,im
            trm(i,j,l,n) = 
     &      am(l,i,j)*dxyp(j)*TR_MM(n)*bymair*4.d-10*ICfactor
          end do; end do; end do

        case('Paraffin')
          select case(PI_run)
          case(1)     ; ICfactor=PIratio_other
          case default; ICfactor=1.d0
          end select
          do l=1,lm; do j=J_0,J_1; do i=1,im
            trm(i,j,l,n) = 
     &      am(l,i,j)*dxyp(j)*TR_MM(n)*bymair*5.d-10*ICfactor
          end do; end do; end do
#endif

        case ('CFC')
#ifdef SHINDELL_STRAT_CHEM
          do l=1,lm; do j=J_0,J_1; do i=1,im
            trm(I,J,L,n) = CFCIC(I,J,L)*4.25d0 ! note factor
          end do   ; end do   ; end do
#endif

        case ('BrONO2','HBr','HOBr')
          do l=1,lm; do j=J_0,J_1; do i=1,im
            if(l.ge.LS1) then
              trm(i,j,l,n) = am(l,i,j)*dxyp(j)*TR_MM(n)*bymair*2.d-13
            else
              trm(i,j,l,n) = am(l,i,j)*dxyp(j)*TR_MM(n)*bymair*1.d-13
            end if
          end do; end do; end do

        case ('HOCl')
          do l=1,lm; do j=J_0,J_1; do i=1,im
            if(l.ge.LS1) then
              trm(i,j,l,n) = am(l,i,j)*dxyp(j)*TR_MM(n)*bymair*5.d-11
            else
              trm(i,j,l,n) = am(l,i,j)*dxyp(j)*TR_MM(n)*bymair*1.d-11
            end if
          end do; end do; end do

        case('DMS')
          do l=1,lm; do j=J_0,J_1; do i=1,im
            trm(i,j,l,n) = am(l,i,j)*dxyp(j)*TR_MM(n)*bymair*5.d-13
          end do; end do; end do

        case('MSA', 'SO2', 'SO4', 'SO4_d1', 'SO4_d2', 'SO4_d3','SO4_d4',
     *         'BCII', 'BCIA', 'BCB', 'OCII', 'OCIA', 'OCB', 'H2O2_s',
     *         'seasalt1', 'seasalt2','SO4_s1','SO4_s2')
          do l=1,lm; do j=J_0,J_1; do i=1,im
            trm(i,j,l,n) = am(l,i,j)*dxyp(j)*TR_MM(n)*bymair*5.d-14
          end do; end do; end do

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
        CASE('Clay','Silt1','Silt2','Silt3','Silt4',
     &       'ClayIlli','ClayKaol','ClaySmec','ClayCalc','ClayQuar',
     &       'Sil1Quar','Sil1Feld','Sil1Calc','Sil1Hema','Sil1Gyps',
     &       'Sil2Quar','Sil2Feld','Sil2Calc','Sil2Hema','Sil2Gyps',
     &       'Sil3Quar','Sil3Feld','Sil3Calc','Sil3Hema','Sil3Gyps',
     &       'Sil1QuHe','Sil2QuHe','Sil3QuHe')
          ! defaults ok
          hbaij=0D0
          ricntd=0D0
#endif

      end select

C**** Initialise pbl profile if necessary
      if (needtrs(n)) then
        do it=1,4
        do j=J_0,J_1
        do ipbl=1,npbl
#ifdef TRACERS_WATER
          if(tr_wd_TYPE(n).eq.nWATER)THEN
            trabl(ipbl,n,:,j,it) = trinit*qabl(ipbl,:,j,it)
          ELSE
#endif
            trabl(ipbl,n,:,j,it) = trm(:,j,1,n)*byam(1,:,j)*bydxyp(j)
#ifdef TRACERS_WATER
          END IF
#endif
        end do
        end do
        end do
      end if

      write(6,*) ' Tracer ',trname(n),' initialized at itime=',itime
      end if
      end do
#endif

#if (defined TRACERS_WATER) || (defined TRACERS_OCEAN)
C**** Initialise ocean tracers if necessary
      call tracer_ic_ocean
#endif
C****
#ifdef TRACERS_AEROSOLS_Koch
c read in DMS source
      if (imAER.eq.0) then
      call openunit('DMS_SEA',mon_unit,.false.)
      DO 8 mm =1,12
      READ(mon_unit,*) mont
      do ir=1,9999
      read(mon_unit,901) ii,jj,dmsconc
c I'm not using the moments
      IF (II.EQ.0) GO TO 8
      DMSinput(ii,jj,mm)=dmsconc
      end do
  8   CONTINUE
        call closeunit(mon_unit)
      else  !AEROCOM run for now these are in
         status=NF_OPEN('DMS_SEA',NCNOWRIT,ncidu)
         status=NF_INQ_VARID(ncidu,'dms',id1)
          start(1)=1
          start(2)=1
          start(3)=1
          count(1)=im
          count(2)=jm
          count(3)=366
          status=NF_GET_VARA_REAL(ncidu,id1,start,count,DMS_AER)
          status=NF_CLOSE('DMS_SEA',NCNOWRIT,ncidu)
      endif
 901  FORMAT(3X,2(I4),E11.3,5F9.2)
c read in AEROCOM seasalt
      if (imAER.eq.1) then
         status=NF_OPEN('SALT1',NCNOWRIT,ncidu)
         status=NF_INQ_VARID(ncidu,'salt',id1)
          start(1)=1
          start(2)=1
          start(3)=1
          count(1)=im
          count(2)=jm
          count(3)=366
          status=NF_GET_VARA_REAL(ncidu,id1,start,count,SS1_AER)
          status=NF_CLOSE('SALT1',NCNOWRIT,ncidu)
         status=NF_OPEN('SALT2',NCNOWRIT,ncidu)
         status=NF_INQ_VARID(ncidu,'salt',id1)
          start(1)=1
          start(2)=1
          start(3)=1
          count(1)=im
          count(2)=jm
          count(3)=366
          status=NF_GET_VARA_REAL(ncidu,id1,start,count,SS2_AER)
          status=NF_CLOSE('SALT2',NCNOWRIT,ncidu)
      endif
c read in AEROCOM SO2 emissions
c Industrial
      SO2_src(:,:,:)=0.d0
      SO2_src_3D(:,:,:,:)=0.d0
      if (imPI.eq.0) then
      call openunit('SO2_IND',iuc,.false.)
      do mm=1,9999
      read(iuc,*) ii,jj,carbstuff
      if (ii.eq.0.) go to 81
      SO2_src(ii,jj,1)=carbstuff/(sday*30.4)/12.d0
      end do
  81  call closeunit(iuc)
      endif
c volcano - continuous 
      call openunit('SO2_VOLCANO',iuc,.false.)
      do mm=1,99999
      read(iuc,*) ii,jj,ll,carbstuff
      if (ii.eq.0.) go to 83
      SO2_src_3D(ii,jj,ll,1)=carbstuff/(sday*30.4d0)/12.d0
      end do
  83  call closeunit(iuc)
      if (imAER.eq.1) then
      SO2_biosrc_3D(:,:,:,:)=0.d0
c Biomass 
      call openunit('SO2_BIOMASS',iuc,.false.)
      do mm=1,99999
      read(iuc,*) ii,jj,mmm,ll,carbstuff
      if (ii.eq.0.) go to 82
      if (imPI.eq.1) carbstuff=carbstuff*0.5d0
      SO2_biosrc_3D(ii,jj,ll,mmm)=carbstuff/(sday*30.4d0)
      end do
  82  call closeunit(iuc)
c volcano - explosive 
      call openunit('SO2_VOLCANO_EXP',iuc,.false.)
      do mm=1,99999
      read(iuc,*) ii,jj,ll,carbstuff
      if (ii.eq.0.) go to 84
      SO2_src_3D(ii,jj,ll,1)=SO2_src_3D(ii,jj,ll,1)
     * +carbstuff/(sday*30.4d0)/12.d0
      end do
  84  call closeunit(iuc)
      
      endif 
c read in BC and OC sources
c Terpenes first
      OCT_src(:,:,:)=0.d0
      if (imAER.eq.0) then
      call openunit('TERPENE',mon_unit,.false.)
      do mm=1,12
      do mmm=1,9999
      read(mon_unit,*) ii,jj,carbstuff
      if (ii.eq.0.) go to 15
      OCT_src(ii,jj,mm)=carbstuff
      end do
 15   continue
      end do
      call closeunit(mon_unit)
c units are mg Terpene/m2/month
      do mm=1,12; do j=J_0,J_1; do i=1,im
       OCT_src(i,j,mm) = OCT_src(i,j,mm)*dxyp(j)*0.1d0
     *  *1.3d0/(1.D6*3600.d0*24.d0*30.4d0)
      end do; end do; end do
c This assumes 10% emission yield (Chin, Penner)
c 1.3 converts OC to OM
      else ! AEROCOM
      call openunit('TERPENE',mon_unit,.false.)
      do mmm=1,99999
      read(mon_unit,*) ii,jj,mm,carbstuff
      if (ii.eq.0.) go to 18
      carbstuff=carbstuff/(sday*30.4d0)
      OCT_src(ii,jj,mm)=carbstuff*1.3d0
      end do
 18   continue
      call closeunit(mon_unit)
      endif
c Now industrial and biomass
      BCI_src(:,:)=0.d0
      if (imPI.eq.0) then
      call openunit('BC_BIOFUEL',iuc,.false.)
      do mm=1,9999
      read(iuc,*) ii,jj,carbstuff
      if (ii.eq.0.) go to 9
      BCI_src(ii,jj)=carbstuff
      end do
  9   call closeunit(iuc)
      call openunit('BC_FOSSIL_FUEL',iuc,.false.)
      do mm=1,9999
      read(iuc,*) ii,jj,carbstuff
      if (ii.eq.0.) go to 10
      BCI_src(ii,jj)=BCI_src(ii,jj)+carbstuff
      end do
 10   call closeunit(iuc)
      OCI_src(:,:)=0.d0
      call openunit('OC_BIOFUEL',iuc,.false.)
      do mm=1,9999
      read(iuc,*) ii,jj,carbstuff
      if (ii.eq.0.) go to 11
      OCI_src(ii,jj)=carbstuff
      end do
 11   call closeunit(iuc)
      call openunit('OC_FOSSIL_FUEL',iuc,.false.)
      do mm=1,9999
      read(iuc,*) ii,jj,carbstuff
      if (ii.eq.0.) go to 12
      OCI_src(ii,jj)=OCI_src(ii,jj)+carbstuff
      end do
 12   call closeunit(iuc)
      do i=1,im
      do j=j_0,J_1
      ccnv=1.d0/(sday*30.4)  !*dxyp(j))
c convert from month to second. dxyp??
      BCI_src(i,j)=BCI_src(i,j)*ccnv/12.d0
      OCI_src(i,j)=OCI_src(i,j)*ccnv/12.d0*1.3d0
      end do
      end do
      endif
      BCB_src(:,:,:,:)=0.d0
      OCB_src(:,:,:,:)=0.d0
      if (imAER.eq.0) then
      call openunit('BC_BIOMASS',iuc,.false.)
      do mm=1,12
      do mmm=1,9999
      read(iuc,*) ii,jj,carbstuff
      if (ii.eq.0.) go to 71
      if (imPI.eq.1) carbstuff=carbstuff*0.5d0
      carbstuff=carbstuff*1000.d0/30.4d0/sday
      BCB_src(ii,jj,1,mm)=carbstuff
      OCB_src(ii,jj,1,mm)=carbstuff*7.9d0
      end do
 71   continue
      end do
      call closeunit(iuc)
      else ! AEROCOM
      call openunit('BC_BIOMASS',iuc,.false.)
      do mm=1,99999
      read(iuc,*) ii,jj,mmm,ll,carbstuff 
      if (ii.eq.0.) go to 13
      if (imPI.eq.1) carbstuff=carbstuff*0.5d0
      BCB_src(ii,jj,ll,mmm)=carbstuff
      end do
 13   call closeunit(iuc)
      call openunit('OC_BIOMASS',iuc,.false.)
      do mm=1,99999
      read(iuc,*) ii,jj,mmm,ll,carbstuff 
      if (ii.eq.0.) go to 14
      if (imPI.eq.1) carbstuff=carbstuff*0.5d0
      OCB_src(ii,jj,ll,mmm)=carbstuff
      end do
 14   call closeunit(iuc)
      do i=1,im
      do j=j_0,J_1
      ccnv=1.d0/(sday*30.4) 
      do m=1,12 
      do l=1,7
      BCB_src(i,j,l,m)=BCB_src(i,j,l,m)*ccnv
      OCB_src(i,j,l,m)=OCB_src(i,j,l,m)*ccnv*1.3d0
      end do
      end do
      end do
      end do
      endif
#endif
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
      IF (ifirst) THEN

c     read in lookup table for calculation of mean surface wind speed from PDF
        CALL openunit('LKTAB1',io_data,.true.,.true.)
        READ (io_data) table1
        CALL closeunit(io_data)
c     index of table for sub grid scale velocity (sigma) from 0.0001 to 50 m/s
        sum=0.d0
        DO j=1,kjm
          IF (j <= 30) THEN
            sum=sum+0.0001d0+float(j-1)*0.00008d0
            x21(j)=sum
          ELSE IF (j > 30) THEN
            sum=sum-0.055254d0+0.005471d0*float(j)-
     &           1.938365d-4*float(j)**2.d0+
     &           3.109634d-6*float(j)**3.d0-
     &           2.126684d-8*float(j)**4.d0+
     &           5.128648d-11*float(j)**5.d0
            x21(j)=sum
          END IF
        END DO
c     index of table for GCM surface wind speed from 0.0001 to 50 m/s 
        x11(:)=x21(:)
c     prescribed AEROCOM dust emissions
        IF (imDust == 1) THEN
          statusd=NF_OPEN('dust_bin1',NCNOWRIT,ncidd1)
          statusd=NF_OPEN('dust_bin2',NCNOWRIT,ncidd2)
          statusd=NF_OPEN('dust_bin3',NCNOWRIT,ncidd3)
          statusd=NF_OPEN('dust_bin4',NCNOWRIT,ncidd4)

          statusd=NF_INQ_VARID(ncidd1,'dust',idd1)
          statusd=NF_INQ_VARID(ncidd2,'dust',idd2)
          statusd=NF_INQ_VARID(ncidd3,'dust',idd3)
          statusd=NF_INQ_VARID(ncidd4,'dust',idd4)

          startd(1)=1
          startd(2)=1

          countd(1)=Im
          countd(2)=Jm
          countd(3)=1

          d_dust=0.D0

          DO k=1,JDperY

            IF (k > 59) THEN
              startd(3)=k+1
            ELSE
              startd(3)=k
            END IF

            statusd=NF_GET_VARA_REAL(ncidd1,idd1,startd,countd,work)
            d_dust(:,:,1,k)=DBLE(work(:,:))
            statusd=NF_GET_VARA_REAL(ncidd2,idd2,startd,countd,work)
            d_dust(:,:,2,k)=DBLE(work(:,:))
            statusd=NF_GET_VARA_REAL(ncidd3,idd3,startd,countd,work)
            d_dust(:,:,3,k)=DBLE(work(:,:))
            statusd=NF_GET_VARA_REAL(ncidd4,idd4,startd,countd,work)
            d_dust(:,:,4,k)=DBLE(work(:,:))

          END DO

          statusd=NF_CLOSE('dust_bin1',NCNOWRIT,ncidd1)
          statusd=NF_CLOSE('dust_bin2',NCNOWRIT,ncidd2)
          statusd=NF_CLOSE('dust_bin3',NCNOWRIT,ncidd3)
          statusd=NF_CLOSE('dust_bin4',NCNOWRIT,ncidd4)

        ELSE IF (imDUST == 0) THEN
c     interactive dust emissions
c**** Read input: threshold speed
          CALL openunit('VTRSH',io_data,.true.,.true.)
          READ (io_data) vtrsh
          CALL closeunit(io_data)
c**** Read input: fraction clay
          CALL openunit('FRCLAY',io_data,.true.,.true.)
          READ (io_data) frclay
          CALL closeunit(io_data)
c**** Read input: fraction silt
          CALL openunit('FRSILT',io_data,.true.,.true.)
          READ (io_data) frsilt
          CALL closeunit(io_data)
c**** Read input: prec-evap data
          CALL openunit('DRYHR',io_data,.true.,.true.)
          READ (io_data) dryhr
          CALL closeunit(io_data)
c**** Read input: ERS data
          call openunit('ERS',io_data,.true.,.true.)
          DO k=1,JMperY
            READ(io_data) ((ers_data(i,j,k),i=1,im),j=1,jm)
          END DO
          call closeunit(io_data)
c**** Read input: GINOUX data
          call openunit('GIN',io_data,.true.,.true.)
          read (io_data) gin_data
          call closeunit(io_data)
c**** Read input: EMISSION LOOKUP TABLE data
          call openunit('LKTAB',io_data,.true.,.true.)
          DO k=1,lkm
            READ(io_data) ((table(i,j,k),i=1,lim),j=1,ljm)
          END DO
          call closeunit(io_data)
c index of table for threshold velocity from 6.5 to 17 m/s
          DO k=1,lkm
            x3(k)=6.d0+0.5d0*k
          END DO
c index of table for sub grid scale velocity (sigma) from .0001 to 30 m/s
          sum=0.d0
          DO j=1,ljm
            IF (j <= 30) THEN
              sum=sum+0.0001d0+float(j-1)*0.00008d0
              x2(j)=sum
            ELSE IF (j > 30) THEN
              sum=sum-0.055254d0+0.005471d0*float(j)-
     &             1.938365d-4*float(j)**2.d0+
     &             3.109634d-6*float(j)**3.d0-
     &             2.126684d-8*float(j)**4.d0+
     &             5.128648d-11*float(j)**5.d0
              x2(j)=sum
            END IF
          END DO
c index of table for GCM surface wind speed from 0.0001 to 30 m/s
          x1(:)=x2(:)
        ELSE
          CALL stop_model
     &     ('Stopped in tracer_IC: parameter imDUST must be 0 or 1',255)
        END IF
#if (defined TRACERS_MINERALS) || (defined TRACERS_QUARZHEM)
        CALL openunit('MINFR',io_data,.true.,.true.)
        READ(io_data) (((minfr(i,j,k),i=1,Im),j=1,Jm),k=1,Mtrac)
        CALL closeunit(io_data)
#endif
        ifirst=.FALSE.
      END IF

#endif

      end subroutine tracer_IC


      subroutine daily_tracer(iact)
!@sum daily_tracer is called once a day for tracers
!@+   SUBROUTINE tracer_IC is called from daily_tracer to allow for
!@+     tracers that 'turn on' on different dates.
!@auth Jean Lerner
C**** Note this routine must always exist (but can be a dummy routine)
      USE MODEL_COM, only: jmon,itime,coupled_chem
      USE DOMAIN_DECOMP, only : grid, get
      USE TRACER_COM, only: ntm,trname,itime_tr0,nLightning,nAircraft
#ifdef TRACERS_SPECIAL_Lerner
      USE TRACERS_MPchem_COM, only: n_MPtable,tcscale,STRATCHEM_SETUP
      USE LINOZ_CHEM_COM, only: LINOZ_SETUP
#endif
#ifdef TRACERS_SPECIAL_Shindell
      USE FLUXES, only: tr3Dsource
      USE TRCHEM_Shindell_COM,only: PI_run
#endif
      IMPLICIT NONE
      INTEGER n,iact,last_month
      data last_month/-1/
      INTEGER J_0, J_1
C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT=J_0, J_STOP=J_1)

#ifdef TRACERS_SPECIAL_Lerner
      if (iact.eq.0) then
C**** Initialize tables for linoz
      do n=1,ntm
        if (trname(n).eq."O3" .and. itime.ge.itime_tr0(n)) then
          call linoz_setup(n)
          exit
        end if
      end do

C**** Initialize tables for Prather StratChem tracers
      do n=1,ntm
        if (trname(n).eq."N2O" .or. trname(n).eq."CH4" .or.
     *      trname(n).eq."CFC11") then
          if (itime.ge.itime_tr0(n))
     *    call stratchem_setup(n_MPtable(n),trname(n))
        end if
      end do
      end if  ! iact=0

C**** Prather StratChem tracers and linoz tables change each month
      IF (JMON.NE.last_month) THEN
        do n=1,ntm
          if ((trname(n).eq."N2O" .or. trname(n).eq."CH4" .or.
     *         trname(n).eq."CFC11") .and. itime.ge.itime_tr0(n)) then
            CALL STRTL  ! one call does all based on n_MPtable_max
            exit
          end if
        end do
        do n=1,ntm
          if (trname(n).eq."O3" .and. itime.ge.itime_tr0(n)) then
            CALL linoz_STRATL
            exit
          end if
        end do
        last_month = JMON
      END IF

C**** Tracer specific call for CO2
      do n=1,ntm
        if (trname(n).eq."CO2") then
          call read_CO2_sources(n,iact)
          exit
        end if
      end do

C**** Tracer specific call for CH4
      do n=1,ntm
        if (trname(n).eq."CH4") then
          call read_CH4_sources(n,iact)
          exit
        end if
      end do
#endif

#ifdef TRACERS_SPECIAL_Shindell
C**** Daily tracer-specific calls to read 2D and 3D sources:
          if (COUPLED_CHEM.ne.1) then
          call get_dms_offline   !not applied directly;used in chemistry.
          call get_so2_offline   !not applied directly;used in chemistry.
          endif
      do n=1,ntm
        select case (trname(n))
        case ('NOx')
          tr3Dsource(:,J_0:J_1,:,nAircraft,n)  = 0.
          if(PI_run.ne.1)call get_aircraft_NOx
          call      read_NOx_sources(n,iact)
C         (lightning called from tracer_3Dsource)
        case ('CO')
          call       read_CO_sources(n,iact)
        case ('CH4')
          call      read_CH4_sources(n,iact)
        case ('Isoprene')
          call read_Isoprene_sources(n,iact)
        case ('Alkenes')
          call  read_Alkenes_sources(n,iact)
        case ('Paraffin')
          call read_Paraffin_sources(n,iact)
        case ('N2O5')
          tr3Dsource(:,J_0:J_1,:,:,n) = 0.
          if (COUPLED_CHEM.ne.1)
     &    call get_sulfate_N2O5 !not applied directly;used in chemistry.
        end select
      end do
#endif
#ifdef TRACERS_AEROSOLS_Koch
c      if (COUPLED_CHEM.ne.1) call get_O3_offline
#endif
C****
C**** Initialize tracers here to allow for tracers that 'turn on'
C**** at the start of any day
      call tracer_IC

      return
      end subroutine daily_tracer

#ifdef TRACERS_ON
      SUBROUTINE set_tracer_2Dsource
!@sum tracer_source calculates non-interactive sources for tracers
!@auth Jean Lerner/Gavin Schmidt
      USE MODEL_COM, only: FEARTH,itime,JDperY,fland,psf,pmtop,jmpery
     *  ,dtsrc
#ifdef TRACERS_AEROSOLS_Koch
     *  ,jmon
#endif
      USE DOMAIN_DECOMP, only : GRID, GET, GLOBALSUM

      USE GEOM, only: dxyp,areag,lat_dg
      USE QUSDEF
      USE DYNAMICS, only: am  ! Air mass of each box (kg/m^2)
      USE TRACER_COM
      USE FLUXES, only: trsource
      USE SEAICE_COM, only: rsi
      USE CONSTANT, only: tf,sday,hrday,bygrav,mair
      USE PBLCOM, only: tsavg
#ifdef TRACERS_SPECIAL_Lerner
      USE CO2_SOURCES, only: co2_src
      USE CH4_SOURCES, only: ch4_src
#endif
#ifdef TRACERS_SPECIAL_Shindell
      USE TRACER_SOURCES, only: CO_src,ch4_src,NOx_src,Isoprene_src,
     &                          Alkenes_src,Paraffin_src
#endif
#ifdef TRACERS_AEROSOLS_Koch
       USE AEROSOL_SOURCES, only: SO2_src,BCI_src,OCI_src,BCB_src,
     * OCB_src,OCT_src
#endif
      implicit none
      integer :: i,j,ns,l,ky,n
      REAL*8 :: source,sarea,steppy,base,steppd,x,airm,anngas,
     *  steph,stepx,stepp,tmon,bydt,tnew
      REAL*8 :: sarea_prt(GRID%J_STRT_HALO:GRID%J_STOP_HALO)

      INTEGER J_0, J_1

C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT=J_0, J_STOP=J_1)

      bydt = 1./DTsrc
C**** All sources are saved as kg/s
      do n=1,ntm
      if (itime.lt.itime_tr0(n)) cycle
      select case (trname(n))

      case default
!     write(6,*) ' Sources for ',trname(n),' are not in this routine'
C****
C**** Surface Sources of SF6 (Same grid as CFC11)
C****
      case ('SF6','CFC11')
        trsource(:,:,:,n)=0
C**** SF6 source increases each year by .3pptv/year
C**** CFC source is the same each year
C**** Distribute source over ice-free land
        steppy = 1./(sday*JDperY)
        if (trname(n).eq.'SF6') then
C         Make sure index KY=1 in year that tracer turns on
          ky = 1 + (itime-itime_tr0(n))/(hrday*JDperY)
          base = (0.3d-12)*tr_mm(n)/mair !pptm
          x = base*ky*steppy
          airm = (psf-pmtop)*100.*bygrav*AREAG !(kg/m**2 X m**2 = kg)
          anngas = x*airm/steppy
        else
          anngas = 310.d6
        endif
C GISS-ESMF EXCEPTIONAL CASE
c   Several regions defined need to be addressed with domain decomposition
c     We can use lat_dg to distinguish regions in resolution independent
c     manner (also lon_dg for the i index), but we get a 'gather' for
c     the area sums (sarea).

C**** Source over United States and Canada
        source = .37d0*anngas*steppy
        sarea_prt(:)  = 0.
        do j=MAX(31,J_0),MIN(35,J_1)
          do i=12,22
            sarea_prt(j) = sarea_prt(j) + dxyp(j)*fearth(i,j)
          enddo
        enddo
        CALL GLOBALSUM(grid, sarea_prt, sarea, all=.true.)
        do j=J_0,J_1     ! 31,35
          if (nint(lat_dg(j,1)).ge.28 .and. nint(lat_dg(j,1)).le.48)
     *         then
            do i=12,22
              trsource(i,j,1,n) = source*dxyp(j)*fearth(i,j)/sarea
            enddo
          end if
        enddo
C**** Source over Europe and Russia
        source = .37d0*anngas*steppy
        sarea_prt(:)  = 0.
        do j=MAX(33,J_0),MIN(39,J_1)
          do i=35,45
            sarea_prt(j) = sarea_prt(j) + dxyp(j)*fearth(i,j)
          enddo
        enddo
        CALL GLOBALSUM(grid, sarea_prt, sarea, all=.true.)
        do j=J_0,J_1     ! 33,39
          if (nint(lat_dg(j,1)).ge.36 .and. nint(lat_dg(j,1)).le.64)
     *         then
          do i=35,45
            trsource(i,j,1,n) = source*dxyp(j)*fearth(i,j)/sarea
          enddo
          end if
        enddo
C**** Source over Far East
        source = .13d0*anngas*steppy
        sarea_prt  = 0.
        do j=MAX(29,J_0),MIN(34,J_1)
          do i=61,66
            sarea_prt(j) = sarea_prt(j) + dxyp(j)*fearth(i,j)
          enddo
        enddo
        CALL GLOBALSUM(grid, sarea_prt, sarea, all=.true.)
        do j=J_0,J_1     ! 29,34
          if (nint(lat_dg(j,1)).ge.20 .and. nint(lat_dg(j,1)).le.44)
     *         then
          do i=61,66
            trsource(i,j,1,n) = source*dxyp(j)*fearth(i,j)/sarea
          enddo
          end if
        enddo
C**** Source over Middle East
        source = .05d0*anngas*steppy
        sarea_prt  = 0.
        do j=MAX(28,J_0),MIN(32,J_1)
          do i=43,51
            sarea_prt(j) = sarea_prt(j) + dxyp(j)*fearth(i,j)
          enddo
        enddo
        CALL GLOBALSUM(grid, sarea_prt, sarea, all=.true.)
        do j=J_0,J_1     ! 28,32
          if (nint(lat_dg(j,1)).ge.16 .and. nint(lat_dg(j,1)).le.36)
     *         then
          do i=43,51
            trsource(i,j,1,n) = source*dxyp(j)*fearth(i,j)/sarea
          enddo
          end if
        enddo
C**** Source over South America
        source = .04d0*anngas*steppy
c        do j=J_0,J_1   ! example coding
c        if (nint(lat_dg(j,1)).ge.-24 .and. nint(lat_dg(j,1)).le.-20)
c     *         then
        i=27; j=18
        IF (j >= J_0 .and. j <= J_1) trsource(i,j,1,n) = 0.5*source
        i=28; j=18
        IF (j >= J_0 .and. j <= J_1) trsource(i,j,1,n) = 0.5*source
c        end if
c        end do
C**** Source over South Africa
        source = .02d0*anngas*steppy
        i=42; j=17
        IF (j >= J_0 .and. j <= J_1) trsource(i,j,1,n) = source
C**** Source over Australia and New Zealand
        source = .02d0*anngas*steppy
        i=66; j=15
        IF (j >= J_0 .and. j <= J_1) trsource(i,j,1,n) = source

C****
C**** Surface Sources for Radon-222
C****
      case ('Rn222')
        trsource(:,J_0:J_1,:,n)=0
C**** ground source
        steppd = 1./sday
        do j=J_0,J_1
          do i=1,im
C**** source from ice-free land
            if(tsavg(i,j).lt.tf) then !composite surface air temperature
              trsource(i,j,1,n) = 1.0d-16*steppd*dxyp(j)*fearth(i,j)
            else  ! 1 atom/cm^2/s
              trsource(i,j,1,n) = 3.2d-16*steppd*dxyp(j)*fearth(i,j)
            end if
C**** source from ice-free ocean
            trsource(i,j,1,n) =trsource(i,j,1,n)+ 1.6d-18*steppd*dxyp(j)
     *           *(1.-fland(i,j))*(1.-rsi(i,j))
          enddo                 !i
        enddo                   !j

#ifdef TRACERS_SPECIAL_Lerner
C****
C**** Sources and sinks for CO2 (kg/s)
C****
      case ('CO2')
        do ns=1,ntsurfsrc(n)
          do j=J_0,J_1
            trsource(:,j,ns,n) = co2_src(:,j,ns)*dxyp(j)
          end do
        end do

C****
C**** Sources and sinks for CH4 (kg/s)
C****
      case ('CH4')
        do ns=1,ntsurfsrc(n)
          do j=J_0,J_1
            trsource(:,j,ns,n) = ch4_src(:,j,ns)*dxyp(j)
          end do
        end do
C****
C**** Sources and sinks for N2O:
C**** First layer is set to a constant 462.2 ppbm. (300 PPB V)
C****
      case ('N2O')
      do j=J_0,J_1
        trsource(:,j,1,n) = (am(1,:,j)*dxyp(j)*462.2d-9
     *   -trm(:,j,1,n))*bydt
      end do
#endif
C****
C**** Sources and sinks for 14C02
C**** NOTE: This tracer is supposed to start on 10/16
C**** The tracer is reset to specific values in layer 1 only if
C****   this results in a sink
C****
      case ('14CO2')
      tmon = (itime-itime_tr0(n))*jmpery/(hrday*jdpery)  !(12./8760.)
      trsource(:,J_0:J_1,1,n) = 0.
      do j=J_0,J_1
         if (j <= jm/2) then
            do i=1,im
               tnew = am(1,i,j)*dxyp(j)*(4.82d-18*46./mair)*
     *          (44.5 + tmon*(1.02535d0 - tmon*
     *                  (2.13565d-2 - tmon*8.61853d-5)))
               if (tnew.lt.trm(i,j,1,n))
     *             trsource(i,j,1,n) = (tnew-trm(i,j,1,n))*bydt
            end do
         else
            do i=1,im
               tnew = am(1,i,j)*dxyp(j)*(4.82d-18*46./mair)*
     *          (73.0 - tmon*(0.27823d0 + tmon*
     *                  (3.45648d-3 - tmon*4.21159d-5)))
               if (tnew.lt.trm(i,j,1,n))
     *             trsource(i,j,1,n) = (tnew-trm(i,j,1,n))*bydt
            end do
         endif
      end do

C****
C**** No non-interactive surface sources of Water
C****
      case ('Water')
        trsource(:,J_0:J_1,:,n)=0

#ifdef TRACERS_SPECIAL_Shindell
      case ('CO')
        do ns=1,ntsurfsrc(n)
         do j=J_0,J_1
            trsource(:,j,ns,n) = CO_src(:,j,ns)*dxyp(j)
         end do
        end do

      case ('CH4')
        do ns=1,ntsurfsrc(n)
          do j=J_0,J_1
            trsource(:,j,ns,n) = ch4_src(:,j,ns)*dxyp(j)
          end do
        end do

      case ('Alkenes')
        do ns=1,ntsurfsrc(n)
         do j=J_0,J_1
            trsource(:,j,ns,n) = Alkenes_src(:,j,ns)*dxyp(j)
         end do
        end do

      case ('Paraffin')
        do ns=1,ntsurfsrc(n)
          do j=J_0,J_1
            trsource(:,j,ns,n) = Paraffin_src(:,j,ns)*dxyp(j)
          end do
        end do

      case ('Isoprene')
        do ns=1,ntsurfsrc(n)
         do j=J_0,J_1
            trsource(:,j,ns,n) = Isoprene_src(:,j,ns)*dxyp(j)
         end do
        end do

      case ('NOx')
        do ns=1,ntsurfsrc(n)
          do j=J_0,J_1
            trsource(:,j,ns,n) = NOx_src(:,j,ns)*dxyp(j)
          end do
        end do
#endif
#ifdef TRACERS_AEROSOLS_Koch
      case ('SO2')
        if (imAER.eq.0) then
        call read_SO2_source(n)
        endif
        do ns=1,ntsurfsrc(n)
         do j=J_0,J_1
            if (ntsurfsrc(n).eq.1) then
            trsource(:,j,ns,n) = SO2_src(:,j,ns)*0.975d0
            else
            trsource(:,j,ns,n) = SO2_src(:,j,ns)
            endif
         end do
        end do
c we assume 97.5% emission as SO2, 2.5% as sulfate (*tr_mm/tr_mm)
      case ('SO4')
         do j=J_0,J_1
            trsource(:,j,1,n) = SO2_src(:,j,1)*0.0375d0
         end do
      case ('BCII')
         do j=J_0,J_1
            trsource(:,j,1,n) = BCI_src(:,j)  !*3.d0
c arbitrary factor of 3
         end do
      case ('OCII')
         do j=J_0,J_1
            trsource(:,j,1,n) = OCI_src(:,j)  !*2.d0
c arbitrary factor of 2
            trsource(:,j,2,n) = OCT_src(:,j,jmon)
         end do
#endif
      end select

      end do
C****
      END SUBROUTINE set_tracer_2Dsource


      SUBROUTINE tracer_3Dsource
!@sum tracer_3Dsource calculates interactive sources for tracers
!@+   Please note that if the generic routine 'apply_tracer_3Dsource'
!@+   is used, all diagnostics and moments are updated automatically.
!@auth Jean Lerner/Greg Faluvegi
!@calls DIAGTCA, masterchem, apply_tracer_3Dsource
      USE DOMAIN_DECOMP, only : GRID, GET
      USE TRACER_COM
      USE FLUXES, only: tr3Dsource
      USE MODEL_COM, only: itime,jmon
      USE GEOM, only: dxyp,bydxyp
      USE DYNAMICS, only: am,byam ! Air mass of each box (kg/m^2)
      USE apply3d, only : apply_tracer_3Dsource
#ifdef TRACERS_COSMO
      USE COSMO_SOURCES, only: be7_src_3d, be7_src_param
#endif
#ifdef TRACERS_AEROSOLS_Koch
      USE AEROSOL_SOURCES, only: SO2_src_3d,BCI_src_3d,BCB_src,
     *     OCB_src,SO2_biosrc_3D,lmAER
      USE PBLCOM, only: dclev
c Laki emissions
      USE DYNAMICS, only: LTROPO
      USE CONSTANT, only: sday
      USE MODEL_COM, only: jday,jyear
      USE LAKI_SOURCE, only: LAKI_MON,LAKI_DAY,LAKI_AMT_T,LAKI_AMT_S
#endif
      implicit none
      INTEGER n,ns,najl,i,j,l,mnow,blay
      INTEGER J_0, J_1
      INTEGER :: ltpp
      
C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT=J_0, J_STOP=J_1)

C**** All sources are saved as kg/s
      do n=1,ntm

      if (itime.lt.itime_tr0(n)) cycle
      select case (trname(n))

      case default
#ifdef TRACERS_SPECIAL_Lerner
C****
      case ('CH4')
      tr3Dsource(:,J_0:J_1,:,:,n) = 0.
      call Trop_chem_CH4(1,n)
      call apply_tracer_3Dsource(1,n)
      call Strat_chem_Prather(2,n)
      call apply_tracer_3Dsource(2,n,.false.)
cc      call DIAGTCA(itcon_3Dsrc(2,n),n)
C****
      case ('O3')
      tr3Dsource(:,J_0:J_1,:,:,n) = 0.
      call Trop_chem_O3(1,n)
      call apply_tracer_3Dsource(1,n)
cc      call DIAGTCA(itcon_3Dsrc(1,n),n)
      call Strat_chem_O3(2,n)
      call apply_tracer_3Dsource(2,n)
cc      call DIAGTCA(itcon_3Dsrc(2,n),n)
C****
      case ('N2O')
      tr3Dsource(:,J_0:J_1,:,:,n) = 0.
      call Strat_chem_Prather(1,n)
      call apply_tracer_3Dsource(1,n,.FALSE.)
cc      call DIAGTCA(itcon_3Dsrc(1,n),n)
C****
      case ('CFC11')
      tr3Dsource(:,J_0:J_1,:,:,n) = 0.
      call Strat_chem_Prather(1,n)
      call apply_tracer_3Dsource(1,n,.FALSE.)
cc      call DIAGTCA(itcon_3Dsrc(1,n),n)
C****
#endif
#ifdef TRACERS_AEROSOLS_Koch
      case ('SO2')
C**** three 3D sources (aircraft, volcanos and biomass) read in from files
c  Laki
      if (JYEAR.eq.1783) then
       do j=1,10
       if (JMON.eq.LAKI_MON(j).and.JDAY.eq.LAKI_DAY(j)) then
       do l=1,5
        SO2_src_3d(33,40,l,1)=SO2_src_3d(33,40,l,1)+LAKI_AMT_T(j)
     *          /sday*1000.d0/5.d0
       end do
       ltpp=LTROPO(33,40)
       do l=ltpp+1,ltpp+3
        SO2_src_3d(33,40,l,1)=SO2_src_3d(33,40,l,1)+LAKI_AMT_S(j)
     *          /sday*1000.d0/3.d0
       end do
       endif
       enddo
      endif
c End Laki code
      tr3Dsource(:,J_0:J_1,:,1,n) = SO2_src_3d(:,J_0:J_1,:,1)*0.975d0
      call apply_tracer_3Dsource(1,n) ! volcanos
      tr3Dsource(:,J_0:J_1,:,2,n) = SO2_src_3d(:,J_0:J_1,:,2)
      if (imAER.eq.0) call apply_tracer_3Dsource(2,n) ! aircraft
      tr3Dsource(:,J_0:J_1,1:lmAER,3,n)=SO2_biosrc_3D(:,J_0:J_1,:,jmon)
     *     *0.975d0
      tr3Dsource(:,J_0:J_1,lmAER+1:lm,3,n) = 0.
      call apply_tracer_3Dsource(3,n) ! biomass

      case ('SO4')
C**** three 3D sources ( volcanos and biomass) read in from files
      tr3Dsource(:,J_0:J_1,:,3,n) = SO2_src_3d(:,J_0:J_1,:,1)*0.0375d0
      call apply_tracer_3Dsource(3,n) ! volcanos
      tr3Dsource(:,J_0:J_1,1:lmAER,4,n)=SO2_biosrc_3D(:,J_0:J_1,:,jmon)
     *     *0.0375d0
      tr3Dsource(:,J_0:J_1,lmAER+1:lm,4,n) = 0.
      call apply_tracer_3Dsource(4,n) ! biomass

       case ('BCIA')
C**** aircraft source for fresh industrial BC 
      tr3Dsource(:,J_0:J_1,:,2,n) = BCI_src_3d(:,J_0:J_1,:)
      call apply_tracer_3Dsource(2,n) ! aircraft

       case ('BCB')
C**** biomass source for BC 
      tr3Dsource(:,J_0:J_1,:,1,n) = 0.
      if (imAER.eq.0) then
      do j=J_0,J_1; do i=1,im
      blay=int(dclev(i,j)+0.5)
      do l=1,blay
      tr3Dsource(i,j,l,1,n) = BCB_src(i,j,1,jmon)/real(blay)
      end do
      end do; end do
      else
      tr3Dsource(:,J_0:J_1,1:lmAER,1,n) = BCB_src(:,J_0:J_1,:,jmon)
      endif
      call apply_tracer_3Dsource(1,n) ! biomass 
       case ('OCB')
C**** biomass source for OC 
      tr3Dsource(:,J_0:J_1,:,1,n) = 0.
      if (imAER.eq.0) then
      do j=J_0,J_1; do i=1,im
      blay=int(dclev(i,j)+0.5)
      do l=1,blay
      tr3Dsource(i,j,l,1,n) = OCB_src(i,j,1,jmon)/real(blay)
      end do
      end do; end do
      else
      tr3Dsource(:,J_0:J_1,1:lmAER,1,n) = OCB_src(:,J_0:J_1,:,jmon)
      endif
      call apply_tracer_3Dsource(1,n) ! biomass 
#endif
#ifdef TRACERS_COSMO
C****
      case ('Be7')
c cosmogenic src
        do l=1,lm; do j=J_0,J_1; do i=1,im
          tr3Dsource(i,j,l,1,n) = be7_src_param * am(l,i,j) *
     *         be7_src_3d(i,j,l)
        end do; end do; end do
        call apply_tracer_3Dsource(1,n)
C****
      case ('Be10')
c 0.52 is ratio of Be10 to Be7 production
c tr_mm(n_Be10)/tr_mm(n_Be7)= 10./7. is ratio of molecular weights
c cosmogenic src
        do l=1,lm; do j=J_0,J_1; do i=1,im
          tr3Dsource(i,j,l,1,n)=0.52d0 * be7_src_param * am(l,i,j)
     *         * be7_src_3d(i,j,l) * tr_mm(n_Be10)/tr_mm(n_Be7)
        end do; end do; end do
        call apply_tracer_3Dsource(1,n)
C****
      case('Pb210')
        call apply_tracer_3Dsource(1,n) !radioactive decay of Rn222
#endif

      end select

      end do

#ifdef TRACERS_AEROSOLS_Koch
       call aerosol_gas_chem
       call apply_tracer_3Dsource(1,n_DMS)  ! DMS chem sink
       call apply_tracer_3Dsource(1,n_MSA)  ! MSA chem source
       call apply_tracer_3Dsource(4,n_SO2)  ! SO2 chem source
       call apply_tracer_3Dsource(5,n_SO2)  ! SO2 chem sink
       call apply_tracer_3Dsource(1,n_SO4)  ! SO4 chem source
       call apply_tracer_3Dsource(1,n_H2O2_s) ! H2O2 chem source
       call apply_tracer_3Dsource(2,n_H2O2_s) ! H2O2 chem sink
       call apply_tracer_3Dsource(1,n_BCII) ! BCII aging sink 
       call apply_tracer_3Dsource(1,n_BCIA) ! BCIA aging source 
       call apply_tracer_3Dsource(1,n_OCII) ! OCII aging sink 
       call apply_tracer_3Dsource(1,n_OCIA) ! OCIA aging source 

       call heter
       call apply_tracer_3Dsource(6,n_SO2) ! SO2 het chem sink
       call apply_tracer_3Dsource(2,n_SO4) ! SO4 het chem source
       call apply_tracer_3Dsource(3,n_H2O2_s) ! H2O2 het chem sink

#ifdef TRACERS_HETCHEM
       call apply_tracer_3Dsource(1,n_SO4_s1) ! SO4 chem prod on ss1
       call apply_tracer_3Dsource(1,n_SO4_s2) ! SO4 chem prod on ss2
       call apply_tracer_3Dsource(1,n_SO4_d1) ! SO4 chem prod on dust
       call apply_tracer_3Dsource(1,n_SO4_d2) ! SO4 chem prod on dust
       call apply_tracer_3Dsource(1,n_SO4_d3) ! SO4 chem prod on dust
       call apply_tracer_3Dsource(1,n_SO4_d4) ! SO4 chem prod on dust
#endif
#endif

#ifdef TRACERS_SPECIAL_Shindell
C Apply non-chemistry 3D sources, so they can be "seen" by chemistry:
C (Note: using this method, tracer moments are changed just like they
C are done for chemistry.  It might be better to do it like surface
C sources are done? -- GSF 11/26/02)
c
      CALL TIMER (MNOW,MTRACE)
      call apply_tracer_3Dsource(nAircraft,n_NOx)
      tr3Dsource(:,J_0:J_1,:,nLightning,n_NOx) = 0.
      call get_lightning_NOx
      call apply_tracer_3Dsource(nLightning,n_NOx)
c
C**** Make sure that these 3D sources for all chem tracers start at 0.:
      tr3Dsource(:,J_0:J_1,:,nChemistry,1:ntm_chem)  = 0.
#ifndef SHINDELL_STRAT_CHEM
      tr3Dsource(:,J_0:J_1,:,nStratwrite,1:ntm_chem) = 0.
#endif
C**** Call the model CHEMISTRY and STRATOSPHERE OVERWRITE:

      CALL masterchem ! does chemistry and stratospheric over-writing.
                      ! tr3Dsource defined within, for both processes

C**** Apply chemistry and stratosphere overwrite changes:

      do n=1,ntm_chem
        call apply_tracer_3Dsource(nChemistry,n)
#ifndef SHINDELL_STRAT_CHEM
        call apply_tracer_3Dsource(nStratwrite,n)
#endif
      end do
      CALL TIMER (MNOW,MCHEM)
#endif

      return
      END SUBROUTINE tracer_3Dsource
#endif

#ifdef TRACERS_WATER
C---SUBROUTINES FOR TRACER WET DEPOSITION-------------------------------

      SUBROUTINE GET_COND_FACTOR(L,N,WMXTR,TEMP,TEMP0,LHX,FCLOUD,FQ0,fq,
     *  TR_CONV,TRWML,TM,THLAW,TR_LEF,pl,ntix,CLDSAVT)
!@sum  GET_COND_FACTOR calculation of condensate fraction for tracers
!@+    within or below convective or large-scale clouds. Gas
!@+    condensation uses Henry's Law if not freezing.
!@auth Dorothy Koch (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on CB436TdsM23 CLOUDCHCC and CLOUDCHEM subroutines)
c
C**** GLOBAL parameters and variables:
      USE CONSTANT, only: BYGASC, MAIR,teeny,lhe,tf,by3
      USE TRACER_COM, only: tr_RKD,tr_DHD,nWATER,nGAS,nPART,tr_wd_TYPE
     *     ,trname,ntm,lm,t_qlimit,fq_aer
#ifdef TRACERS_SPECIAL_O18
     &     ,supsatfac
#endif
#ifdef TRACERS_HETCHEM
     *     ,trm ,n_SO4_d1, n_SO4_d2, n_SO4_d3, n_SO4_d4

      USE MODEL_COM, only  : dtsrc
#endif
c      USE CLOUDS, only: PL, NTIX
c
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
c
!@param BY298K unknown meaning for now (assumed= 1./298K)
!@var Ppas pressure at current altitude (in Pascal=kg/s2/m)
!@var TFAC exponential coeffiecient of tracer condensation temperature
!@+   dependence (mole/joule)
!@var FCLOUD fraction of cloud available for tracer condensation
!@var SSFAC dummy variable (assumed units= kg water?)
!@var FQ            fraction of tracer that goes into condensate
!@var FQ0 default fraction of water tracer that goes into condensate
!@var L index for altitude loop
!@var N index for tracer number loop
!@var WMXTR mixing ratio of water available for tracer condensation?
!@var SUPSAT super-saturation ratio for cloud droplets
!@var LHX latent heat flag for whether condensation is to ice or water
!@var RKD dummy variable (= tr_RKD*EXP[ ])
      REAL*8, PARAMETER :: BY298K=3.3557D-3
      REAL*8 Ppas, tfac, ssfac, RKD,CLDINC
#ifdef TRACERS_SPECIAL_O18
      real*8 tdegc,alph,fracvs,fracvl,kin_cond_ice,fqi,gint
      integer i
!@param nstep no. of steps for integration of Rayleigh condensation
      integer, parameter :: nstep=6   !8
!@param wgt weightings for simpson's rule integration
      real*8, parameter, dimension(nstep+1) ::
     *     wgt = (/ by3, 4*by3, 2*by3, 4*by3, 2*by3, 4*by3, by3 /)
c     *     wgt = (/ by3, 4*by3, 2*by3, 4*by3, 2*by3, 4*by3, 2*by3, 4*by3
c     *     , by3 /)
#endif
      REAL*8,  INTENT(IN) :: fq0, FCLOUD, WMXTR, TEMP, TEMP0,LHX, TR_LEF
     *     , pl,CLDSAVT
      REAL*8,  INTENT(IN), DIMENSION(ntm,lm) :: trwml
      REAL*8,  INTENT(IN), DIMENSION(lm,ntm) :: TM
      REAL*8,  INTENT(OUT):: fq,thlaw
      INTEGER, INTENT(IN) :: L, N, ntix(ntm)
      LOGICAL TR_CONV
      REAL*8 :: SUPSAT
c
C**** CALCULATE the fraction of tracer mass that becomes condensate:
c
      thlaw=0.
      SELECT CASE(tr_wd_TYPE(NTIX(N)))
        CASE(nGAS)                            ! gas tracer
          fq = 0.D0                           ! frozen and default case
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_SPECIAL_Shindell)
          IF(LHX.eq.LHE) THEN                 ! if not frozen then:
            Ppas = PL*1.D2                 ! pressure to pascals
            tfac = (1.D0/TEMP - BY298K)*BYGASC
            IF(tr_DHD(NTIX(N)).ne.0.D0) THEN
              RKD=tr_RKD(NTIX(N))*DEXP(-tr_DHD(NTIX(N))*tfac)
            ELSE
              RKD=tr_RKD(NTIX(N))
            END IF
c           clwc=WMXTR*MAIR*1.D-3*Ppas*BYGASC/(TEMP*FCLOUD)
c           ssfac=RKD*GASC*TEMP*clwc   ! Henry's Law
            ssfac=RKD*WMXTR*MAIR*1.D-3*Ppas/(CLDSAVT+teeny)
            if (.not.tr_conv) then  !stratiform
              thlaw=(ssfac*tr_lef*tm(l,NTIX(N))
     *        -TRWML(NTIX(N),L))/(1.D0+ssfac)
              if (thlaw.lt.0.) thlaw=0.d0
              if (thlaw.gt.tm(l,NTIX(N))) thlaw=tm(l,NTIX(N))
            else  !if convection
              fq=ssfac / (1.D0 + ssfac)
              if (fq.ge.1.) fq=1.d0
              thlaw=0.
            endif
            if (FCLOUD.LT.1.D-16) fq=0.d0
            if (FCLOUD.LT.1.D-16) thlaw=0.d0
#ifdef TRACERS_SPECIAL_Shindell
            if(t_qlimit(NTIX(N)).and.fq.gt.1.)fq=1.!no negative tracers
#endif
          END IF
#endif
        CASE(nWATER)                          ! water tracer
#ifdef TRACERS_SPECIAL_O18
          if (fq0.gt.0. .and. fq0.lt.1.) then
C**** If process occurs at constant temperature, calculate condensate
C**** in equilibrium with source vapour. Otherwise, use mid-point
C**** temperature and estimate instantaneous fractionation. This gives
C**** a very good estimate to complete integral
C****
            if (abs(temp-temp0).gt.1d-14) then  ! use instantaneous frac
              tdegc=0.5*(temp0 + temp) -tf
C**** Calculate alpha (fractionation coefficient)
                if (LHX.eq.LHE) then ! cond to water
                  alph=1./fracvl(tdegc,trname(ntix(n)))
                else            ! cond to ice
                  alph=1./fracvs(tdegc,trname(ntix(n)))
C**** kinetic fractionation can occur as a function of supersaturation
C**** this is a parameterisation from Georg Hoffmann
                  supsat=1d0-supsatfac*tdegc
                  if (supsat .gt. 1.) alph=kin_cond_ice(alph,supsat
     *                 ,trname(ntix(n)))
                end if
                fq = 1.- (1.-fq0)**alph
            else
C**** assume condensate in equilibrium with vapour at temp
              tdegc=temp -tf
              if (LHX.eq.LHE) then ! cond to water
                alph=1./fracvl(tdegc,trname(ntix(n)))
              else              ! cond to ice
                alph=1./fracvs(tdegc,trname(ntix(n)))
C**** kinetic fractionation can occur as a function of supersaturation
C**** this is a parameterisation from Georg Hoffmann
                supsat=1d0-supsatfac*tdegc
                if (supsat .gt. 1.) alph=kin_cond_ice(alph,supsat
     *               ,trname(ntix(n)))
              end if
              fq = alph * fq0/(1.+(alph-1.)*fq0)
            end if
          else
            fq = fq0
          end if
#else
          fq = fq0
#endif
        CASE(nPART)                           ! particulate tracer
          fq = 0.D0                           ! defaults to zero.
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_COSMO) ||\
    (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
c only dissolve if the cloud has grown
c#ifdef TRACERS_HETCHEM
c      select case(trname(ntix(n)))
c      case('Clay')
c         if (trm(ii,j,l,n_SO4_d1)  > 0.1 ) then
c            fq_aer(NTIX(N))  = 1.
c         else
c            fq_aer(NTIX(N))  = 0.
c         endif
    
c      case('Silt1')
c         if (trm(ii,j,l,n_SO4_d2)  > 0.1 ) then
c            fq_aer(NTIX(N))  = 1.
c         else
c            fq_aer(NTIX(N))  = 0.
c         endif
   
c      case('Silt2')
c         if (trm(ii,j,l,n_SO4_d3)  > 0.1 ) then
c            fq_aer(NTIX(N))  = 1.
c         else
c            fq_aer(NTIX(N))  = 0.
c         endif
c      case('Silt3')
c         if (trm(ii,j,l,n_SO4_d4)  > 0.1 ) then
c            fq_aer(NTIX(N))  = 1.
c         else
c            fq_aer(NTIX(N))  = 0.
c         endif
c      end select
c#endif
           CLDINC=CLDSAVT-FCLOUD 
          if (fq0.gt.0.and.CLDINC.gt.0.) then 
          if(LHX.EQ.LHE) then !liquid cloud
            fq = fq_aer(NTIX(N))*CLDINC
           else ! ice cloud - small dissolution
            fq = fq_aer(NTIX(N))*CLDINC*0.05d0
           endif
          endif
c complete dissolution in convective clouds
c with double dissolution if partially soluble
          if (TR_CONV) then
           if (LHX.EQ.LHE) then !liquid cloud
               fq=(1.d0+fq_aer(ntix(n)))/2.d0
c              fq=(1.d0+3.d0*fq_aer(ntix(n)))/4.d0
           else
               fq=(1.d0+fq_aer(ntix(n)))/2.d0*0.05d0
c              fq=(1.d0+3.d0*fq_aer(ntix(n)))/4.d0*0.05d0
           endif
          endif
          if (FCLOUD.LT.1.D-16) fq=0.d0
          if (fq.ge.1.d0) fq=0.9999
#endif

        CASE DEFAULT                                ! error
          call stop_model(
     &    'tr_wd_TYPE(NTIX(N)) out of range in GET_COND_FACTOR',255)
      END SELECT
c
      RETURN
      END SUBROUTINE GET_COND_FACTOR


      SUBROUTINE GET_PREC_FACTOR(N,BELOW_CLOUD,CM,FCLD,FQ0,fq,ntix)
!@sum  GET_PREC_FACTOR calculation of the precipitation scavenging
!@+    fraction for tracers WITHIN large scale clouds. Current version
!@+    uses the first order removal rate based on [Giorgi and
!@+    Chameides, 1986], for gaseous and particulate tracers.
!@auth Dorothy Koch (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on CB436TdsM23 RAINOUT subroutine)
c
C**** GLOBAL parameters and variables:
      USE MODEL_COM, only: dtsrc
      USE TRACER_COM, only: nWATER, nGAS, nPART, tr_wd_TYPE,ntm
c      USE CLOUDS, only: NTIX
c
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
!@var FQ            tracer fraction scavenged into precipitation
!@var FQ0 [default] tracer fraction scavenged into precipitation
!@var FCLD cloud fraction
!@var N index for tracer number loop
!@var CM conversion rate for large cloud water content
!@var BELOW_CLOUD logical- is the current level below cloud?
      LOGICAL, INTENT(IN) :: BELOW_CLOUD
      INTEGER, INTENT(IN) :: N,ntix(ntm)
      REAL*8,  INTENT(IN) :: FQ0, FCLD, CM
      REAL*8,  INTENT(OUT):: FQ
c
      SELECT CASE(tr_wd_TYPE(NTIX(N)))
        CASE(nGAS)                                ! gas
          IF(BELOW_CLOUD) THEN
            fq = 0.D0
          ELSE
c           minus preserves FPRT sign convention in LSCOND
c           fq = -(CM*DTsrc*(EXP(-CM*DTsrc)- 1.0)*FCLD)
            fq=FQ0
          END IF
        CASE(nWATER)                          ! water/original method
          fq = FQ0                            ! no fractionation
        CASE(nPART)                           ! aerosols
          IF(BELOW_CLOUD) THEN
            fq = 0.D0
          ELSE
c           minus preserves FPRT sign convention in LSCOND
c           fq = -(CM*DTsrc*(EXP(-CM*DTsrc)- 1.0)*FCLD)
c We could try the following, to make removal amount proportional
c   to autoconversion (CM*DTsrc). The 2nd factor of CM*DTsrc
c   approximates the fraction of the cloud that is precipitating
            fq =  FQ0   ! CM*DTsrc*CM*DTsrc
          END IF
        CASE DEFAULT                          ! error
          call stop_model(
     &    'tr_wd_TYPE(NTIX(N)) out of range in GET_FPRT',255)
      END SELECT
c
      RETURN
      END SUBROUTINE GET_PREC_FACTOR


      SUBROUTINE GET_WASH_FACTOR(N,b_beta_DT,PREC,fq
     * ,TEMP,LHX,WMXTR,FCLOUD,L,TM,TRPR,THLAW,pl,ntix)
!@sum  GET_WASH_FACTOR calculation of the fraction of tracer
!@+    scavanged by precipitation below convective clouds ("washout").
!@auth Dorothy Koch (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on CB436TdsM23 CWASH and WASH_EVAP routines)
c
C**** GLOBAL parameters and variables:
      USE TRACER_COM, only: nWATER, nGAS, nPART, tr_wd_TYPE,
     * tr_RKD,tr_DHD,LM,NTM
#ifdef TRACERS_AEROSOLS_Koch
     * ,trname,n_seasalt1,n_seasalt2
c     USE PBLCOM, only: wsavg
#endif
c      USE CLOUDS, only: NTIX,PL
      USE CONSTANT, only: BYGASC,LHE,MAIR,teeny
c
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
!@var FQ fraction of tracer scavenged by below-cloud precipitation
!@param rc_wash aerosol washout rate constant (mm-1)
!@var PREC precipitation amount from layer above for washout (mm)
!@var b_beta_DT precipitating grid box fraction from lowest
!@+   percipitating layer.
!@+   The name was chosen to correspond to Koch et al. p. 23,802.
!@var N index for tracer number loop
      INTEGER, INTENT(IN) :: N,L,ntix(ntm)
      REAL*8, INTENT(OUT):: FQ,THLAW
      REAL*8, INTENT(IN) :: PREC,b_beta_DT,TEMP,LHX,WMXTR,FCLOUD,
     *  TM(LM,NTM),pl
      REAL*8, PARAMETER :: rc_wash = 1.D-1, BY298K=3.3557D-3
      REAL*8 Ppas, tfac, ssfac, RKD, TRPR(ntm)
C
      thlaw=0.
      SELECT CASE(tr_wd_TYPE(NTIX(N)))
        CASE(nGAS)                            ! gas
          fq = 0.D0                           ! frozen and default case
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_SPECIAL_Shindell)
          IF(LHX.EQ.LHE) THEN                 ! if not frozen then:
            Ppas = PL*1.D2                 ! pressure to pascals
            tfac = (1.D0/TEMP - BY298K)*BYGASC
            IF(tr_DHD(NTIX(N)).ne.0.D0) THEN
              RKD=tr_RKD(NTIX(N))*DEXP(-tr_DHD(NTIX(N))*tfac)
            ELSE
              RKD=tr_RKD(NTIX(N))
            END IF
            ssfac=RKD*WMXTR*MAIR*1.D-3*Ppas/(FCLOUD+teeny)
            thlaw=(ssfac*tm(l,NTIX(N))-TRPR(NTIX(N)))
     *            /(1.D0+ssfac)
            if (thlaw.lt.0.) thlaw=0.d0
            if (thlaw.gt.tm(l,NTIX(N))) thlaw=tm(l,NTIX(N))
            if (FCLOUD.lt.1.D-16) fq=0.d0
            if (FCLOUD.LT.1.D-16) thlaw=0.d0
          ENDIF
#endif
        CASE(nWATER)                          ! water/original method
          fq = 0.D0
        CASE(nPART)                           ! aerosols
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_COSMO) ||\
    (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
          fq = -b_beta_DT*(EXP(-PREC*rc_wash)-1.D0)
          if (FCLOUD.lt.1.D-16) fq=0.d0
          if (fq.lt.0.) fq=0.d0
c         if (wsavg(i,j).gt.10.and.PREC.gt.0.and.l.eq.1) then
c         select case (trname(n))
c          case('seasalt1')
c          fq=0.
c          case('seasalt2')
c          fq=0.
c         end select
c         endif
#endif
        CASE DEFAULT                          ! error
          call stop_model(
     &    'tr_wd_TYPE(NTIX(N)) out of range in WASHOUT_TRACER',255)
      END SELECT
c
      RETURN
      END SUBROUTINE GET_WASH_FACTOR

      SUBROUTINE GET_EVAP_FACTOR(N,TEMP,LHX,QBELOW,HEFF,FQ0,fq,ntix)
!@sum  GET_EVAP_FACTOR calculation of the evaporation fraction
!@+    for tracers.
!@auth Dorothy Koch (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on CB436TdsM23 EVAPD and WASH_EVAP routines)
c
C**** GLOBAL parameters and variables:
      USE CONSTANT, only : tf,lhe
      USE TRACER_COM, only: ntm,tr_evap_fact, tr_wd_TYPE,nwater,trname
c      USE CLOUDS, only: NTIX
c
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
!@var FQ            fraction of tracer evaporated
!@var FQ0 [default] fraction of tracer evaporated
!@var N index for tracer number loop
      INTEGER, INTENT(IN) :: N,ntix(ntm)
      REAL*8,  INTENT(OUT):: FQ
      REAL*8,  INTENT(IN) :: FQ0,TEMP,LHX
!@var QBELOW true if evap is occuring below cloud
      LOGICAL, INTENT(IN) :: QBELOW
!@var HEFF effective relative humidity for evap occuring below cloud
      REAL*8, INTENT(IN) :: HEFF
#ifdef TRACERS_SPECIAL_O18
      real*8 tdegc,alph,fracvl,fracvs,kin_evap_prec
#endif
c
      select case (tr_wd_TYPE(NTIX(N)))
      case default
        fq=FQ0*tr_evap_fact(tr_wd_TYPE(NTIX(N)))
        if(FQ0.ge.1.) fq=1.D0 ! total evaporation
c
      case (nWater)
#ifdef TRACERS_SPECIAL_O18
          tdegc=temp-tf
          if (lhx.eq.lhe) then
            alph=fracvl(tdegc,trname(ntix(n)))
C**** below clouds kinetic effects with evap into unsaturated air
            if (QBELOW.and.heff.lt.1.) alph=kin_evap_prec(alph,heff
     *           ,trname(ntix(n)))
          else
C**** no fractionation for ice evap
            alph=1.
          end if
          if (fq0.ne.1.) then
             fq = 1. - (1.-fq0)**alph
          else
            fq = fq0
          end if
#else
          fq=FQ0*tr_evap_fact(tr_wd_TYPE(NTIX(N)))
          if(FQ0.ge.1.) fq=1.D0 ! total evaporation
#endif
      end select
      RETURN
      END SUBROUTINE GET_EVAP_FACTOR
#endif

#if (defined TRACERS_SPECIAL_Shindell) || (defined TRACERS_AEROSOLS_Koch)
      SUBROUTINE GET_SULF_GAS_RATES
!@sum  GET_SULF_GAS_RATES calculation of rate coefficients for
!@+    gas phase sulfur oxidation chemistry 
!@auth Bell
!@ver  1.0
      USE MODEL_COM, only: im,jm,lm,t,ls1
      USE DYNAMICS, only: pmid,am,pk,LTROPO
      USE GEOM, only: dxyp,imaxj
      USE TRACER_COM, only: rsulf1,rsulf2,rsulf3,rsulf4 
#ifdef TRACERS_SPECIAL_Shindell
      USE TRCHEM_Shindell_COM, only: which_trop
#endif

      real*8 ppres,te,tt,mm,dmm,rk4,ek4
#ifdef TRACERS_SPECIAL_Shindell
!@var maxl chosen tropopause 0=LTROPO(I,J), 1=LS1-1
      integer maxl
#endif
 
C Initialise
      rsulf1(:,:,:)=0.0
      rsulf2(:,:,:)=0.0
      rsulf3(:,:,:)=0.0
      rsulf4(:,:,:)=0.0

C Reactions
C***1.DMS + OH -> 0.75SO2 + 0.25MSA
C***2.DMS + OH -> SO2
C***3.DMS + NO3 -> HNO3 + SO2
C***4.SO2 + OH -> SO4 + HO2

      do l=1,LM
      do j=1,jm
      do i=1,imaxj(j)
c
      maxl = ltropo(i,j)
#ifdef TRACERS_SPECIAL_Shindell
      if(which_trop.eq.1)maxl=ls1-1
#endif
      if(l.le.maxl) then

C Calculate effective temperature

        ppres=pmid(l,i,j)*9.869d-4 !in atm
        te=pk(l,i,j)*t(i,j,l)
        mm=am(l,i,j)*dxyp(j)
        tt = 1.d0/te

c DMM is number density of air in molecules/cm3
        dmm=ppres/(.082d0*te)*6.02d20

        rsulf1(i,j,l) = 1.7d-22*dmm*0.21d0*1.d-20*exp(7810.d0*tt)/
     *       (1.d0+5.5d-20*exp(7460.d0*tt)*dmm*0.21d0*1.d-11)
        
        rsulf2(i,j,l) = 9.6d-12*exp(-234.d0*tt)
        
        rsulf3(i,j,l) = 1.9d-13*exp(520.d0*tt)
        
        rk4 = 4.0d-20 *((tt*300.d0)**(3.3d0))*dmm*1.d-11
        ek4 = 1.d0/(1.d0 + ((log10(rk4/2.0d-12))**2.d0))
        
        rsulf4(i,j,l) = (rk4/(1.d0 + rk4/2.0d-12))*(0.45d0**ek4)

      endif

      end do
      end do
      end do
      
      END SUBROUTINE GET_SULF_GAS_RATES
#endif
