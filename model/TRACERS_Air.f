#include "rundeck_opts.h"

#ifdef TRACERS_ON
!@sum  TRACERS_DRV: tracer-dependent routines for air mass tracers
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
!@calls sync_param, SET_TCON
      USE CONSTANT, only: mair,mwat,sday
      USE MODEL_COM, only: dtsrc,byim
      USE DAGCOM, only: ia_src,ia_12hr,ir_log2
      USE TRACER_COM
      USE TRACER_DIAG_COM
      USE PARAM
#ifdef TRACERS_SPECIAL_Lerner
      USE TRACER_MPchem_COM, only: n_MPtable,tcscale
!@dbparam dsol describes portion of solar cycle being modeled for linoz
!@+      +1.0 = solar max, 0.0 = neutral, -1.0 = solar min
      USE LINOZ_CHEM_COM, only: dsol
#endif
#ifdef TRACERS_WATER
      USE LANDICE_COM, only : trli0    ! should these be in tracer_com?
      USE SEAICE_COM, only : trsi0
#endif      
      implicit none
      integer :: l,k,n
      character*20 sum_unit(ntm),inst_unit(ntm)   ! for conservation
      character*10 CMR
      logical :: qcon(KTCON-1), qsum(KTCON-1), T=.TRUE. , F=.FALSE.
      character*50 :: unit_string
#ifdef TRACERS_WATER
      real*8 fracls
#endif      

C**** Set defaults for tracer attributes (all dimensioned ntm)
      itime_tr0 = 0
      t_qlimit = .true.
      needtrs = .false.
      trdecay = 0.
      ntsurfsrc = 0
      trpdens = 0.
      trradius = 0.
#ifdef TRACERS_SPECIAL_Lerner
      n_MPtable = 0
      tcscale = 0.
#endif
#ifdef TRACERS_WATER
      tr_wd_TYPE = nGas         ! or  nPART  or nWATER
      tr_DHD = 0.
      tr_RKD = 0.
      trw0 = 0.
      trli0 = 0.
      trsi0 = 0.
      tr_H2ObyCH4 = 0.
#ifdef TRACERS_OCEAN
      trglac = 0.
#endif
#endif
C**** Define individual tracer characteristics
      do n=1,ntm
      select case (trname(n))
      
      case ('Air')
      n_Air = n
          itime_tr0(n) = 0.
          ntm_power(n) = -2
          tr_mm(n) = mair

      case ('SF6')
      n_SF6 = n
          itime_tr0(n) = 0.
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
          ntsurfsrc(n) = 1
#ifdef TRACERS_SPECIAL_Lerner
          n_MPtable(n) = 1
          tcscale(n_MPtable(n)) = 1.
#endif

      case ('CFC11')   !!! should start April 1
      n_CFC11 = n
          itime_tr0(n) = 0.
          ntm_power(n) = -12
          tr_mm(n) = 137.4d0
          ntsurfsrc(n) = 1
#ifdef TRACERS_SPECIAL_Lerner
          n_MPtable(n) = 2
          tcscale(n_MPtable(n)) = 1.
#endif

      case ('14CO2')   !!! should start 10/16
      n_14CO2 = n
          itime_tr0(n) = 0.
          ntm_power(n) = -18
          tr_mm(n) = 46.d0
          ntsurfsrc(n) = 1

      case ('CH4')
      n_CH4 = n
          ntm_power(n) = -9
          tr_mm(n) = 16.d0
          ntsurfsrc(n) = 14
#ifdef TRACERS_SPECIAL_Lerner
          n_MPtable(n) = 3
          tcscale(n_MPtable(n)) = 1.
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

#ifdef TRACERS_WATER
      case ('Water')
      n_Water = n
          ntm_power(n) = -4
          tr_mm(n) = mwat
          needtrs(n) = .true.
          tr_wd_TYPE(n) = nWater
          trw0(n) = 1.
          trli0(n) = 1.
          trsi0(n) = 1.
          tr_H2ObyCH4(n) = 1.
#ifdef TRACERS_OCEAN
          trglac(n) = 1.
#endif

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
#ifdef TRACERS_OCEAN
          trglac(n) = 0.
#endif
#endif
#endif

      end select
      end do

C**** Tracer sources and sinks
C**** Defaults for jls (sources, sinks, etc.)
C**** These need to be 'hand coded' depending on circumstances
      do k=1,ktajls  ! max number of sources and sinks
        jgrid_jls(k) = 1
        ia_jls(k) = ia_src
        scale_jls(k) = 1./DTsrc
      end do
      jls_index(:) = 0

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
        jls_index(k) = n
        jls_ltop(k) = 1
        jls_power(k) = -3.
        units_jls(k) = unit_string(jls_power(k),'kg/s')

      case ('Rn222')
        k = k + 1
        jls_decay(n) = k   ! special array for all radioactive sinks
        sname_jls(k) = 'Decay_of_'//trname(n)
        lname_jls(k) = 'LOSS OF RADON-222 BY DECAY'
        jls_index(k) = n
        jls_ltop(k) = lm
        jls_power(k) = -12.
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_source(1,n) = k
        sname_jls(k) = 'Ground_Source_of_'//trname(n)
        lname_jls(k) = 'RADON-222 SOURCE, LAYER 1'
        jls_index(k) = n
        jls_ltop(k) = 1
        jls_power(k) = -10.
        units_jls(k) = unit_string(jls_power(k),'kg/s')

! keep AIJ and AJL CO2 sources in same order !!
      case ('CO2')
        k = k + 1
        jls_source(1,n) = k
        sname_jls(k) = 'Fossil_fuel_source_'//trname(n)
        lname_jls(k) = 'CO2 Fossil fuel source (Marland)'
        jls_index(k) = n
        jls_ltop(k) = 1
        jls_power(k) = 3
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_source(2,n) = k
        sname_jls(k) = 'fertilization_sink_'//trname(n)
        lname_jls(k) = 'CO2 fertilization sink (Friedlingstein)'
        jls_index(k) = n
        jls_ltop(k) = 1
        jls_power(k) = 3
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_source(3,n) = k
        sname_jls(k) = 'Northern_forest_regrowth_'//trname(n)
        lname_jls(k) = 'CO2 Northern forest regrowth sink'
        jls_index(k) = n
        jls_ltop(k) = 1
        jls_power(k) = 3
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_source(4,n) = k
        sname_jls(k) = 'Land_Use_Modification_'//trname(n)
        lname_jls(k) = 'CO2 from Land use modification (Houton)'
        jls_index(k) = n
        jls_ltop(k) = 1
        jls_power(k) = 3
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_source(5,n) = k
        sname_jls(k) = 'Ecosystem_exchange_'//trname(n)
        lname_jls(k) = 'CO2 Ecosystem exchange (Matthews)'
        jls_index(k) = n
        jls_ltop(k) = 1
        jls_power(k) = 3
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_source(6,n) = k
        sname_jls(k) = 'Ocean_exchange_'//trname(n)
        lname_jls(k) = 'CO2 Ocean exchange'
        jls_index(k) = n
        jls_ltop(k) = 1
        jls_power(k) = 3
        units_jls(k) = unit_string(jls_power(k),'kg/s')

      case ('N2O')
        k = k + 1
        jls_source(1,n) = k
        sname_jls(k) = 'L1_sink_'//trname(n)
        lname_jls(k) = 'CHANGE OF N20 BY RESETTING TO 462.2d-9, L1'
        jls_index(k) = n
        jls_ltop(k) = 1
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_3Dsource(1,n) = k
        sname_jls(k) = 'Stratos_chem_change_'//trname(n)
        lname_jls(k) = 'CHANGE OF N2O BY CHEMISTRY IN STRATOS'
        jls_index(k) = n
        jls_ltop(k) = lm
        jls_power(k) = -1
        units_jls(k) = unit_string(jls_power(k),'kg/s')

      case ('CFC11')   !!! should start April 1
        k = k + 1
        jls_source(1,n) = k
        sname_jls(k) = 'L1_sink_'//trname(n)
        lname_jls(k) = 'CHANGE OF CFC-11 BY SOURCE, L1'
        jls_index(k) = n
        jls_ltop(k) = 1
        jls_power(k) = -1
        units_jls(k) = unit_string(jls_power(k),'kg/s')
        k = k + 1
        jls_3Dsource(1,n) = k
        sname_jls(k) = 'Stratos_chem_change_'//trname(n)
        lname_jls(k) = 'CHANGE OF CFC-11 BY CHEMISTRY IN STRATOS'
        jls_index(k) = n
        jls_ltop(k) = lm
        jls_power(k) = -3
        units_jls(k) = unit_string(jls_power(k),'kg/s')

      case ('14CO2')   !!! should start 10/16
        k = k + 1
        jls_source(1,n) = k
        sname_jls(k) = 'L1_sink_'//trname(n)
        lname_jls(k) = 'CHANGE OF 14CO2 by SINK, L1'
        jls_index(k) = n
        jls_ltop(k) = 1
        jls_power(k) = -4
        units_jls(k) = unit_string(jls_power(k),'kg/s')

      case ('CH4')
       k = k + 1
        jls_source(1,n) = k
        sname_jls(k) = 'Animal_source_'//trname(n)
        lname_jls(k) = 'CH4 Animal source'
        jls_index(k) = n
        jls_ltop(k) = 1
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),'kg/s')
       k = k + 1
        jls_source(2,n) = k
        sname_jls(k) = 'Coal_Mine_source_'//trname(n)
        lname_jls(k) = 'CH4 Coal Mine source'
        jls_index(k) = n
        jls_ltop(k) = 1
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),'kg/s')
       k = k + 1
        jls_source(3,n) = k
        sname_jls(k) = 'Gas_Leak_source_'//trname(n)
        lname_jls(k) = 'CH4 Gas Leak source'
        jls_index(k) = n
        jls_ltop(k) = 1
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),'kg/s')
       k = k + 1
        jls_source(4,n) = k
        sname_jls(k) = 'Gas_Venting_source_'//trname(n)
        lname_jls(k) = 'CH4 Gas Venting source'
        jls_index(k) = n
        jls_ltop(k) = 1
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),'kg/s')
       k = k + 1
        jls_source(5,n) = k
        sname_jls(k) = 'Municipal_solid_waste_source_'//trname(n)
        lname_jls(k) = 'CH4 Municipal solid waste source'
        jls_index(k) = n
        jls_ltop(k) = 1
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),'kg/s')
       k = k + 1
        jls_source(6,n) = k
        sname_jls(k) = 'Soil_sink_'//trname(n)
        lname_jls(k) = 'CH4 sink due to soil absorption'
        jls_index(k) = n
        jls_ltop(k) = 1
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),'kg/s')
       k = k + 1
        jls_source(7,n) = k
        sname_jls(k) = 'Termite_source_'//trname(n)
        lname_jls(k) = 'CH4 Termite source'
        jls_index(k) = n
        jls_ltop(k) = 1
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),'kg/s')
       k = k + 1
        jls_source(8,n) = k
        sname_jls(k) = 'Coal_combustion_source_'//trname(n)
        lname_jls(k) = 'CH4 Coal combustion source'
        jls_index(k) = n
        jls_ltop(k) = 1
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),'kg/s')
       k = k + 1
        jls_source(9,n) = k
        sname_jls(k) = 'Ocean_source_'//trname(n)
        lname_jls(k) = 'CH4 Ocean source'
        jls_index(k) = n
        jls_ltop(k) = 1
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),'kg/s')
       k = k + 1
        jls_source(10,n) = k
        sname_jls(k) = 'Fresh_Water_lake_source_'//trname(n)
        lname_jls(k) = 'CH4 Fresh Water lake source'
        jls_index(k) = n
        jls_ltop(k) = 1
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),'kg/s')
       k = k + 1
        jls_source(11,n) = k
        sname_jls(k) = 'Misc_Ground_source_'//trname(n)
        lname_jls(k) = 'CH4 Misc_Ground source'
        jls_index(k) = n
        jls_ltop(k) = 1
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),'kg/s')
       k = k + 1
        jls_source(12,n) = k
        sname_jls(k) = 'Biomass_burning_source_'//trname(n)
        lname_jls(k) = 'CH4 Biomass burning source'
        jls_index(k) = n
        jls_ltop(k) = 1
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),'kg/s')
       k = k + 1
        jls_source(13,n) = k
        sname_jls(k) = 'Rice Cultivation_source_'//trname(n)
        lname_jls(k) = 'CH4 Rice Cultivation source'
        jls_index(k) = n
        jls_ltop(k) = 1
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),'kg/s')
       k = k + 1
        jls_source(14,n) = k
        sname_jls(k) = 'Wetlands+Tundra_source_'//trname(n)
        lname_jls(k) = 'CH4 Wetlands+Tundra source'
        jls_index(k) = n
        jls_ltop(k) = 1
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),'kg/s')
       k = k + 1
        jls_3Dsource(2,n) = k
        sname_jls(k) = 'Stratos_Chem_change_'//trname(n)
        lname_jls(k) = 'CHANGE OF CH4 BY CHEMISTRY IN STRATOS'
        jls_index(k) = n
        jls_ltop(k) = lm
        jls_power(k) = -1
        units_jls(k) = unit_string(jls_power(k),' kg/s')
       k = k + 1
        jls_3Dsource(1,n) = k
        sname_jls(k) = 'Tropos_Chem_change_'//trname(n)
        lname_jls(k) = 'CHANGE OF CH4 BY CHEMISTRY IN TROPOSPHERE'
        jls_index(k) = n
        jls_ltop(k) = lm
        jls_power(k) = -1
        units_jls(k) = unit_string(jls_power(k),'kg/s')
       k = k + 1
        jls_3Dsource(3,n) = k
        sname_jls(k) = 'Total_Chem_change'//trname(n)
        lname_jls(k) = 'TOTAL CHANGE OF CH4 BY CHEMISTRY'
        jls_index(k) = n
        jls_ltop(k) = lm
        jls_power(k) = -1
        units_jls(k) = unit_string(jls_power(k),'kg/s')

      case ('O3')
       k = k + 1
        jls_3Dsource(2,n) = k
        sname_jls(k) = 'Stratos_Chem_change_'//trname(n)
        lname_jls(k) = 'CHANGE OF O3 BY CHEMISTRY IN STRATOS'
        jls_index(k) = n
        jls_ltop(k) = lm
        jls_power(k) = 1
        units_jls(k) = unit_string(jls_power(k),' kg/s')
       k = k + 1
        jls_3Dsource(1,n) = k
        sname_jls(k) = 'Tropos_Chem_change_'//trname(n)
        lname_jls(k) = 'CHANGE OF O3 BY CHEMISTRY IN TROPOSPHERE'
        jls_index(k) = n
        jls_ltop(k) = lm
        jls_power(k) = 1
        units_jls(k) = unit_string(jls_power(k),'kg/s')
       k = k + 1
        jls_3Dsource(3,n) = k
        sname_jls(k) = 'Total_Chem_change'//trname(n)
        lname_jls(k) = 'TOTAL CHANGE OF O3 BY CHEMISTRY'
        jls_index(k) = n
        jls_ltop(k) = lm
        jls_power(k) = 1
        units_jls(k) = unit_string(jls_power(k),'kg/s')

#ifdef TRACERS_WATER
C**** generic ones for many water tracers
      case ('Water', 'H2O18', 'HDO')
       k = k + 1
        jls_source(1,n)=k
        sname_jls(k) = 'Evap_'//trname(n)
        lname_jls(k) = 'EVAPORATION OF '//trname(n)
        jls_index(k) = n
        jls_ltop(k) = 1
        jls_power(k) = ntm_power(n)+4
        scale_jls(k) = SDAY*byim/DTsrc 
        units_jls(k) = unit_string(jls_power(k),'mm/day')
       k = k + 1
        jls_source(2,n)=k
        sname_jls(k) = 'Ocn_Evap_'//trname(n)
        lname_jls(k) = 'OCEAN EVAP OF '//trname(n)
        jls_index(k) = n
        jls_ltop(k) = 1
        jls_power(k) = ntm_power(n)+4
        scale_jls(k) = SDAY*byim/DTsrc 
        units_jls(k) = unit_string(jls_power(k),'mm/day')
       k = k + 1
        jls_source(3,n)=k
        sname_jls(k) = 'Precip_'//trname(n)
        lname_jls(k) = 'PRECIPITATION OF '//trname(n)
        jls_index(k) = n
        jls_ltop(k) = 1
        jls_power(k) = ntm_power(n)+4
        scale_jls(k) = SDAY*byim/DTsrc 
        units_jls(k) = unit_string(jls_power(k),'mm/day')
       k = k + 1
        jls_source(4,n)=k
        sname_jls(k) = 'Ocn_Precip_'//trname(n)
        lname_jls(k) = 'OCEAN PRECIP OF '//trname(n)
        jls_index(k) = n
        jls_ltop(k) = 1
        jls_power(k) = ntm_power(n)+4
        scale_jls(k) = SDAY*byim/DTsrc 
        units_jls(k) = unit_string(jls_power(k),'mm/day')

C**** special unique to HTO
      case ('HTO')
       k = k + 1
        jls_source(1,n)=k
        sname_jls(k) = 'Evap_'//trname(n)
        lname_jls(k) = 'EVAPORATION OF '//trname(n)
        jls_index(k) = n
        jls_ltop(k) = 1
        jls_power(k) = ntm_power(n)+4
        scale_jls(k) = SDAY*byim/DTsrc 
        units_jls(k) = unit_string(jls_power(k),'mm/day')
       k = k + 1
        jls_source(2,n)=k
        sname_jls(k) = 'Ocn_Evap_'//trname(n)
        lname_jls(k) = 'OCEAN EVAP OF '//trname(n)
        jls_index(k) = n
        jls_ltop(k) = 1
        jls_power(k) = ntm_power(n)+4
        scale_jls(k) = SDAY*byim/DTsrc 
        units_jls(k) = unit_string(jls_power(k),'mm/day')
       k = k + 1
        jls_source(3,n)=k
        sname_jls(k) = 'Precip_'//trname(n)
        lname_jls(k) = 'PRECIPITATION OF '//trname(n)
        jls_index(k) = n
        jls_ltop(k) = 1
        jls_power(k) = ntm_power(n)+4
        scale_jls(k) = SDAY*byim/DTsrc 
        units_jls(k) = unit_string(jls_power(k),'mm/day')
       k = k + 1
        jls_source(4,n)=k
        sname_jls(k) = 'Ocn_Precip_'//trname(n)
        lname_jls(k) = 'OCEAN PRECIP OF '//trname(n)
        jls_index(k) = n
        jls_ltop(k) = 1
        jls_power(k) = ntm_power(n)+4
        scale_jls(k) = SDAY*byim/DTsrc 
        units_jls(k) = unit_string(jls_power(k),'mm/day')
       k = k + 1
        jls_decay(n) = k   ! special array for all radioactive sinks
        sname_jls(k) = 'Decay_of_'//trname(n)
        lname_jls(k) = 'LOSS OF '//TRIM(trname(n))//' BY DECAY'
        jls_index(k) = n
        jls_ltop(k) = lm
        jls_power(k) = ntm_power(n)+8
        scale_jls(k) = 1./DTsrc 
        units_jls(k) = unit_string(jls_power(k),'kg/s')
#endif

C**** Here are some more examples of generalised diag. configuration
c      n = n_dust
c        k = k + 1 
c        jls_grav(n) = k   ! special array grav. settling sinks
c        sname_jls(k) = 'Grav_Settle_of_'//trname(n)
c        lname_jls(k) = 'LOSS OF DUST BY SETTLING'
c        jls_index(k) = n
c        jls_ltop(k) = lm
c        jls_power(k) = -11.
c        units_jls(k) = unit_string(jls_power(k),'kg/s')

C**** Checks
      if (ntsurfsrc(n).gt.ntsurfsrcmax) then
!       write(6,*) ' ntsurfsrc too large for ',trname(n)
        write(6,*) ' Increase ntsurfsrcmax to at least',ntsurfsrc(n)
        STOP ' Ntsurfsrc too large.  Increase ntsurfsrcmax'
      end if

      end select
      end do

      if (k.gt. ktajls) then
        write (6,*) 
     &   'tjl_defs: Increase ktajls=',ktajls,' to at least ',k
        stop 'ktajls too small'
      end if


C**** Tracer sources and sinks
C**** Defaults for ijts (sources, sinks, etc.)
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
        sname_ijts(k) = 'SF6_CFC-GRID_SOURCE,_LAYER_1'
        ijts_power(k) = -15
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc

      case ('Rn222')
      k = k+1
        ijts_source(1,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src   
        lname_ijts(k) = 'Rn222 L 1 SOURCE'
        sname_ijts(k) = 'Radon-222_SOURCE,_Layer_1'
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
      k = k + 1
        ijts_source(1,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'N2O CHANGE IN L 1'
        sname_ijts(k) = 'N2O_CHANGE_IN_L_1'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc

      case ('CFC11')
      k = k + 1
        ijts_source(1,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CFC_11 L 1 SOURCE'
        sname_ijts(k) = 'CFC_11_SOURCE,_LAYER_1'
        ijts_power(k) = -15
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

      case ('CH4')
      k = k + 1
        ijts_source(1,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 Animal source'
        sname_ijts(k) = 'CH4_Animal_source'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(2,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 Coal Mine source'
        sname_ijts(k) = 'CH4_Coal Mine source'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(3,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 Gas Leak source'
        sname_ijts(k) = 'CH4_Gas Leak source'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(4,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 Gas Venting source'
        sname_ijts(k) = 'CH4_Gas Venting source'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(5,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 Municipal solid waste src'
        sname_ijts(k) = 'CH4_Municipal solid waste src'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(6,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 sink due to soil absorp.'
        sname_ijts(k) = 'CH4_sink due to soil absorp.'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(7,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 Termite source'
        sname_ijts(k) = 'CH4_Termite source'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(8,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 Coal combustion source'
        sname_ijts(k) = 'CH4_Coal combustion source'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(9,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 Ocean source'
        sname_ijts(k) = 'CH4_Ocean_source'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(10,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 Fresh Water lake source'
        sname_ijts(k) = 'CH4_Fresh Water lake source'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(11,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 Misc. Ground source'
        sname_ijts(k) = 'CH4_Misc. Ground source'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(12,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 Biomass burning source'
        sname_ijts(k) = 'CH4_Biomass burning source'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(13,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 Rice cultivation source'
        sname_ijts(k) = 'CH4_Rice cultivation source'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(14,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 Wetlands+Tundra source'
        sname_ijts(k) = 'CH4_Wetlands+Tundra source'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),'kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc

#ifdef TRACERS_WATER
        case ('Water')
          ! nothing I can think of....
        case ('H2O18')
        case ('HDO')
        case ('HTO')
#endif

      end select
      end do

      if (k .gt. ktaijs) then
        write (6,*)'ijt_defs: Increase ktaijs=',ktaijs,' to at least ',k
        stop 'ktaijs too small'
      end if

C**** Initialize conservation diagnostics
C**** To add a new conservation diagnostic:
C****       Set up a QCON, and call SET_TCON to allocate array numbers,
C****       set up scales, titles, etc. 
C**** QCON denotes when the conservation diags should be accumulated
C**** QSUM says whether that diag is to be used in summation (if the
C****      routine DIAGCTB is used, this must be false).
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
      itcon_surf(1,N) = 13
      qcon(itcon_surf(1,N)) = .true.; conpts(1) = 'Reset in L1'
      itcon_3Dsrc(1,N) = 14
      qcon(itcon_3Dsrc(1,N)) = .true.; conpts(2) = 'Strat. Chem'
      qsum(itcon_3Dsrc(1,N)) = .true.
      CALL SET_TCON(QCON,TRNAME(N),QSUM,inst_unit(n),
     *     sum_unit(n),scale_inst(n),scale_change(n), N,CONPTs)
      qsum(13:) = .false.  ! reset to defaults for next tracer

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

#ifdef TRACERS_WATER
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
#endif

C**** Here are some more examples of conservation diag configuration
C**** Gravitional settling:
c      n=n_dust
c      itcon_grav(n) = xx
c      qcon(itcon_grav(n)) = .true.; conpts(yy) = 'SETTLING'
c      qsum(itcon_grav(n)) = .true.
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
      WRITE (6,'(A14,2I8)') "KTACC=",KTACC

      return
      end subroutine init_tracer


      SUBROUTINE tracer_IC
!@sum tracer_IC initializes tracers when they are first switched on
!@auth Jean Lerner
      USE CONSTANT, only: mair,rhow
      USE MODEL_COM, only: itime,im,jm,lm,ls1
#ifdef TRACERS_WATER
     *  ,q,wm,flice,fearth
#endif
      USE TRACER_COM, only: ntm,trm,trmom,itime_tr0,trname,needtrs,
     *   tr_mm
#ifdef TRACERS_WATER
     *  ,trwm,trw0
      USE SOMTQ_COM, only : qmom
      USE LANDICE, only : ace1li,ace2li
      USE LANDICE_COM, only : trli0,trsnowli,trlndi,snowli
      USE SEAICE, only : xsi,ace1i
      USE SEAICE_COM, only : rsi,msi,snowi,trsi,trsi0,ssi
      USE LAKES_COM, only : trlake,mwl,mldlk,flake
      USE GHYCOM, only : trbare,trvege,trsnowbv,wbare,wvege,snowbv
      USE FLUXES, only : gtracer
#endif
      USE GEOM, only: dxyp,bydxyp
      USE DYNAMICS, only: am,byam  ! Air mass of each box (kg/m^2)
      USE PBLCOM, only: npbl,trabl
#ifdef TRACERS_WATER
     *  ,qabl
#endif
#ifdef TRACERS_SPECIAL_Lerner
      USE LINOZ_CHEM_COM, only: tlt0m,tltzm, tltzzm
      USE PRATHER_CHEM_COM, only: nstrtc
      USE QUSDEF, only : mz,mzz
#endif
      USE FILEMANAGER, only: openunit,closeunit
      IMPLICIT NONE
      INTEGER i,n,l,j,iu_data,ipbl,it,lr
      CHARACTER*80 title
      REAL*8 CFC11ic,ic14CO2(im,jm,lm),conv
      REAL*8 :: trinit =1., tmominit=0.
      REAL*4 CO2ic(im,jm,lm),N2Oic(jm,lm)
      REAL*8 CH4ic(jm,lm)
      EQUIVALENCE (CO2ic,N2Oic,ic14CO2,CH4ic)

      do n=1,ntm

      if (itime.eq.itime_tr0(n)) then

#ifdef TRACERS_WATER
C**** set default atmospheric liquid water amount to zero for most tracers        
      trwm(:,:,:,n)=0.
#endif
      select case (trname(n)) 
          
        case default            
c          write(6,*) 'In TRACER_IC:',trname(n),' does not exist '
          stop "TRACER_IC"
  
        case ('Air')
          do l=1,lm
          do j=1,jm
            trm(:,j,l,n) = am(l,:,j)*dxyp(j)
          end do; enddo
          do i=2,im
            trm(i,1,:,n) =  trm(1,1,:,n) !poles
            trm(i,jm,:,n) = trm(1,jm,:,n) !poles
          enddo
          trmom(:,:,:,:,n) = 0.

        case ('SF6')
          trm(:,:,:,n) = 0.
          trmom(:,:,:,:,n) = 0.

        case ('Rn222')
          do l=1,lm
            do j=1,jm
              trm(:,j,l,n) = am(l,:,j)*dxyp(j)*1.d-22
          end do; end do
          do i=2,im
            trm(i,1,:,n) =  trm(1,1,:,n)   !poles
            trm(i,jm,:,n) = trm(1,jm,:,n)   !poles
          enddo
          trmom(:,:,:,:,n) = 0.
 
        case ('CO2')
          trm(:,:,:,n) = 0.
          trmom(:,:,:,:,n) = 0.
          call openunit('CO2_IC',iu_data,.true.,.true.)
          read (iu_data) title,co2ic
          call closeunit(iu_data)
          write(6,*) title,' read from CO2_IC'
          do l=1,lm         !ppmv==>ppmm
          do j=1,jm
            trm(:,j,l,n) = co2ic(:,j,l)*am(l,:,j)*dxyp(j)*1.54d-6
          enddo; enddo

        case ('N2O')
          trmom(:,:,:,:,n) = 0.
          call openunit('N2O_IC',iu_data,.true.,.true.)
          read (iu_data) title,N2Oic     ! unit is PPMM/(M*DXYP)
          call closeunit(iu_data)
          write(6,*) title,' read from N2O_IC'
          do l=1,lm         !ppmv==>ppmm
          do j=1,jm
            trm(:,j,l,n) = am(l,:,j)*dxyp(j)*N2Oic(j,l)
          enddo; enddo

        case ('CFC11')   !!! should start April 1
          trmom(:,:,:,:,n) = 0.
          CFC11ic = 268.D-12*136.5/29.029    !268 PPTV
          do l=1,lm
          do j=1,jm
            trm(:,j,l,n) = am(l,:,j)*dxyp(j)*CFC11ic
          enddo; enddo

        case ('14CO2')   !!! this tracer is supposed to start 10/16
          trmom(:,:,:,:,n) = 0.
#ifdef TRACERS_SPECIAL_Lerner
          call get_14CO2_IC(ic14CO2)
          do l=1,lm         !ppmv==>ppmm
          do j=1,jm
            trm(:,j,l,n) = am(l,:,j)*dxyp(j)*ic14CO2(:,j,l)*1.d-18
          enddo; enddo
#endif

        case ('CH4')
          trmom(:,:,:,:,n) = 0.
#ifdef TRACERS_SPECIAL_Lerner
          call get_wofsy_gas_IC(trname(n),CH4ic)
          do l=1,lm         !ppbv==>ppbm
          do j=1,jm                          
            trm(:,j,l,n) = am(l,:,j)*dxyp(j)*CH4ic(j,l)*0.552d-9
          enddo; enddo
#endif

        case ('O3')
          trmom(:,:,:,:,n) = 0.
          do l=1,lm
          do j=1,jm
            trm(:,j,l,n) = am(l,:,j)*dxyp(j)*20.d-9*tr_mm(n)/mair
          enddo; enddo
#ifdef TRACERS_SPECIAL_Lerner
          do l=lm,lm+1-nstrtc,-1
          lr = lm+1-l
            do j=1,jm
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
          tmominit = 0. 
        case ('HDO')   ! dD=-650
          trinit=0.35d0*trw0(n)
          tmominit = 0. 
        case ('HTO')
          trinit=0.
          tmominit = 0. 
        end select

        do l=1,lm
        do j=1,jm
          trm(:,j,l,n) =  q(:,j,l)*am(l,:,j)*dxyp(j)*trinit
          trwm(:,j,l,n)= wm(:,j,l)*am(l,:,j)*dxyp(j)*trinit
          do i=1,im
            trmom(:,i,j,l,n) = qmom(:,i,j,l)*am(l,i,j)*dxyp(j)*tmominit
          end do
        end do
        end do
        do i=2,im
          trm(i, 1,:,n) = trm(1, 1,:,n) !poles
          trm(i,jm,:,n) = trm(1,jm,:,n) !poles
          trwm(i, 1,:,n)= trwm(1, 1,:,n) !poles
          trwm(i,jm,:,n)= trwm(1,jm,:,n) !poles
          trmom(:,i, 1,:,n)=0.
          trmom(:,i,jm,:,n)=0.
        enddo
        if (trname(n).eq."HTO") then ! initialise bomb source
          do l=ls1-1,ls1+1      ! strat. source 
            do j=35,37          ! lat 44 N - 56 N
              trm(:,j,l,n)= q(:,j,l)*am(l,:,j)*dxyp(j)*1d10*1d-18
            end do
          end do
        end if
              
        do j=1,jm
          do i=1,im
C**** lakes
            if (flake(i,j).gt.0) then
              trlake(n,1,i,j)=trw0(n)*mldlk(i,j)*rhow*flake(i,j)*dxyp(j)
              trlake(n,2,i,j)=trw0(n)*mwl(i,j)-trlake(n,1,i,j)
              gtracer(n,1,i,j)=trw0(n)
            end if
c**** ice
            if (rsi(i,j).gt.0) then
              trsi(n,1,i,j)=trsi0(n)*
     *             (xsi(1)*(snowi(i,j)+ace1i)-ssi(1,i,j))
              trsi(n,2,i,j)=trsi0(n)*
     *             (xsi(2)*(snowi(i,j)+ace1i)-ssi(2,i,j))
              trsi(n,3,i,j)=trsi0(n)*(xsi(3)*msi(i,j)-ssi(3,i,j))
              trsi(n,4,i,j)=trsi0(n)*(xsi(4)*msi(i,j)-ssi(4,i,j))
              gtracer(n,2,i,j)=trsi0(n)
            else
              gtracer(n,2,i,j)=0.
            end if
c**** landice
            if (flice(i,j).gt.0) then
              trlndi(n,i,j)=trli0(n)*(ace1li+ace2li)
              trsnowli(n,i,j)=trli0(n)*snowli(i,j)
              gtracer(n,3,i,j)=trli0(n)
            else
              gtracer(n,3,i,j)=0.
            end if
c**** earth
            if (fearth(i,j).gt.0) then
              conv=rhow         ! convert from m to kg/m^2
              trbare  (n,:,i,j)=trw0(n)*wbare (:,i,j)*conv
              trvege  (n,:,i,j)=trw0(n)*wvege (:,i,j)*conv
              trsnowbv(n,1,i,j)=trw0(n)*snowbv(1,i,j)*conv
              trsnowbv(n,2,i,j)=trw0(n)*snowbv(2,i,j)*conv
              gtracer (n,4,i,j)=trw0(n)
            else
              gtracer(n,4,i,j)=0.
            end if
          end do
        end do
#endif
        
      end select

C**** Initialise pbl profile if necessary
      if (needtrs(n)) then
        do it=1,4
        do j=1,jm
        do ipbl=1,npbl
#ifndef TRACERS_WATER
          trabl(ipbl,n,:,j,it) = trm(:,j,1,n)*byam(1,:,j)*bydxyp(j)
#else
          trabl(ipbl,n,:,j,it) = trinit*qabl(ipbl,:,j,it)
#endif
        end do
        end do
        end do
      end if

      write(6,*) ' Tracer ',trname(n),' initialized at itime=',itime
      end if
      end do
#ifdef TRACERS_WATER
C**** Initialise ocean tracers if necessary
      call tracer_ic_ocean
#endif
C****
      end subroutine tracer_IC


      subroutine daily_tracer(iact)
!@sum daily_tracer is called once a day for tracers
!@+   SUBROUTINE tracer_IC is called from daily_tracer to allow for
!@+     tracers that 'turn on' on different dates.
!@auth Jean Lerner
C**** Note this routine must always exist (but can be a dummy routine)
      USE MODEL_COM, only: jmon,itime
      USE TRACER_COM, only: ntm,trname,itime_tr0
#ifdef TRACERS_SPECIAL_Lerner
      USE TRACER_MPchem_COM, only: n_MPtable,tcscale,STRATCHEM_SETUP
      USE LINOZ_CHEM_COM, only: LINOZ_SETUP
#endif
      IMPLICIT NONE
      INTEGER n,iact,last_month
      data last_month/-1/

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

C****
C**** Initialize tracers here to allow for tracers that 'turn on'
C****  at any time
      call tracer_IC

      return
      end subroutine daily_tracer


      SUBROUTINE set_tracer_2Dsource
!@sum tracer_source calculates non-interactive sources for tracers
!@auth Jean Lerner/Gavin Schmidt
      USE MODEL_COM, only: FEARTH,itime,JDperY,fland,psf,pmtop,jmpery
     *  ,dtsrc
      USE GEOM, only: dxyp,areag
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
      implicit none
      integer :: i,j,ns,l,ky,n
      REAL*8 :: source,sarea,steppy,base,steppd,x,airm,anngas,
     *  steph,stepx,stepp,tmon,bydt,tnew

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
C**** Source over United States and Canada
        source = .37d0*anngas*steppy
        sarea  = 0.
        do j=31,35
          do i=12,22
            sarea = sarea + dxyp(j)*fearth(i,j)
          enddo
        enddo
        do j=31,35
          do i=12,22
            trsource(i,j,1,n) = source*dxyp(j)*fearth(i,j)/sarea
          enddo
        enddo
C**** Source over Europe and Russia
        source = .37d0*anngas*steppy
        sarea  = 0.
        do j=33,39
          do i=35,45
            sarea = sarea + dxyp(j)*fearth(i,j)
          enddo
        enddo
        do j=33,39
          do i=35,45
            trsource(i,j,1,n) = source*dxyp(j)*fearth(i,j)/sarea
          enddo
        enddo
C**** Source over Far East
        source = .13d0*anngas*steppy
        sarea  = 0.
        do j=29,34
          do i=61,66
            sarea = sarea + dxyp(j)*fearth(i,j)
          enddo
        enddo
        do j=29,34
          do i=61,66
            trsource(i,j,1,n) = source*dxyp(j)*fearth(i,j)/sarea
          enddo
        enddo
C**** Source over Middle East
        source = .05d0*anngas*steppy
        sarea  = 0.
        do j=28,32
          do i=43,51
            sarea = sarea + dxyp(j)*fearth(i,j)
          enddo
        enddo
        do j=28,32
          do i=43,51
            trsource(i,j,1,n) = source*dxyp(j)*fearth(i,j)/sarea
          enddo
        enddo
C**** Source over South America
        source = .04d0*anngas*steppy
        i=27; j=18
        trsource(i,j,1,n) = 0.5*source
        i=28; j=18
        trsource(i,j,1,n) = 0.5*source
C**** Source over South Africa
        source = .02d0*anngas*steppy
        i=42; j=17
        trsource(i,j,1,n) = source
C**** Source over Australia and New Zealand
        source = .02d0*anngas*steppy
        i=66; j=15
        trsource(i,j,1,n) = source

C****
C**** Surface Sources for Radon-222
C****
      case ('Rn222')
        trsource(:,:,:,n)=0
C**** ground source
        steppd = 1./sday
        do j=1,jm
          do i=1,im
C**** source from ice-free land
            if(tsavg(i,j).lt.tf) then !composite surface air temperature
              trsource(i,j,1,n) = 1.0d-16*steppd*dxyp(j)*fearth(i,j)
            else
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
          do j=1,jm
            trsource(:,j,ns,n) = co2_src(:,j,ns)*dxyp(j)
          end do
        end do

C****
C**** Sources and sinks for CH4 (kg/s)
C****
      case ('CH4')
        do ns=1,ntsurfsrc(n)
          do j=1,jm
            trsource(:,j,ns,n) = ch4_src(:,j,ns)*dxyp(j)
          end do
        end do
#endif
C****
C**** Sources and sinks for N2O:
C**** First layer is set to a constant 462.2 ppbm. (300 PPB V)
C****
      case ('N2O')
      do j=1,jm
        trsource(:,j,1,n) = (am(1,:,j)*dxyp(j)*462.2d-9
     *   -trm(:,j,1,n))*bydt
      end do

C****
C**** Sources and sinks for 14C02
C**** NOTE: This tracer is supposed to start on 10/16
C**** The tracer is reset to specific values in layer 1 only if
C****   this results in a sink
C****
      case ('14CO2')
      tmon = (itime-itime_tr0(n))*jmpery/(hrday*jdpery)  !(12./8760.)
      trsource(:,:,1,n) = 0.
      do j=1,jm/2
      do i=1,im
        tnew = am(1,i,j)*dxyp(j)*(4.82d-18*46./mair)*
     *   (44.5 + tmon*(1.02535d0 - tmon*(2.13565d-2 - tmon*8.61853d-5)))
        if (tnew.lt.trm(i,j,1,n)) 
     *         trsource(i,j,1,n) = (tnew-trm(i,j,1,n))*bydt
      end do
      end do
      do j=1+jm/2,jm
      do i=1,im
        tnew = am(1,i,j)*dxyp(j)*(4.82d-18*46./mair)*
     *   (73.0 - tmon*(0.27823d0 + tmon*(3.45648d-3 - tmon*4.21159d-5)))
        if (tnew.lt.trm(i,j,1,n)) 
     *         trsource(i,j,1,n) = (tnew-trm(i,j,1,n))*bydt
      end do
      end do

C****
C**** No non-interactive surface sources of Water
C****
      case ('Water')
        trsource(:,:,:,n)=0

      end select

      end do
C****
      END SUBROUTINE set_tracer_2Dsource


      SUBROUTINE tracer_3Dsource
!@sum tracer_3Dsource calculates interactive sources for tracers
!@auth Jean Lerner
      USE MODEL_COM, only: itime
      USE TRACER_COM
      USE TRACER_DIAG_COM, only: itcon_3Dsrc,tajls,jls_3Dsource
      USE FLUXES, only: tr3Dsource
      implicit none
      integer n

C**** All sources are saved as kg/s
      do n=1,ntm

      if (itime.lt.itime_tr0(n)) cycle
      select case (trname(n))

      case default
#ifdef TRACERS_SPECIAL_Lerner
C****
      case ('CH4')
      tr3Dsource(:,:,:,:,n) = 0.
      call Trop_chem_CH4(1,n)
      call apply_tracer_3Dsource(1,n)
      call Strat_chem_Prather(2,n)
      call DIAGTCA(itcon_3Dsrc(2,n),n)
      tajls(:,:,jls_3Dsource(3,n)) = 
     *  tajls(:,:,jls_3Dsource(1,n))+tajls(:,:,jls_3Dsource(2,n))
C****
      case ('O3')
      tr3Dsource(:,:,:,:,n) = 0.
      call Trop_chem_O3(1,n)
   !  call apply_tracer_3Dsource(1,n)
      call DIAGTCA(itcon_3Dsrc(1,n),n)
      call Strat_chem_O3(2,n)
   !  call apply_tracer_3Dsource(2,n)
      call DIAGTCA(itcon_3Dsrc(2,n),n)
      tajls(:,:,jls_3Dsource(3,n)) = 
     *  tajls(:,:,jls_3Dsource(1,n))+tajls(:,:,jls_3Dsource(2,n))
C****
      case ('N2O')
      tr3Dsource(:,:,:,:,n) = 0.
      call Strat_chem_Prather(1,n)
      call DIAGTCA(itcon_3Dsrc(1,n),n)
C****
      case ('CFC11')
      tr3Dsource(:,:,:,:,n) = 0.
      call Strat_chem_Prather(1,n)
      call DIAGTCA(itcon_3Dsrc(1,n),n)
C****
#endif
      end select

      end do
      return
      END SUBROUTINE tracer_3Dsource

#endif


#ifdef TRACERS_WATER
C---SUBROUTINES FOR TRACER WET DEPOSITION----------------------------

      SUBROUTINE GET_COND_FACTOR(L,N,WMXTR,TEMP,FCLOUD,FQ0,fq)
!@sum  GET_COND_FACTOR calculation of condensate fraction for tracers
!@+    within or below convective or large-scale clouds. Gas 
!@+    condensation uses Henry's Law if not freezing.
!@auth Dorothy Koch (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on CB436TdsM23 CLOUDCHCC and CLOUDCHEM subroutines)
c
C**** GLOBAL parameters and variables:
      USE CONSTANT, only: TF, BYGASC, MAIR,teeny
      USE TRACER_COM, only: tr_RKD,tr_DHD,nWATER,nGAS,nPART,tr_wd_TYPE
     *     ,trname
      USE CLOUDS, only: PL, NTIX
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
!@var FQ0 [default] fraction of tracer that goes into condensate
!@var L index for altitude loop
!@var N index for tracer number loop
!@var WMXTR mixing ratio of water available for tracer condensation?
!@var RKD dummy variable (= tr_RKD*EXP[ ])
      REAL*8, PARAMETER :: BY298K=3.3557D-3
      REAL*8 Ppas, tfac, ssfac, RKD
#ifdef TRACERS_SPECIAL_O18
      real*8 tdegc,alph,fracvs,fracvl
#endif
      REAL*8,  INTENT(IN) :: fq0, FCLOUD, WMXTR, TEMP
      REAL*8,  INTENT(OUT):: fq
      INTEGER, INTENT(IN) :: L, N
c
C**** CALCULATE the fraction of tracer mass that becomes condensate:
c
      SELECT CASE(tr_wd_TYPE(NTIX(N)))                 
        CASE(nGAS)                            ! gas tracer
          fq = 0.D0                           ! frozen case
          IF(TEMP.ge.TF) THEN                ! if not frozen then: 
            Ppas = PL(L)*1.D2                 ! pressure to pascals
            tfac = (1.D0/TEMP - BY298K)*BYGASC
            IF(tr_DHD(NTIX(N)).ne.0.D0) THEN
              RKD=tr_RKD(NTIX(N))*DEXP(-tr_DHD(NTIX(N))*tfac)
            ELSE  
              RKD=tr_RKD(NTIX(N))
            END IF
c           clwc=WMXTR*MAIR*1.D-3*Ppas*BYGASC/(TEMP*FCLOUD)
c           ssfac=RKD*GASC*TEMP*clwc   ! Henry's Law
            ssfac=RKD*WMXTR*MAIR*1.D-3*Ppas/(FCLOUD+teeny)
            fq=ssfac / (1.D0 + ssfac)
          END IF
        CASE(nWATER)                          ! water tracer
#ifdef TRACERS_SPECIAL_O18
C**** calculate condensate in equilibrium with source vapour
          if (fq0.gt.0.) then
            tdegc=temp-tf
            if (tdegc.ge.0) then
              alph=1./fracvl(tdegc,trname(ntix(n)))
            else
              alph=1./fracvs(tdegc,trname(ntix(n)))
            end if
            if (fq0.ne.1.) then ! just to be safe
              fq = alph * fq0/(1.+(alph-1.)*fq0)
            else
              fq = fq0
            end if
          else
            fq = 0.
          end if
#else
          fq = fq0                                  
#endif
        CASE(nPART)                           ! particulate tracer
          fq = 0.D0                           ! temporarily zero.
c NOTE 1: Dorothy has some code that will be put here to condense 
c aerosols. GSF 1/4/02
c
c NOTE 2:  Really, any aerosol 'formation', (meaning the production of
c an aerosol tracer due to cloud chemistry, or any flux among tracers),
c should be done elsewhere, like in a chemistry section of the model.
c But if it is impossible to supply that section with the variables
c needed from the wet deposition code, then the aerosol formation code
c should probably go here... If you add this, please make appropriate
c changes in the subroutine's name/summary above. GSF 1/4/02. 
        CASE DEFAULT                                ! error
          STOP 'tr_wd_TYPE(NTIX(N)) out of range in SCAVENGE_TRACER'
      END SELECT
c      
      RETURN
      END SUBROUTINE GET_COND_FACTOR


      SUBROUTINE GET_PREC_FACTOR(N,BELOW_CLOUD,CM,FCLD,FQ0,fq)
!@sum  GET_PREC_FACTOR calculation of the precipitation scavenging
!@+    fraction for tracers WITHIN large scale clouds. Current version
!@+    uses the first order removal rate based on [Giorgi and
!@+    Chameides, 1986], for gaseous and particulate tracers.
!@auth Dorothy Koch (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on CB436TdsM23 RAINOUT subroutine)
c
C**** GLOBAL parameters and variables:
      USE MODEL_COM, only: dtsrc
      USE TRACER_COM, only: nWATER, nGAS, nPART, tr_wd_TYPE
      USE CLOUDS, only: NTIX
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
      INTEGER, INTENT(IN) :: N
      REAL*8,  INTENT(IN) :: FQ0, FCLD, CM
      REAL*8,  INTENT(OUT):: FQ
c
      SELECT CASE(tr_wd_TYPE(NTIX(N)))                 
        CASE(nGAS)                                ! gas
          IF(BELOW_CLOUD) THEN
            fq = 0.D0
          ELSE
c           minus preserves FPRT sign convention in LSCOND
            fq = -(CM*DTsrc*(EXP(-CM*DTsrc)- 1.0)*FCLD)
          END IF
        CASE(nWATER)                          ! water/original method
          fq = FQ0                            ! no fractionation
        CASE(nPART)                           ! aerosols
          IF(BELOW_CLOUD) THEN
            fq = 0.D0
          ELSE
c           minus preserves FPRT sign convention in LSCOND
            fq = -(CM*DTsrc*(EXP(-CM*DTsrc)- 1.0)*FCLD)
          END IF
        CASE DEFAULT                          ! error
          STOP 'tr_wd_TYPE(NTIX(N)) out of range in GET_FPRT'
      END SELECT
c
      RETURN
      END SUBROUTINE GET_PREC_FACTOR


      SUBROUTINE GET_WASH_FACTOR(N,b_beta_DT,PREC,fq)
!@sum  GET_WASH_FACTOR calculation of the fraction of tracer 
!@+    scavanged by precipitation below convective clouds ("washout"). 
!@auth Dorothy Koch (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on CB436TdsM23 CWASH and WASH_EVAP routines)
c
C**** GLOBAL parameters and variables:
      USE TRACER_COM, only: nWATER, nGAS, nPART, tr_wd_TYPE
      USE CLOUDS, only: NTIX
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
      INTEGER, INTENT(IN) :: N
      REAL*8, INTENT(OUT):: FQ
      REAL*8, INTENT(IN) :: PREC, b_beta_DT
      REAL*8, PARAMETER :: rc_wash = 1.D-1
C
      SELECT CASE(tr_wd_TYPE(NTIX(N)))                 
        CASE(nGAS)                            ! gas
          fq = 0.D0
        CASE(nWATER)                          ! water/original method
          fq = 0.D0                              
        CASE(nPART)                           ! aerosols
          fq = b_beta_DT*(DEXP(-PREC*rc_wash)-1.)
        CASE DEFAULT                          ! error
          STOP 'tr_wd_TYPE(NTIX(N)) out of range in WASHOUT_TRACER'
      END SELECT   
c      
      RETURN
      END SUBROUTINE GET_WASH_FACTOR


      SUBROUTINE GET_EVAP_FACTOR(N,TEMP,FQ0,fq)
!@sum  GET_EVAP_FACTOR calculation of the evaporation fraction
!@+    for tracers.
!@auth Dorothy Koch (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on CB436TdsM23 EVAPD and WASH_EVAP routines)
c
C**** GLOBAL parameters and variables:
      USE CONSTANT, only : tf
      USE TRACER_COM, only: tr_evap_fact, tr_wd_TYPE,nwater,trname
      USE CLOUDS, only: NTIX
c
      IMPLICIT NONE
c      
C**** Local parameters and variables and arguments:
!@var FQ            fraction of tracer evaporated
!@var FQ0 [default] fraction of tracer evaporated
!@var N index for tracer number loop
      INTEGER, INTENT(IN) :: N
      REAL*8,  INTENT(OUT):: FQ
      REAL*8,  INTENT(IN) :: FQ0,TEMP
#ifdef TRACERS_SPECIAL_O18
      real*8 tdegc,alph,fracvl,fracvs
#endif
c
      select case (tr_wd_TYPE(NTIX(N)))
      case default
        fq=FQ0*tr_evap_fact(tr_wd_TYPE(NTIX(N)))
c
      case (nWater)
#ifdef TRACERS_SPECIAL_O18
          tdegc=temp-tf
          if (tdegc.ge.0) then
            alph=fracvl(tdegc,trname(ntix(n)))
          else
            alph=fracvs(tdegc,trname(ntix(n)))
          end if
          if (fq0.ne.1.) then 
c           if (fq0.lt.0.9) then ! approximate 
              fq = alph * fq0
c           else ! calculate actual rayleigh curve (necessary?)
c             fq = 1. - (1.-fq0)**alph
c           end if
          else
            fq = fq0
          end if
#else
          fq=FQ0*tr_evap_fact(tr_wd_TYPE(NTIX(N)))
#endif
      end select
      RETURN
      END SUBROUTINE GET_EVAP_FACTOR 
#endif
