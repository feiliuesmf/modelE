#include "rundeck_opts.h"

#ifdef TRACERS_ON
!@sum  TRACERS_Air: tracer-dependent routines for air mass tracers
!@+    Routines included:
!@+      Those that MUST EXIST for all tracers: 
!@+        Diagnostic specs: init_tracer
!@+        Tracer initialisation + sources: tracer_ic, set_tracer_source
!@+        Entry points: daily_tracer
!@+      Those that are unique for specific tracers: 
!@+        read_co2_sources, read_ch4_sources
!@auth Jean Lerner/Gavin Schmidt
!@ver  1.0

      SUBROUTINE init_tracer
!@sum init_tracer initializes trace gas attributes and diagnostics
!@auth J. Lerner
!@calls sync_param, SET_TCON
      use DAGCOM, only: ia_src,ia_12hr,ir_log2
      USE MODEL_COM, only: dtsrc
      USE TRACER_COM
      USE TRACER_DIAG_COM
      USE CONSTANT, only: mair,sday
      USE PARAM
      implicit none
      integer :: l,k,n
      character*20 sum_unit(ntm),inst_unit(ntm)   ! for conservation
      character*10 CMR
      logical :: qcon(KTCON-1), qsum(KTCON-1), T=.TRUE. , F=.FALSE.
!@var to_volume_MixRat: For printout of tracer concentration
!@+   to_volume_MixRat=1: printout is in Volume Mixing Ratio
!@+   to_volume_MixRat=0: printout is in Mass Mixing Ratio
      INTEGER :: to_volume_MixRat=0 
      character*50 :: unit_string

C**** Tracer sources and sinks
C**** Defaults for jls (sources, sinks, etc.)
C**** These need to be 'hand coded' depending on circumstances
      do k=1,ktajls  ! max number of sources and sinks
        jgrid_jls(k) = 1
        ia_jls(k) = ia_src
        scale_jls(k) = 1./DTsrc
      end do
      jls_index(:) = 0

C**** Define individual tracer characteristics
      do n=1,ntm
      select case (trname(n))
      
      case ('Air')
      n_Air = n
          itime_tr0(n) = 0.
          ntm_power(n) = -2
          tr_mm(n) = mair
          t_qlimit(n) = .true.
          needtrs(n) = .false.
          trdecay(n) = 0.d0
          nt3Dsrc(n) = 0;  ntsurfsrc(n) = 0
          trpdens(n) = 0;  trradius(n) = 0

      case ('SF6')
      n_SF6 = n
          itime_tr0(n) = 0.
          ntm_power(n) = -14
          tr_mm(n) = 146.01d0
          t_qlimit(n) = .true.
          needtrs(n) = .false.
          trdecay(n) = 0.d0
          nt3Dsrc(n) = 0;  ntsurfsrc(n) = 1
          trpdens(n) = 0;  trradius(n) = 0

      case ('Rn222')
      n_Rn222 = n
          itime_tr0(n) = 0.
          ntm_power(n) = -21
          tr_mm(n) = 222.d0
          t_qlimit(n) = .true.
          needtrs(n) = .false.
          trdecay(n) = 2.1d-6
          nt3Dsrc(n) = 0;  ntsurfsrc(n) = 1
          trpdens(n) = 0;  trradius(n) = 0

      case ('CO2')
      n_CO2 = n
          itime_tr0(n) = 0.
          ntm_power(n) = -6
          tr_mm(n) = 44.d0
          t_qlimit(n) = .false.
          needtrs(n) = .false.
          trdecay(n) = 0.d0
          nt3Dsrc(n) = 0;  ntsurfsrc(n) = 6
          trpdens(n) = 0;  trradius(n) = 0

      case ('N2O')
      n_N2O = n
          itime_tr0(n) = 0.
          ntm_power(n) = -9
          tr_mm(n) = 44.d0
          t_qlimit(n) = .true.
          needtrs(n) = .false.
          trdecay(n) = 0.d0
          nt3Dsrc(n) = 1;  ntsurfsrc(n) = 1
          trpdens(n) = 0;  trradius(n) = 0

      case ('CFC11')   !!! should start April 1
      n_CFC11 = n
          itime_tr0(n) = 0.
          ntm_power(n) = -12
          tr_mm(n) = 137.4d0
          t_qlimit(n) = .true.
          needtrs(n) = .false.
          trdecay(n) = 0.d0
          nt3Dsrc(n) = 1;  ntsurfsrc(n) = 1
          trpdens(n) = 0;  trradius(n) = 0

      case ('14CO2')   !!! should start 10/16
      n_14CO2 = n
          itime_tr0(n) = 0.
          ntm_power(n) = -18
          tr_mm(n) = 46.d0
          t_qlimit(n) = .true.
          needtrs(n) = .false.
          trdecay(n) = 0.d0
          nt3Dsrc(n) = 0;  ntsurfsrc(n) = 1
          trpdens(n) = 0;  trradius(n) = 0
      case ('CH4')
      n_CH4 = n
          itime_tr0(n) = 0.
          ntm_power(n) = -9
          tr_mm(n) = 16.d0
          t_qlimit(n) = .true.
          needtrs(n) = .false.
          trdecay(n) = 0.d0
          nt3Dsrc(n) = 2;  ntsurfsrc(n) = 14
          trpdens(n) = 0;  trradius(n) = 0
      case ('O3')
      n_O3 = n
          itime_tr0(n) = 0.
          ntm_power(n) = -8
          tr_mm(n) = 48.d0
          t_qlimit(n) = .true.
          needtrs(n) = .false.
          trdecay(n) = 0.d0
          nt3Dsrc(n) = 2;  ntsurfsrc(n) = 0
          trpdens(n) = 0;  trradius(n) = 0

      end select
      end do

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
        units_jls(k) = unit_string(jls_power(k),' kg/s')

      case ('Rn222')
        k = k + 1
        jls_decay(n) = k   ! special array for all radioactive sinks
        sname_jls(k) = 'Decay_of_'//trname(n)
        lname_jls(k) = 'LOSS OF RADON-222 BY DECAY'
        jls_index(k) = n
        jls_ltop(k) = lm
        jls_power(k) = -11.
        units_jls(k) = unit_string(jls_power(k),' kg/s')
        k = k + 1
        jls_source(1,n) = k
        sname_jls(k) = 'Ground_Source_of_'//trname(n)
        lname_jls(k) = 'RADON-222 SOURCE, LAYER 1'
        jls_index(k) = n
        jls_ltop(k) = 1
        jls_power(k) = -10.
        units_jls(k) = unit_string(jls_power(k),' kg/s')

! keep AIJ and AJL CO2 sources in same order !!
      case ('CO2')
        k = k + 1
        jls_source(1,n) = k
        sname_jls(k) = 'Fossil_fuel_source_'//trname(n)
        lname_jls(k) = 'CO2 Fossil fuel source (Marland)'
        jls_index(k) = n
        jls_ltop(k) = 1
        jls_power(k) = 3
        units_jls(k) = unit_string(jls_power(k),' kg/s')
        k = k + 1
        jls_source(2,n) = k
        sname_jls(k) = 'fertilization_sink_'//trname(n)
        lname_jls(k) = 'CO2 fertilization sink (Friedlingstein)'
        jls_index(k) = n
        jls_ltop(k) = 1
        jls_power(k) = 3
        units_jls(k) = unit_string(jls_power(k),' kg/s')
        k = k + 1
        jls_source(3,n) = k
        sname_jls(k) = 'Northern_forest_regrowth_'//trname(n)
        lname_jls(k) = 'CO2 Northern forest regrowth sink'
        jls_index(k) = n
        jls_ltop(k) = 1
        jls_power(k) = 3
        units_jls(k) = unit_string(jls_power(k),' kg/s')
        k = k + 1
        jls_source(4,n) = k
        sname_jls(k) = 'Land_Use_Modification_'//trname(n)
        lname_jls(k) = 'CO2 from Land use modification (Houton)'
        jls_index(k) = n
        jls_ltop(k) = 1
        jls_power(k) = 3
        units_jls(k) = unit_string(jls_power(k),' kg/s')
        k = k + 1
        jls_source(5,n) = k
        sname_jls(k) = 'Ecosystem_exchange_'//trname(n)
        lname_jls(k) = 'CO2 Ecosystem exchange (Matthews)'
        jls_index(k) = n
        jls_ltop(k) = 1
        jls_power(k) = 3
        units_jls(k) = unit_string(jls_power(k),' kg/s')
        k = k + 1
        jls_source(6,n) = k
        sname_jls(k) = 'Ocean_exchange_'//trname(n)
        lname_jls(k) = 'CO2 Ocean exchange'
        jls_index(k) = n
        jls_ltop(k) = 1
        jls_power(k) = 3
        units_jls(k) = unit_string(jls_power(k),' kg/s')

      case ('N2O')
        k = k + 1
        jls_source(1,n) = k
        sname_jls(k) = 'L1_sink_'//trname(n)
        lname_jls(k) = 'CHANGE OF N20 BY RESETTING TO 462.2d-9, L1'
        jls_index(k) = n
        jls_ltop(k) = 1
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),' kg/s')
        k = k + 1
        jls_3Dsource(1,n) = k
        sname_jls(k) = 'Stratos_chem_change_'//trname(n)
        lname_jls(k) = 'CHANGE OF N2O BY CHEMISTRY IN STRATOS'
        jls_index(k) = n
        jls_ltop(k) = lm
        jls_power(k) = -1
        units_jls(k) = unit_string(jls_power(k),' kg/s')

      case ('CFC11')   !!! should start April 1
        k = k + 1
        jls_source(1,n) = k
        sname_jls(k) = 'L1_sink_'//trname(n)
        lname_jls(k) = 'CHANGE OF CFC-11 BY SOURCE, L1'
        jls_index(k) = n
        jls_ltop(k) = 1
        jls_power(k) = -1
        units_jls(k) = unit_string(jls_power(k),' kg/s')
        k = k + 1
        jls_3Dsource(1,n) = k
        sname_jls(k) = 'Stratos_chem_change_'//trname(n)
        lname_jls(k) = 'CHANGE OF CFC-11 BY CHEMISTRY IN STRATOS'
        jls_index(k) = n
        jls_ltop(k) = lm
        jls_power(k) = -14
        units_jls(k) = unit_string(jls_power(k),' kg/s')

      case ('14CO2')   !!! should start 10/16
        k = k + 1
        jls_source(1,n) = k
        sname_jls(k) = 'L1_sink_'//trname(n)
        lname_jls(k) = 'CHANGE OF 14CO2 by SINK, L1'
        jls_index(k) = n
        jls_ltop(k) = 1
        jls_power(k) = -4
        units_jls(k) = unit_string(jls_power(k),' kg/s')

      case ('CH4')
       k = k + 1
        jls_source(1,n) = k
        sname_jls(k) = 'Animal_source_'//trname(n)
        lname_jls(k) = 'CH4 Animal source'
        jls_index(k) = n
        jls_ltop(k) = 1
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),' kg/s')
       k = k + 1
        jls_source(2,n) = k
        sname_jls(k) = 'Coal_Mine_source_'//trname(n)
        lname_jls(k) = 'CH4 Coal Mine source'
        jls_index(k) = n
        jls_ltop(k) = 1
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),' kg/s')
       k = k + 1
        jls_source(3,n) = k
        sname_jls(k) = 'Gas_Leak_source_'//trname(n)
        lname_jls(k) = 'CH4 Gas Leak source'
        jls_index(k) = n
        jls_ltop(k) = 1
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),' kg/s')
       k = k + 1
        jls_source(4,n) = k
        sname_jls(k) = 'Gas_Venting_source_'//trname(n)
        lname_jls(k) = 'CH4 Gas Venting source'
        jls_index(k) = n
        jls_ltop(k) = 1
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),' kg/s')
       k = k + 1
        jls_source(5,n) = k
        sname_jls(k) = 'Municipal_solid_waste_source_'//trname(n)
        lname_jls(k) = 'CH4 Municipal solid waste source'
        jls_index(k) = n
        jls_ltop(k) = 1
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),' kg/s')
       k = k + 1
        jls_source(6,n) = k
        sname_jls(k) = 'Soil_sink_'//trname(n)
        lname_jls(k) = 'CH4 sink due to soil absorption'
        jls_index(k) = n
        jls_ltop(k) = 1
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),' kg/s')
       k = k + 1
        jls_source(7,n) = k
        sname_jls(k) = 'Termite_source_'//trname(n)
        lname_jls(k) = 'CH4 Termite source'
        jls_index(k) = n
        jls_ltop(k) = 1
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),' kg/s')
       k = k + 1
        jls_source(8,n) = k
        sname_jls(k) = 'Coal_combustion_source_'//trname(n)
        lname_jls(k) = 'CH4 Coal combustion source'
        jls_index(k) = n
        jls_ltop(k) = 1
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),' kg/s')
       k = k + 1
        jls_source(9,n) = k
        sname_jls(k) = 'Ocean_source_'//trname(n)
        lname_jls(k) = 'CH4 Ocean source'
        jls_index(k) = n
        jls_ltop(k) = 1
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),' kg/s')
       k = k + 1
        jls_source(10,n) = k
        sname_jls(k) = 'Fresh_Water_lake_source_'//trname(n)
        lname_jls(k) = 'CH4 Fresh Water lake source'
        jls_index(k) = n
        jls_ltop(k) = 1
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),' kg/s')
       k = k + 1
        jls_source(11,n) = k
        sname_jls(k) = 'Misc_Ground_source_'//trname(n)
        lname_jls(k) = 'CH4 Misc_Ground source'
        jls_index(k) = n
        jls_ltop(k) = 1
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),' kg/s')
       k = k + 1
        jls_source(12,n) = k
        sname_jls(k) = 'Biomass_burning_source_'//trname(n)
        lname_jls(k) = 'CH4 Biomass burning source'
        jls_index(k) = n
        jls_ltop(k) = 1
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),' kg/s')
       k = k + 1
        jls_source(13,n) = k
        sname_jls(k) = 'Rice Cultivation_source_'//trname(n)
        lname_jls(k) = 'CH4 Rice Cultivation source'
        jls_index(k) = n
        jls_ltop(k) = 1
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),' kg/s')
       k = k + 1
        jls_source(14,n) = k
        sname_jls(k) = 'Wetlands+Tundra_source_'//trname(n)
        lname_jls(k) = 'CH4 Wetlands+Tundra source'
        jls_index(k) = n
        jls_ltop(k) = 1
        jls_power(k) = 0
        units_jls(k) = unit_string(jls_power(k),' kg/s')
       k = k + 1
        jls_3Dsource(1,n) = k
        sname_jls(k) = 'Tropos_Chem_change_'//trname(n)
        lname_jls(k) = 'CHANGE OF CH4 BY CHEMISTRY IN TROPOSHPERE'
        jls_index(k) = n
        jls_ltop(k) = lm
        jls_power(k) = -1
        units_jls(k) = unit_string(jls_power(k),' kg/s')
       k = k + 1
        jls_3Dsource(2,n) = k
        sname_jls(k) = 'Stratos_Chem_change_'//trname(n)
        lname_jls(k) = 'CHANGE OF CH4 BY CHEMISTRY IN STRATOS'
        jls_index(k) = n
        jls_ltop(k) = lm
        jls_power(k) = -1
        units_jls(k) = unit_string(jls_power(k),' kg/s')

      case ('O3')
       k = k + 1
        jls_3Dsource(1,n) = k
        sname_jls(k) = 'Tropos_Chem_change_'//trname(n)
        lname_jls(k) = 'CHANGE OF O3 BY CHEMISTRY IN TROPOSHPERE'
        jls_index(k) = n
        jls_ltop(k) = lm
        jls_power(k) = 1
        units_jls(k) = unit_string(jls_power(k),' kg/s')
       k = k + 1
        jls_3Dsource(2,n) = k
        sname_jls(k) = 'Stratos_Chem_change_'//trname(n)
        lname_jls(k) = 'CHANGE OF O3 BY CHEMISTRY IN STRATOS'
        jls_index(k) = n
        jls_ltop(k) = lm
        jls_power(k) = 1
        units_jls(k) = unit_string(jls_power(k),' kg/s')

C**** Here are some more examples of generalised diag. configuration
c      n = n_dust
c        k = k + 1 
c        jls_grav(n) = k   ! special array grav. settling sinks
c        sname_jls(k) = 'Grav_Settle_of_'//trname(n)
c        lname_jls(k) = 'LOSS OF DUST BY SETTLING'
c        jls_index(k) = n
c        jls_ltop(k) = lm
c        jls_power(k) = -11.
c        units_jls(k) = unit_string(jls_power(k),' kg/s')

C**** Checks
      if (ntsurfsrc(n).gt.ntsurfsrcmax) then
!       write(6,*) ' ntsurfsrc too large for ',trname(n)
        write(6,*) ' Increase ntsurfsrcmax to at least',ntsurfsrc(n)
        STOP ' Ntsurfsrc too large.  Increase ntsurfsrcmax'
      end if
      if (nt3Dsrc(n).gt.nt3Dsrcmax) then
!       write(6,*) ' nt3Dsrc too large for ',trname(n)
        write(6,*) ' Increase nt3Dsrcmax to at least',nt3Dsrc(n)
        STOP ' nt3Dsrc too large.  Increase nt3Dsrcmax'
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
        units_ijts(k) = unit_string(ijts_power(k),' kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc

      case ('Rn222')
      k = k+1
        ijts_source(1,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src   
        lname_ijts(k) = 'Rn222 L 1 SOURCE'
        sname_ijts(k) = 'Radon-222_SOURCE,_Layer_1'
        ijts_power(k) = -21.
        units_ijts(k) = unit_string(ijts_power(k),' kg/s*m^2')
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
        units_ijts(k) = unit_string(ijts_power(k),' kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1 
        ijts_source(2,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src   
        sname_ijts(k) = 'fertilization_sink_'//trname(n)
        lname_ijts(k) = 'CO2 fertilization'
        ijts_power(k) = -11.
        units_ijts(k) = unit_string(ijts_power(k),' kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1 
        ijts_source(3,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src   
        sname_ijts(k) = 'Northern_forest_regrowth_'//trname(n)
        lname_ijts(k) = 'CO2 North forest regrowth'
        ijts_power(k) = -11.
        units_ijts(k) = unit_string(ijts_power(k),' kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1 
        ijts_source(4,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src   
        sname_ijts(k) = 'Land_Use_Modification_'//trname(n)
        lname_ijts(k) = 'CO2 from Land use mods'
        ijts_power(k) = -11.
        units_ijts(k) = unit_string(ijts_power(k),' kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1 
        ijts_source(5,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src   
        sname_ijts(k) = 'Ecosystem_exchange_'//trname(n)
        lname_ijts(k) = 'CO2 Ecosystem exch'
        ijts_power(k) = -11.
        units_ijts(k) = unit_string(ijts_power(k),' kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1 
        ijts_source(6,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src   
        sname_ijts(k) = 'Ocean_exchange_'//trname(n)
        lname_ijts(k) = 'CO2 Ocean exchange'
        ijts_power(k) = -11.
        units_ijts(k) = unit_string(ijts_power(k),' kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc

      case ('N2O')
      k = k + 1
        ijts_source(1,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'N2O CHANGE IN L 1'
        sname_ijts(k) = 'N2O_CHANGE_IN_L_1'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),' kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc

      case ('CFC11')
      k = k + 1
        ijts_source(1,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CFC_11 L 1 SOURCE'
        sname_ijts(k) = 'CFC_11_SOURCE,_LAYER_1'
        ijts_power(k) = -15
        units_ijts(k) = unit_string(ijts_power(k),' kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc

      case ('14CO2')
      k = k + 1
        ijts_source(1,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = '14CO2 L 1 Sink'
        sname_ijts(k) = '14CO2_L1_Sink'
        ijts_power(k) = -21
        units_ijts(k) = unit_string(ijts_power(k),' kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc

      case ('CH4')
      k = k + 1
        ijts_source(1,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 Animal source'
        sname_ijts(k) = 'CH4_Animal_source'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),' kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(2,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 Coal Mine source'
        sname_ijts(k) = 'CH4_Coal Mine source'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),' kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(3,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 Gas Leak source'
        sname_ijts(k) = 'CH4_Gas Leak source'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),' kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(4,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 Gas Venting source'
        sname_ijts(k) = 'CH4_Gas Venting source'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),' kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(5,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 Municipal solid waste src'
        sname_ijts(k) = 'CH4_Municipal solid waste src'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),' kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(6,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 sink due to soil absorp.'
        sname_ijts(k) = 'CH4_sink due to soil absorp.'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),' kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(7,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 Termite source'
        sname_ijts(k) = 'CH4_Termite source'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),' kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(8,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 Coal combustion source'
        sname_ijts(k) = 'CH4_Coal combustion source'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),' kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(9,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 Ocean source'
        sname_ijts(k) = 'CH4_Ocean_source'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),' kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(10,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 Fresh Water lake source'
        sname_ijts(k) = 'CH4_Fresh Water lake source'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),' kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(11,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 Misc. Ground source'
        sname_ijts(k) = 'CH4_Misc. Ground source'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),' kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(12,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 Biomass burning source'
        sname_ijts(k) = 'CH4_Biomass burning source'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),' kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(13,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 Rice cultivation source'
        sname_ijts(k) = 'CH4_Rice cultivation source'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),' kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1
        ijts_source(14,n) = k
        ijts_index(k) = n
        ia_ijts(k) = ia_src
        lname_ijts(k) = 'CH4 Wetlands+Tundra source'
        sname_ijts(k) = 'CH4_Wetlands+Tundra source'
        ijts_power(k) = -12
        units_ijts(k) = unit_string(ijts_power(k),' kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc

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
        inst_unit(n) = unit_string(kt_power_inst(n),  ' kg/m^2)')
         sum_unit(n) = unit_string(kt_power_change(n),' kg/m^2 s)')
      end do

      N = n_air
      CALL SET_TCON(QCON,TRNAME(N),QSUM,inst_unit(n),
     *     sum_unit(n),scale_inst(n),scale_change(n), N,CONPTs)

      N = n_SF6
      CALL SET_TCON(QCON,TRNAME(N),QSUM,inst_unit(n),
     *     sum_unit(n),scale_inst(n),scale_change(n), N,CONPTs)

      N = n_Rn222
      itcon_decay(n) = 13
      qcon(itcon_decay(n)) = .true.; conpts(1) = 'DECAY'
      qsum(itcon_decay(n)) = .true.
      CALL SET_TCON(QCON,TRNAME(N),QSUM,inst_unit(n),
     *     sum_unit(n),scale_inst(n),scale_change(n), N,CONPTs)

      N = n_CO2
      itcon_surf(1,N) = 13
      qcon(itcon_surf(1,N)) = .true.; conpts(1) = 'FossilFuel'
!     qsum(itcon_surf(1,N)) = .false.
      itcon_surf(2,N) = 14
      qcon(itcon_surf(2,N)) = .true.; conpts(2) = 'Fertilization'
!     qsum(itcon_surf(2,N)) = .false.
      itcon_surf(3,N) = 15
      qcon(itcon_surf(3,N)) = .true.; conpts(3) = 'Forest Regrowth'
!     qsum(itcon_surf(3,N)) = .false.
      itcon_surf(4,N) = 16
      qcon(itcon_surf(4,N)) = .true.; conpts(4) = 'Land Use'
!     qsum(itcon_surf(4,N)) = .false.
      itcon_surf(5,N) = 17
      qcon(itcon_surf(5,N)) = .true.; conpts(5) = 'Ecosystem Exch'
!     qsum(itcon_surf(5,N)) = .false.
      itcon_surf(6,N) = 18
      qcon(itcon_surf(6,N)) = .true.; conpts(6) = 'Ocean Exch'
!     qsum(itcon_surf(6,N)) = .false.
      qsum(itcon_surf(1:6,N)) = .false.
      CALL SET_TCON(QCON,TRNAME(N),QSUM,inst_unit(n),
     *     sum_unit(n),scale_inst(n),scale_change(n), N,CONPTs)
      qcon(13:) = .false.  ! reset to defaults for next tracer
      qsum(13:) = .false.  ! reset to defaults for next tracer

      N = n_N2O
      itcon_surf(1,N) = 13
      qcon(itcon_surf(1,N)) = .true.; conpts(1) = 'Reset in L1'
      itcon_3Dsrc(1,N) = 14
      qcon(itcon_3Dsrc(1,N)) = .true.; conpts(2) = 'Strat. Chem'
      qsum(itcon_3Dsrc(1,N)) = .true.
      CALL SET_TCON(QCON,TRNAME(N),QSUM,inst_unit(n),
     *     sum_unit(n),scale_inst(n),scale_change(n), N,CONPTs)

      N = n_CFC11
      itcon_surf(1,N) = 13
      qcon(itcon_surf(1,N)) = .true.; conpts(1) = 'L1 Source'
      itcon_3Dsrc(1,N) = 14
      qcon(itcon_surf(1,N)) = .true.; conpts(2) = 'Strat. Chem'
      qsum(itcon_3Dsrc(1,N)) = .true.
      CALL SET_TCON(QCON,TRNAME(N),QSUM,inst_unit(n),
     *     sum_unit(n),scale_inst(n),scale_change(n), N,CONPTs)

      N = n_14CO2
      itcon_surf(1,N) = 13
      qcon(itcon_surf(1,N)) = .true.; conpts(1) = 'Observed drift'
      CALL SET_TCON(QCON,TRNAME(N),QSUM,inst_unit(n),
     *     sum_unit(n),scale_inst(n),scale_change(n), N,CONPTs)

      N = n_CH4
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
      qsum(itcon_surf(1:14,N)) = .false.
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

      N = n_O3
      itcon_3Dsrc(1,N) = 1
      qcon(itcon_3Dsrc(1,N)) = .true.; conpts(1) = 'Tropos. Chem'
      itcon_3Dsrc(2,N) = 2
      qcon(itcon_3Dsrc(2,N)) = .true.; conpts(2) = 'Stratos. Chem'
      CALL SET_TCON(QCON,TRNAME(N),QSUM,inst_unit(n),
     *     sum_unit(n),scale_inst(n),scale_change(n), N,CONPTs)

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

C**** print out total tracer diagnostic array size
      WRITE (6,'(A14,2I8)') "KTACC=",KTACC

      return
      end subroutine init_tracer


      SUBROUTINE tracer_IC
!@sum tracer_IC initializes tracers when they are first switched on
!@auth Jean Lerner
      USE MODEL_COM, only: itime,im,jm,lm,ls1
      USE GEOM, only: dxyp,bydxyp
      USE CONSTANT, only: mair
      USE DYNAMICS, only: am,byam  ! Air mass of each box (kg/m^2)
      USE PBLCOM, only: npbl,trabl
      USE TRACER_COM, only: ntm,trm,trmom,itime_tr0,trname,needtrs,
     *   tr_mm
      USE FILEMANAGER, only: openunit,closeunit
      IMPLICIT NONE
      INTEGER i,n,l,j,iu_data,ipbl,it,lr
      CHARACTER*80 title
      REAL*8 CFC11ic,ic14CO2(im,jm,lm)
      REAL*4 CO2ic(im,jm,lm),N2Oic(jm,lm),CH4ic(jm,lm)
      EQUIVALENCE (CO2ic,N2Oic,ic14CO2,CH4ic)

      do n=1,ntm

      if (itime.eq.itime_tr0(n)) then
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
          call get_14CO2_IC(ic14CO2)
          do l=1,lm         !ppmv==>ppmm
          do j=1,jm
            trm(:,j,l,n) = am(l,:,j)*dxyp(j)*ic14CO2(:,j,l)*1.d-18
          enddo; enddo

        case ('CH4')
          trmom(:,:,:,:,n) = 0.
          call openunit('CH4_IC',iu_data,.true.,.true.)
          read (iu_data) title,CH4ic
          call closeunit(iu_data)
          write(6,*) title,' read from CH4_ic'
          do l=1,lm         !ppbv==>ppbm
          do j=1,jm                          
            trm(:,j,l,n) = am(l,:,j)*dxyp(j)*CH4ic(j,l)*0.552d-9
          enddo; enddo

        case ('O3')
          trmom(:,:,:,:,n) = 0.
          do l=1,lm
          do j=1,jm
            trm(:,j,l,n) = am(l,:,j)*dxyp(j)*20.d-9*tr_mm(n)/mair
          enddo; enddo
          do l=lm,ls1+1,-1
          lr = lm+1-l
c           do j=1,jm
c           trm(:,j,l,n) = 
c    *          tlt0m(j,lr,1)*am(l,:,j)*dxyp(j)*tr_mm(n)/mair
c           trmom(mz,:,j,l,n)  =
c    *          tltzm(j,lr,1)*am(l,:,j)*dxyp(j)*tr_mm(n)/mair
c           trmom(mzz,:,j,l,n)  =
c    *         tltzzm(j,lr,1)*am(l,:,j)*dxyp(j)*tr_mm(n)/mair
c           end do
          end do

      end select

C**** Initialise pbl profile if necessary
      if (needtrs(n)) then
        do it=1,4
        do j=1,jm
        do ipbl=1,npbl
          trabl(ipbl,n,:,j,it) = trm(:,j,1,n)*byam(1,:,j)*bydxyp(j)
        end do
        end do
        end do
      end if

      write(6,*) ' Tracer ',trname(n),' initialized at itime=',itime
      end if
      end do
C****
      end subroutine tracer_IC


      subroutine daily_tracer(iact)
!@sum daily_tracer is called once a day for tracers
!@auth Jean Lerner
C**** Note this routine must always exist (but can be a dummy routine)
      USE TRACER_COM, only: ntm,trname
      IMPLICIT NONE
      INTEGER n,iact

C**** Initialize tracers here to allow for tracers that 'turn on'
C****  at any time
      call tracer_IC

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

      return
      end subroutine daily_tracer


      MODULE CH4_SOURCES
      USE TRACER_COM
!@var CH4_src CH4 surface sources and sinks (kg/s)
      integer, parameter :: nch4src=14
      real*8 CH4_src(im,jm,nch4src)
      END MODULE CH4_SOURCES


      subroutine read_CH4_sources(nt,iact)
!@sum reads in CH4 sources and sinks
!@auth Jean Lerner
C****
C**** There are 3 monthly sources and 11 annual sources
C**** Annual sources are read in at start and re-start of run only
C**** Monthly sources are interpolated each day
      USE MODEL_COM, only: itime,JDperY,im,jm,jday,focean,fearth
      USE LAKES_COM, only: flake
      USE CONSTANT, only: sday
      USE FILEMANAGER, only: openunit,closeunit, openunits,closeunits
      USE TRACER_COM, only: itime_tr0,trname
      use CH4_SOURCES, only: src=>ch4_src,nsrc=>nch4src
      implicit none
      character*80 title
      logical :: ifirst=.true.
      real*8 adj(nsrc)
      data adj/1.3847,1.0285,3.904,1.659,1.233,1.194,0.999,
     *  7.2154,3.7247d-5,3.1399d-4,5.4838d-5,
     *  0.4369,0.7533,0.9818/
!@var nanns,nmons: number of annual and monthly input files
      integer, parameter :: nanns=11,nmons=3
      integer ann_units(nanns-3),mon_units(nmons-3)
      character*12 :: ann_files(nanns-3) =
     *  (/'CH4_ANIMALS ','CH4_COALMINE','CH4_GASLEAK ','CH4_GASVENT ',
     *    'CH4_CITYDUMP','CH4_SOIL_ABS','CH4_TERMITES','CH4_COALBURN'/)
      logical :: ann_bins(nanns-3)=(/.true.,.true.,.true.,.true.,
     *        .true.,.true.,.true.,.true./)
      character*8 :: mon_files(nmons) =
     *   (/'CH4_BURN','CH4_RICE','CH4_WETL'/)
      real*8 adj_wet(jm)
      data adj_wet/15*0.6585,16*1.761,15*0.6585/ !zonal adj FOR WETLANDS
      data iwet/14/  !!! position of wetlands array in src
      logical :: mon_bins(nmons)=(/.true.,.true.,.true./)
      real*8 tlca(im,jm,nmons),tlcb(im,jm,nmons)  ! for monthly sources
      integer i,j,nt,iact,iu,k,jdlast,iwet
      data jdlast /0/
      save ifirst,jdlast,tlca,tlcb

      if (itime.lt.itime_tr0(nt)) return
C****
C**** Annual Sources and sinks
C**** Apply adjustment factors to bring sources into balance
C**** Annual sources are in KG C/M2/Y
C**** Sources need to be kg/m^2 s; convert /year to /s
C****
      if (ifirst) then
        k = 0
        call openunits(ann_files,ann_units,ann_bins,nanns-3)
        ann_units(nanns-1:nanns) = ann_units(nanns-3) 
        do iu = ann_units(1),ann_units(nanns-3)
          k = k+1
          call readt (iu,0,src(1,1,k),im*jm,src(1,1,k),1)
          src(:,:,k) = src(:,:,k)*adj(k)/(sday*JDperY)
        end do
        call closeunits(ann_units,nanns-3)
        ! 3 miscellaneous sources
          k = k+1
          src(:,:,k) = focean(:,:)*adj(k)/(sday*JDperY)
          k = k+1
          src(:,:,k) =  flake(:,:)*adj(k)/(sday*JDperY)
          k = k+1
          src(:,:,k) = fearth(:,:)*adj(k)/(sday*JDperY)
      endif
C****
C**** Monthly sources are interpolated to the current day
C****
C**** Also, Apply adjustment factors to bring sources into balance
C**** Monthly sources are in KG C/M2/S => src in kg/m^2 s
C****
      if (iact.eq.1 .and. .not.ifirst) return
      ifirst = .false.
      call openunits(mon_files,mon_units,mon_bins,nmons)
      j = 0
      do k = nanns+1,nsrc
        j = j+1
        call read_monthly_sources(mon_units(j),jdlast,
     *    tlca(1,1,j),tlcb(1,1,j),src(1,1,k))
        src(:,:,k) = src(:,:,k)*adj(k)
      end do
      jdlast = jday

      call closeunits(mon_units,nmons)
      write(6,*) trname(nt),'Sources interpolated to current day'
C****
C**** Zonal adjustment for combined wetlands and tundra
C****
      do j=1,jm
        src(:,j,iwet) = src(:,j,iwet)*adj_wet(j)
      end do
      return
      end subroutine read_CH4_sources


      MODULE CO2_SOURCES
      USE TRACER_COM
!@var co2_src C02 surface sources and sinks (kg/s)
      integer, parameter :: nco2src=6
      real*8 co2_src(im,jm,nco2src)
      END MODULE CO2_SOURCES


      subroutine read_CO2_sources(nt,iact)
!@sum reads in CO2 sources and sinks
!@auth Jean Lerner
C****
C**** There are two monthly sources and 4 annual sources
C**** Annual sources are read in at start and re-start of run only
C**** Monthly sources are interpolated each day
      USE MODEL_COM, only: itime,jday,JDperY,im,jm
      USE CONSTANT, only: sday
      USE FILEMANAGER, only: openunit,closeunit, openunits,closeunits
      USE TRACER_COM, only: itime_tr0,trname
      use CO2_SOURCES, only: src=>co2_src,nsrc=>nco2src
      implicit none
      character*80 title
      logical :: ifirst=.true.
      real*8 adj(nsrc)
c     data adj/1.038,1.,1.,5.33,1.,1.749/          ! c
      data adj/3.81d0,3.67d0,3.67d0,19.54d0,3.67d0,6.42d0/     ! co2
!@var nanns,nmons: number of annual and monthly input files
      integer, parameter :: nanns=4,nmons=2
      integer ann_units(nanns),mon_units(nmons)
      character*12 :: ann_files(nanns) = 
     *  (/'CO2_FOS_FUEL','CO2_FERT    ','CO2_REGROWTH','CO2_LAND_USE'/)
      logical :: ann_bins(nanns)=(/.true.,.true.,.true.,.true./)
      character*9 :: mon_files(nmons) = (/'CO2_VEG  ','CO2_OCEAN'/)
      logical :: mon_bins(nmons)=(/.true.,.true./)
      real*8 tlca(im,jm,nmons),tlcb(im,jm,nmons)  ! for monthly sources
      integer i,j,nt,iact,iu,k,jdlast
      data jdlast /0/
      save ifirst,jdlast,tlca,tlcb

      if (itime.lt.itime_tr0(nt)) return
C****
C**** Annual Sources and sink
C**** Apply adjustment factors to bring sources into balance
C**** Annual sources are in KG C/M2/Y
C**** Sources need to be kg/m^2 s; convert /year to /s
C****
      if (ifirst) then
        call openunits(ann_files,ann_units,ann_bins,nanns)
        k = 0
        do iu = ann_units(1),ann_units(nanns)
          k = k+1
          call readt (iu,0,src(1,1,k),im*jm,src(1,1,k),1)
          src(:,:,k) = src(:,:,k)*adj(k)/(sday*JDperY)
        end do
        call closeunits(ann_units,nanns)
      endif
C****
C**** Monthly sources are interpolated to the current day
C****
C**** Also, Apply adjustment factors to bring sources into balance
C**** Monthly sources are in KG C/M2/S => src in kg/m^2 s
      if (iact.eq.1 .and. .not.ifirst) return
      ifirst = .false.
      call openunits(mon_files,mon_units,mon_bins,nmons)
      j = 0
      do k=nanns+1,nsrc
        j = j+1
        call read_monthly_sources(mon_units(j),jdlast,
     *    tlca(1,1,j),tlcb(1,1,j),src(1,1,k))
        src(:,:,k) = src(:,:,k)*adj(k)
      end do
      jdlast = jday
      call closeunits(mon_units,nmons)
      write(6,*) trname(nt),'Sources interpolated to current day'
C****
      return
      end subroutine read_CO2_sources


      SUBROUTINE set_tracer_source
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
      USE CO2_SOURCES, only: co2_src
      USE CH4_SOURCES, only: ch4_src
      implicit none
      integer :: i,j,ns,l,ky,n
      double precision :: source,sarea,steppy,base,steppd,x,airm,anngas,
     *  steph,stepx,stepp,tmon,by_dt

      by_dt = 1./DTsrc
C**** All sources are saved as kg/s
      do n=1,ntm

      if (itime.lt.itime_tr0(n)) cycle
      select case (trname(n))

      case default
!     write(6,*) ' Sources for ',trname(n),' are not in this routine'
C****
C**** Surface Sources of SF6 (Same grid as CFC)
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
            WRITE (6,'(A,I2,A,I9)')
     *        ' >> KY FOR SF6 IS',KY,' AT itime',itime
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

C****
C**** Sources and sinks for N2O:
C**** First layer is set to a constant 462.2 ppbm. (300 PPB V)
C****
      case ('N2O')
      do j=1,jm
        trsource(:,j,1,n) = (am(1,:,j)*dxyp(j)*462.2d-9
     *   -trm(:,j,1,n))*by_dt
      end do

C****
C**** Sources and sinks for 14C02
C**** NOTE: This tracer is supposed to start on 10/16
C**** The tracer is reset to specific values in layer 1
C****
      case ('14CO2')
      tmon = (itime-itime_tr0(n))*jmpery/(hrday*jdpery)  !(12./8760.)
      do j=1,jm/2
        trsource(:,j,1,n) = (am(1,:,j)*dxyp(j)*(4.82d-18*46./mair)*
     *   (44.5 + tmon*(1.02535d0 - tmon*(2.13565d-2 - tmon*8.61853d-5)))
     *   -trm(:,j,1,n))*by_dt
      end do
      do j=1+jm/2,jm
        trsource(:,j,1,n) = (am(1,:,j)*dxyp(j)*(4.82d-18*46./mair)*
     *   (73.0 - tmon*(0.27823d0 + tmon*(3.45648d-3 - tmon*4.21159d-5)))
     *   -trm(:,j,1,n))*by_dt
      end do

      end select

      end do
C****
      END SUBROUTINE set_tracer_source


      SUBROUTINE set_tracer_3Dsource
!@sum set_tracer_3Dsource calculates interactive sources for tracers
!@auth Jean Lerner
      USE TRACER_COM
      USE MODEL_COM, only: itime
      USE FLUXES, only: tr3Dsource
      implicit none
      integer n

C**** All sources are saved as kg/s
      do n=1,ntm

      if (itime.lt.itime_tr0(n)) cycle
      select case (trname(n))

      case default
C****
      case ('CH4')
      tr3Dsource(:,:,:,:,n) = 0.
      call Trop_chem_CH4(n,1)
      call Strat_chem_CH4(n,2)
C****
      case ('O3')
      tr3Dsource(:,:,:,:,n) = 0.
      call Trop_chem_O3(n,1)
      call Strat_chem_O3(n,2)
C****
      case ('N2O')
      tr3Dsource(:,:,:,:,n) = 0.
      call Strat_chem_N2O(n,1)
C****
      case ('CFC')
      tr3Dsource(:,:,:,:,n) = 0.
      call Strat_chem_CFC(n,1)
C****
      end select
      end do
      return
      END SUBROUTINE set_tracer_3Dsource


      SUBROUTINE get_14CO2_IC(CO2IJL)
!@sum GET_14CO2_IC Calculates initial distribution for 14CO2 tracer
!@auth J.Lerner (modified from program by G. Russell)
C**** NOTE: tracer is supposed to start on 10/16
C**** October 1963 14CO2 Concentrations for GCM  2/26/99
C**** 2/2/2: generalized code for modelE
C****           
      USE MODEL_COM, ONLY: im,jm,lm,ls1,sige,psf,ptop
      USE DYNAMICS, only: pedn
      USE FILEMANAGER, only: openunit,closeunit
      IMPLICIT NONE
      PARAMETER (kmw=60)
      REAL*4 CO2W(37,0:30) ! ,AMASS(IM,JM)
      real*8 p(0:60),CO2JK(JM,0:kmw),CO2IJL(IM,JM,LM)
      CHARACTER*80 TITLE
      integer i,j,jw,k,l,n,iu_in,iu_out,kmw
      real*8 pup,cup,pdn,cdn,psum,csum,psurf,ptrop,sig,w,zk !,stratm

C****
C**** Read in CO2 concentrations from workshop
C****
c     OPEN (1,FILE='workshop.14co2',STATUS='OLD')
      call openunit('14CO2_IC_DATA',iu_in,.false.)
      DO 110 N=1,10
  110 READ (iu_in,911)
      READ (iu_in,911) ((CO2W(J,K),J=1,13),K=0,30)
      READ (iu_in,911)
      READ (iu_in,911)
      READ (iu_in,912) ((CO2W(J,K),J=14,25),K=0,30)
      READ (iu_in,912)
      READ (iu_in,912)
      READ (iu_in,912) ((CO2W(J,K),J=26,37),K=0,30)
      call closeunit(iu_in)
C**** Calculate workshop pressure levels (Pa)
      DO K=0,kmw
        ZK = 2.*K
        P(K) = 100000.*10**(-ZK/16.)
      end do
C****
C**** Interpolate workshop CO2W to CO2JK on GCM latitudes
C****
      DO K=31,kmw       ! above 17 Pa set all equal
        CO2JK( 1,K) = 250.
        CO2JK(JM,K) = 250.
        DO J=2,JM-1
          CO2JK(J,K) = 250.
      end do; end do
      DO K=0,30
        CO2JK( 1,K) = CO2W( 1,K)
        CO2JK(JM,K) = CO2W(37,K)
        DO J=2,JM-1
          W = 1. + (J-1)*36./(JM-1)
          JW=W
          CO2JK(J,K) = CO2W(JW,K)*(JW+1-W) + CO2W(JW+1,K)*(W-JW)
      end do; end do
C****
C**** Read in GCM atmospheric mass for October
C****
c     OPEN (2,FILE='AIJX011.O',FORM='UNFORMATTED',STATUS='OLD')
   !  call openunit('OCTAirMass',iu_in)
   !  READ  (iu_in)
   !  READ  (iu_in) TITLE,AMASS  !! in kg/m**2
   !  call closeunit(iu_in)
   !  WRITE (6,*) 'Read: ',TITLE
C****
C**** Interpoate CO2J to CO2IJL on GCM grid boxes conserving vertical
C**** means
C****
C**** psf, ptrop, pdn ..... in pascals (mb*100)
      ptrop = ptop*100.
    ! STRATM = 101.9368   ! (kg/m**2)(10 mb x 100)/grav
      DO 440 J=1,JM
      DO 440 I=1,IM
    ! PDN = (AMASS(I,J)*SIGE(1)+STRATM)*GRAV
      PDN = pedn(1,i,j)*100.
      psurf = pdn
      CDN = CO2JK(J,0)
      K=1
      DO 430 L=1,LM
      PSUM = 0.
      CSUM = 0.
      if (l.eq.ls1) psurf = psf*100.
      PUP  =  SIGE(L+1)*(psurf-ptrop)+ptrop
  410 IF(P(K).LE.PUP)  GO TO 420
      PSUM = PSUM +  PDN-P(K)
      CSUM = CSUM + (PDN-P(K))*(CDN+CO2JK(J,K))/2.
      PDN  = P(K)
      CDN  = CO2JK(J,K)
      K=K+1
      if (k.gt.kmw) stop ' Please increase kmw in get_14CO2_IC'
      GO TO 410
C****
  420 CUP  = CO2JK(J,K) + (CO2JK(J,K-1)-CO2JK(J,K))*(PUP-P(K))/
     /       (P(K-1)-P(K))
      PSUM = PSUM +  PDN-PUP
      CSUM = CSUM + (PDN-PUP)*(CDN+CUP)/2.
      CO2IJL(I,J,L) = CSUM/PSUM
      PDN = PUP
  430 CDN = CUP
  440 CONTINUE
C****
C**** Scale data to proper units (10**-18 kg 14CO2/kg air)
C****
      CO2IJL(:,:,:) = CO2IJL(:,:,:)*4.82d0*(14.+16.*2.)/29.029d0
      RETURN
C****
  911 FORMAT (5X,13F5.0)
  912 FORMAT (5X,12F5.0)
      END SUBROUTINE get_14CO2_IC


      SUBROUTINE Trop_chem_CH4(n,ns)
!@sum Trop_chem_CH4 calculates tropospheric chemistry for CH4
!@+     by applying a pre-determined chemical loss rate
!@auth Jean Lerner
      USE TRACER_COM
      USE DYNAMICS, only: ltropo
      USE MODEL_COM, only: im,jm,lm,ls1,byim,jyear,jhour,jday,itime
      USE GEOM, only: imaxj
      USE FLUXES, only: tr3Dsource
      USE FILEMANAGER, only: openunit,closeunit
      implicit none
      integer n,ns,i,j,l,ifile
      REAL*4 frqlos(im,jm,ls1),taux
      real*8 tauy,tune
      parameter (tune = 445./501.)
      logical, save :: ifirst=.true.
      save frqlos,tauy

C**** Read chemical loss rate dataset (5-day frequency)
      if (mod(jday,5).gt.0 .and. .not.ifirst) go to 550
c     if (tofday.ne.0. .and. .not.ifirst) go to 550
      if (jhour.ne.0 .and. .not.ifirst) go to 550
      ifirst = .false.
      call openunit('CH4_TROP_FRQ',ifile,.true.,.true.)
  510 continue
      read (ifile,end=515) taux,frqlos
      tauy = nint(taux)+(jyear-1950)*8760.
      IF (itime+ 60..gt.tauy+120.) go to 510
      IF (itime+180..le.tauy+120.) then
        write(6,*)' PROBLEM MATCHING itime ON FLUX FILE',TAUX,TAUY,JYEAR
        STOP 'PROBLEM MATCHING itime ON FLUX FILE'
      end if
      rewind (ifile)
      go to 518
C**** FOR END OF YEAR, USE FIRST RECORD
  515 continue
      rewind (ifile)
      read (ifile,end=515) taux,frqlos
      rewind (ifile)
      tauy = nint(taux)+(jyear-1950)*8760.
  518 continue
      call closeunit(ifile)
      WRITE(6,'(A,2F10.0,2I10)')
     * ' *** CHEMICAL LOSS RATES READ FOR TAUX,TAUY,itime,JYEAR=',
     * TAUX,TAUY,itime,JYEAR
C**** AVERAGE POLES
      do l=1,ls1
        frqlos(1, 1,l) = sum(frqlos(:, 1,l))*byim
        frqlos(1,jm,l) = sum(frqlos(:,jm,l))*byim
      end do
C**** APPLY AN AD-HOC FACTOR TO BRING INTO BALANCE
      frqlos(:,:,:) = frqlos(:,:,:)*tune
C**** Apply the chemistry
  550 continue
      do j=1,jm
      do i=1,imaxj(j)
        do l = 1,ltropo(i,j)
          if (l.le.ls1) tr3Dsource(i,j,l,ns,n) = frqlos(i,j,l)
        end do
      end do
      end do
      IF (jhour.eq.0.) write(6,'(a,i10,F10.0)') 
     *  ' Tropo chem for CH4 performed',itime,TAUY
      return
      END SUBROUTINE Trop_chem_CH4


      SUBROUTINE Strat_chem_CH4(n,ns)
!@sum Trop_chem_CH4 calculates stratospheric chemistry for CH4
!@auth Jean Lerner
      USE TRACER_COM
      USE FLUXES, only: tr3Dsource
      implicit none
      integer n,ns
      return
      END SUBROUTINE Strat_chem_CH4


      SUBROUTINE Trop_chem_O3(n,ns)
!@sum Trop_chem_O3 calculates tropospheric chemistry for O3
!@auth Jean Lerner
      USE TRACER_COM
      USE FLUXES, only: tr3Dsource
      implicit none
      integer n,ns
      return
      END SUBROUTINE Trop_chem_O3


      SUBROUTINE Strat_chem_O3(n,ns)
!@sum Trop_chem_O3 calculates stratospheric chemistry for O3
!@auth Jean Lerner
      USE TRACER_COM
      USE FLUXES, only: tr3Dsource
      implicit none
      integer n,ns
      return
      END SUBROUTINE Strat_chem_O3


      SUBROUTINE Strat_chem_N2O(n,ns)
!@sum Trop_chem_N2O calculates stratospheric chemistry for N2O
!@auth Jean Lerner
      USE TRACER_COM
      USE FLUXES, only: tr3Dsource
      implicit none
      integer n,ns
      return
      END SUBROUTINE Strat_chem_N2O


      SUBROUTINE Strat_chem_CFC(n,ns)
!@sum Trop_chem_CFC calculates stratospheric chemistry for CFC
!@auth Jean Lerner
      USE TRACER_COM
      USE FLUXES, only: tr3Dsource
      implicit none
      integer n,ns
      return
      END SUBROUTINE Strat_chem_CFC
#endif


#ifdef TRACERS_WATER
C---SUBROUTINES FOR TRACER WET DEPOSITION-------------------------------

      SUBROUTINE GET_COND_FACTOR(L,N,WMXTR,FCLOUD,FQ0,fq)
!@sum  GET_COND_FACTOR calculation of condensate fraction for tracers
!@+    within or below convective or large-scale clouds. Gas 
!@+    condensation uses Henry's Law if not freezing.
!@auth Dorothy Koch (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on CB436TdsM23 CLOUDCHCC and CLOUDCHEM subroutines)
c
C**** GLOBAL parameters and variables:
      USE CLOUDS, only: PL, TL, NTIX
      USE TRACER_COM, only: tr_RKD,tr_DHD,nWATER,nGAS,nPART,tr_wd_TYPE
      USE CONSTANT, only: TF, BYGASC, MAIR
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
!@var WMXTR mixing ratio of water available for tracer condensation ( )?
!@var RKD dummy variable (= tr_RKD*EXP[ ])
      REAL*8, PARAMETER :: BY298K=3.3557D-3
      REAL*8 Ppas, tfac, ssfac, RKD
      REAL*8,  INTENT(IN) :: fq0, FCLOUD, WMXTR
      REAL*8,  INTENT(OUT):: fq
      INTEGER, INTENT(IN) :: L, N
c
C**** CALCULATE the fraction of tracer mass that becomes condensate:
c
      SELECT CASE(tr_wd_TYPE(NTIX(N)))                 
        CASE(nGAS)                                  ! gas tracer
          fq = 0.D0                                   ! frozen case
          IF(TL(L).ge.TF) THEN                        ! if not frozen then: 
            Ppas = PL(L)*1.D2                         ! pressure to pascals
            tfac = (1.D0/TL(L) - BY298K)*BYGASC
            IF(tr_DHD(NTIX(N)).ne.0.D0) THEN
              RKD=tr_RKD(NTIX(N))*DEXP(-tr_DHD(NTIX(N))*tfac)
            ELSE  
              RKD=tr_RKD(NTIX(N))
            END IF
c           clwc=WMXTR*MAIR*1.D-3*Ppas*BYGASC/(TL(L)*FCLOUD)
c           ssfac=RKD*GASC*TL(L)*clwc   ! Henry's Law
            ssfac=RKD*WMXTR*MAIR*1.D-3*Ppas/FCLOUD
            fq=ssfac / (1.D0 + ssfac)
          END IF
        CASE(nWATER)                                ! water tracer
          fq = FQ0                                  
        CASE(nPART)                                 ! particulate tracer
          fq = 0.D0                                   ! temporarily zero.
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
!@sum  GET_PREC_FACTOR calculation of the precipitation scavenging fraction
!@+    for tracers WITHIN large scale clouds. Current version uses the
!@+    first order removal rate based on [Giorgi and Chameides, 1986], for
!@+    gaseous and particulate tracers. 
!@auth Dorothy Koch (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on CB436TdsM23 RAINOUT subroutine)
c
C**** GLOBAL parameters and variables:
      USE CLOUDS, only: NTIX
      USE MODEL_COM, only: dtsrc
      USE TRACER_COM, only: nWATER, nGAS, nPART, tr_wd_TYPE
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
        CASE(nWATER)                              ! water/original method
          fq = FQ0                               
        CASE(nPART)                               ! aerosols
          IF(BELOW_CLOUD) THEN
            fq = 0.D0
          ELSE
c           minus preserves FPRT sign convention in LSCOND
            fq = -(CM*DTsrc*(EXP(-CM*DTsrc)- 1.0)*FCLD)
          END IF
        CASE DEFAULT                              ! error
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
!@var b_beta_DT precipitating grid box fraction from lowest percipitating
!@+   layer. The name was chosen to correspond to Koch et al. p. 23,802.
!@var N index for tracer number loop
      INTEGER, INTENT(IN) :: N
      REAL*8, INTENT(OUT):: FQ
      REAL*8, INTENT(IN) :: PREC, b_beta_DT
      REAL*8, PARAMETER :: rc_wash = 1.D-1
C
      SELECT CASE(tr_wd_TYPE(NTIX(N)))                 
        CASE(nGAS)                                ! gas
          fq = 0.D0
        CASE(nWATER)                              ! water/original method
          fq = 0.D0                              
        CASE(nPART)                               ! aerosols
          fq = b_beta_DT*(DEXP(-PREC*rc_wash)-1.)
        CASE DEFAULT                              ! error
          STOP 'tr_wd_TYPE(NTIX(N)) out of range in WASHOUT_TRACER'                              
      END SELECT   
c      
      RETURN
      END SUBROUTINE GET_WASH_FACTOR


      SUBROUTINE GET_EVAP_FACTOR(N,FQ0,fq)
!@sum  GET_EVAP_FACTOR calculation of the evaporation fraction for tracers.
!@auth Dorothy Koch (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on CB436TdsM23 EVAPD and WASH_EVAP routines)
c
C**** GLOBAL parameters and variables:
      USE CLOUDS, only: NTIX
      USE TRACER_COM, only: tr_evap_fact, tr_wd_TYPE
c
      IMPLICIT NONE
c      
C**** Local parameters and variables and arguments:
!@var FQ            fraction of tracer evaporated
!@var FQ0 [default] fraction of tracer evaporated
!@var N index for tracer number loop
      INTEGER, INTENT(IN) :: N
      REAL*8,  INTENT(OUT):: FQ
      REAL*8,  INTENT(IN) :: FQ0
c
      fq=FQ0*tr_evap_fact(tr_wd_TYPE(NTIX(N)))
c
      RETURN
      END SUBROUTINE GET_EVAP_FACTOR 
#endif
