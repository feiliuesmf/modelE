!@sum  TRACERS: tracer routines that are tracer-dependent
!@+    Routines included: 
!@+      Diagnostic specs: TRACER_DIAG_COM, INIT_TRACER, SET_TCON
!@+      Tracer definitions: TRACER_IC, TRACER_SOURCE
!@auth Jean Lerner/Gavin Schmidt
!@ver  1.0

      MODULE TRACER_DIAG_COM
!@sum Tracer diagnostic arrays
!@auth Jean Lerner
!ver   1.0
      USE TRACER_COM, only: ntm 
      USE MODEL_COM, only: im,jm,lm
      USE DAGCOM, only: npts !npts are conservation quantities
      IMPLICIT NONE
      SAVE

C**** TAIJS  <<<< KTAIJS and IJTS_xx are Tracer-Dependent >>>>
C**** TAJLS  <<<< KTAJLS and JLS_xx are Tracer-Dependent >>>>

!@var MMR_to_VMR: converts tracer mass mixing ratio to volume mr 
      REAL*8, DIMENSION(NTM) :: MMR_to_VMR

C**** TAIJLN
!@parm KTAIJLN number of 3D diagnostics for each tracer
      integer, parameter :: ktaijln=ntm*lm
!@var TAIJLN 3D tracer diagnostics (all tracers)
      real*8, dimension(im,jm,lm,ntm) :: taijln
!@var SNAME_IJT, UNITS_IJT: Names and units of lat-sigma tracer IJ diags
      character(len=30), dimension(lm,ntm) :: sname_ijt,units_ijt
!@var LNAME_IJT: descriptions of tracer IJ diagnostics
      character(len=80), dimension(lm,ntm) :: lname_ijt
!@var SCALE_IJT: printout scaling factor for tracer IJ diagnostics
      double precision, dimension(lm,ntm) :: scale_ijt
!@var IR_IJT: range index of IJ diagnosts
      integer, dimension(ntm) :: ir_ijt
!@var IA_IJT: idacc-number for tracer IJ diags (same for all tracers)
      integer ia_ijt
!@var IJTC_POWER: power of 10 used for tracer IJ concentration diags
      integer, dimension(ntm) :: ijtc_power
!@var IJTM_POWER: power of 10 used for tracer IJ mass unit diags
      integer, dimension(ntm) :: ijtm_power
!@var IJT_XX Names for TAIJLN diagnostics (Not currently used)

C**** TAIJKN
!@parm KTAIJK number of 3D diagnostics for each tracer
!@parm KTAIJKN total number of 3D diagnostics for tracers (ntm*ktaijk)
      integer, parameter :: ktaijk=3,ktaijkn=ntm*ktaijk
!@var IJT_XX names for taijkn diagnostics 
      integer ijkn_conc,ijkn_surf,ijkn_mass
!@var TAIJKN lat/lon tracer diagnostics (all tracers)
      real*8, dimension(im,jm,ktaijk,ntm) :: taijkn
!@var SCALE_IJKN: printout scaling factor for tracer IJK diagnostics
      double precision, dimension(ktaijk,ntm) :: scale_ijkn
!@var SNAME_IJKN,UNITS_IJKN: Names and units of lat-sigma tracer diags
      character(len=30), dimension(ktaijk,ntm) :: sname_ijkn,units_ijkn
!@var LNAME_IJKN: descriptions of tracer IJK diags
      character(len=80), dimension(ktaijk,ntm) :: lname_ijkn

C**** TAIJS  <<<< KTAIJS and IJTS_xx are Tracer-Dependent >>>>
!@parm KTAIJS number of special lat/lon tracer diagnostics
      integer, parameter :: ktaijs=8
!@var TAIJS  lat/lon special tracer diagnostics; sources, sinks, etc.
      REAL*8, DIMENSION(IM,JM,ktaijs) :: TAIJS
!@var ijts_xx Names for TAIJS diagnostics
      INTEGER ijts_SF6_source,ijts_Rn222_source,
     *  ijts_CO2_fossil_fuel,ijts_CO2_fertilization,ijts_CO2_regrowth,
     *  ijts_CO2_land_use,ijts_CO2_Ecosystem,ijts_CO2_Ocean
!@var SNAME_IJTS, UNITS_IJTS: Names & units of lat-sigma tracer diags
      character(len=30), dimension(ktaijs) :: sname_ijts,units_ijts
!@var LNAME_IJTS: descriptions of tracer IJTS diags
      character(len=80), dimension(ktaijs) :: lname_ijts
!@var SCALE_IJTS: printout scaling factor for tracer IJTS diagnostics
      double precision, dimension(ktaijs) :: scale_ijts
!@var IA_IJTS: idacc-number for tracer source/sink IJ diags 
      integer ia_ijts(ktaijs)
!@var ijts_power: power of 10 used for tracer IJ source/sink diags
      INTEGER, DIMENSION(ktaijs) :: ijts_power
!@var ijts_source: tracer index associated with a TAIJS diagnostic
      INTEGER, DIMENSION(ktaijs) :: ijts_source

C**** TAJLN
!@parm ktajl,ktajlx number of TAJL tracer diagnostics; 
!@+          ktajlx includes composits
      INTEGER, PARAMETER :: ktajl=8, ktajlx=ktajl+2
!@var TAJLN  vertical tracer diagnostics (all tracers)
      REAL*8, DIMENSION(JM,LM,ktajlx,NTM) :: TAJLN
!@var jlnt_xx Names for TAJLN diagnostics
      INTEGER jlnt_conc,jlnt_nt_tot,jlnt_nt_mm,jlnt_vt_tot,
     &  jlnt_vt_mm,jlnt_mc,jlnt_turb,jlnt_lscond
!@var SNAME_JLN: Names of lat-sigma tracer JL diagnostics
      character(len=30), dimension(ktajlx,ntm) :: sname_jln
!@var LNAME_JLN,UNITS_JLN: descriptions/units of tracer JL diagnostics
      character(len=50), dimension(ktajlx,ntm) :: lname_jln,units_jln
!@var SCALE_JLQ: printout scaling factor for tracer JL diagnostics
      double precision, dimension(ktajlx) :: scale_jlq
!@var IA_JLQ,JGRID_JLQ: idacc-numbers,gridtypes for tracer JL diags
      integer, dimension(ktajlx) :: ia_jlq,jgrid_jlq
!@var JLQ_POWER: power of 10 used for tracer JL diagnostics
!@+        It is associated with a specific physical process
      integer jlq_power(ktajlx)
      integer ktajlp !@var helps keeps track of
!@var scale_jln: Coefficient to print MMR or VMR
      REAL*8, DIMENSION(ntm) :: scale_jln

C**** TAJLS  <<<< KTAJLS and JLS_xx are Tracer-Dependent >>>>
!@parm ktajls number of source/sink TAJLS tracer diagnostics;
      INTEGER, PARAMETER :: ktajls=9
!@var TAJLS  JL special tracer diagnostics for sources, sinks, etc
      REAL*8, DIMENSION(JM,LM,ktajls) :: TAJLS
!@var jls_xx Names for TAJLS diagnostics
      INTEGER jls_SF6_source,jls_RN222_source,jls_RN222_sink,
     *  jls_CO2_fossil_fuel,jls_CO2_fertilization,jls_CO2_regrowth,
     *  jls_CO2_land_use,jls_CO2_Ecosystem,jls_CO2_Ocean
!@var jls_source: Number of source/sink tracer JL diagnostics/tracer 
      integer, dimension(ktajls) :: jls_source
!@var SNAME_JLS: Names of lat-sigma tracer JL sources/sinks
      character(len=30), dimension(ktajls) :: sname_jls
!@var LNAME_JLS,UNITS_JLS: descriptions/units of tracer JLS diags
      character(len=50), dimension(ktajls) :: lname_jls,units_jls
!@var SCALE_JLS: printout scaling factor for tracer JLS diagnostics
      double precision, dimension(ktajls) :: scale_jls
!@var IA_JLS,JGRID_JLS: idacc-numbers,gridtypes for tracer JL diags
      integer, dimension(ktajls) :: ia_jls,jgrid_jls
!@var JLS_POWER: power of 10 used for tracer JLS diagnostics
      integer jls_power(ktajls)
!@var jls_ltop: Top layer for this diagnostic
      integer jls_ltop(ktajls)

C**** TCONSRV
!@parm NTCONS Maximum Number of special tracer conservation points
      INTEGER, PARAMETER :: ntcons=20
!@parm KTCON total number of conservation diagnostics for tracers
      INTEGER, PARAMETER :: KTCON=npts+ntcons+2
!@var TCONSRV conservation diagnostics for tracers
      DOUBLE PRECISION, DIMENSION(JM,ktcon,ntm) :: TCONSRV
!@var SCALE_TCON scales for tracer conservation diagnostics
      DOUBLE PRECISION, DIMENSION(ktcon,ntm) :: SCALE_TCON
!@var TITLE_TCON titles for tracer conservation diagnostics
      CHARACTER*38, DIMENSION(ktcon,ntm) :: TITLE_TCON
!@var IA_TCON IDACC numbers for tracer conservation diagnostics
      INTEGER, DIMENSION(ktcon) ::  IA_TCON
!@var NSUM_TCON indices for summation of conservation diagnostics
      INTEGER, DIMENSION(ktcon,ntm) :: NSUM_TCON
!@var NOFMT indices for TCONSRV array
      INTEGER, DIMENSION(ktcon,ntm) :: NOFMT
!@var CONPTS names of special processes for tracer conservation diags
      CHARACTER*16, DIMENSION(ntcons) :: CONPTS
!@var kt_power_inst,kt_power_change: Exponents for tracer conservation
      INTEGER, DIMENSION(ntm) :: kt_power_inst,kt_power_change
!@var name_tconsrv,lname_tconsrv,units_tconsrv: for tracer conservation
      character(len=20), dimension(ktcon,ntm) :: 
     &   name_tconsrv,units_tconsrv
      character(len=80), dimension(ktcon,ntm) :: lname_tconsrv
!@var SCALE_INST,SCALE_CHANGE: Scale factors for tracer conservation
      double precision, dimension(ntm) :: SCALE_INST,SCALE_CHANGE

C----------------------------------------------------
!@param KTACC total number of tracer diagnostic words
!@var TACC: Contains all tracer diagnostic accumulations
      INTEGER, PARAMETER :: 
     *  ktacc=IM*JM*LM*NTM + IM*JM*ktaijk*NTM + IM*JM*ktaijs + 
     *        JM*LM*ktajlx*NTM + JM*LM*ktajls + JM*NTM*ktcon
      COMMON /TACCUM/ TAIJLN,TAIJKN,TAIJS,TAJLN,TAJLS,TCONSRV
      DOUBLE PRECISION, DIMENSION(KTACC) :: TACC
      EQUIVALENCE (TACC,TAIJLN)
C----------------------------------------------------
      END MODULE TRACER_DIAG_COM


      SUBROUTINE init_tracer
!@sum init_tracer initializes trace gas attributes and diagnostics
!@auth J. Lerner
!@calls sync_param, SET_TCON
      use DAGCOM, only: ia_src
      USE MODEL_COM, only: dtsrc
      use BDIJ, only: ir_log2
      USE TRACER_COM
      USE TRACER_DIAG_COM
      USE CONSTANT, only: mair
      USE PARAM
      implicit none
      integer :: l,k,n
      character*20 sum_unit(ntm),inst_unit(ntm)   ! for conservation
      character*10 CMR
      logical :: qcon(KTCON-1), T=.TRUE. , F=.FALSE.
!@var to_volume_MixRat: For printout of tracer concentration
!@+   to_volume_MixRat=1: printout is in Volume Mixing Ratio
!@+   to_volume_MixRat=0: printout is in Mass Mixing Ratio
      INTEGER :: to_volume_MixRat=0 
      character*50 :: unit_string

C**** Get itime_tr0 from rundeck if it exists
      call sync_param("itime_tr0",itime_tr0,ntm)
C**** Get to_volume_MixRat from rundecks if it exists
      call sync_param("to_volume_MixRat",to_volume_MixRat)

C**** Get factor to convert from mass mixing ratio to volume mr
      if (to_volume_MixRat .eq.1) then
        MMR_to_VMR(:) = mair/tr_mm(:)
        cmr = ' V/V air'
      else
        MMR_to_VMR(:) = 1.d0
        cmr = ' kg/kg air'
      endif
C****
C**** TAJLN(J,L,KQ,N)  (SUM OVER LONGITUDE AND TIME OF)
C****
C**** jlq_power Exponent associated with a physical process 
C****      (for printing tracers).   (ntm_power+jlq_power)
C**** jls_power Exponent associated with a source/sink (for printing)
C**** Defaults for JLN
      scale_jlq(:) = 1./DTsrc
      jgrid_jlq(:) = 1
      ia_jlq(:) = ia_src

C**** Tracer concentration
      do n=1,ntm
        k = 1        ! <<<<< Be sure to do this
        jlnt_conc = k
        sname_jln(k,n) = trim(trname(n))//'_CONCENTRATION' 
        lname_jln(k,n) = trim(trname(n))//' CONCENTRATION' 
        jlq_power(k) = 0.
        units_jln(k,n) = unit_string(ntm_power(n)+jlq_power(k),cmr)
        scale_jlq(k) = 1.d0
        scale_jln(n) = MMR_to_VMR(n)

C**** Physical processes affecting tracers
C****   F (TOTAL NORTHWARD TRANSPORT OF TRACER MASS)  (kg)
        k = k + 1 
        jlnt_nt_tot = k
        sname_jln(k,n) = 'tr_nt_tot_'//trname(n)
        lname_jln(k,n) = 'TOTAL NORTHWARD TRANSPORT OF '//
     &     trim(trname(n))//' MASS'
        jlq_power(k) = 11.
        jgrid_jlq(k) = 2
C****   STM/SM (MEAN MERIDIONAL N.T. OF TRACER MASS)  (kg)
        k = k + 1 
        jlnt_nt_mm = k
        sname_jln(k,n) = 'tr_nt_mm_'//trname(n)
        lname_jln(k,n) = 'NORTHWARD TRANS. OF '//
     &     trim(trname(n))//' MASS BY MERIDIONAL CIRC.'
        jlq_power(k) = 10.
        jgrid_jlq(k) = 2
C****   F (TOTAL VERTICAL TRANSPORT OF TRACER MASS)  (kg)
        k = k + 1
        jlnt_vt_tot = k
        sname_jln(k,n) = 'tr_vt_tot_'//trname(n)
        lname_jln(k,n) = 'TOTAL VERTICAL TRANSPORT OF '//
     &     trim(trname(n))//' MASS'
        jlq_power(k) = 11.
C****   STM/SM (MEAN MERIDIONAL V.T. OF TRACER MASS)  (kg)
        k = k + 1
        jlnt_vt_mm = k
        sname_jln(k,n) = 'tr_vt_mm_'//trname(n)
        lname_jln(k,n) = 'VERTICAL TRANS. OF '//
     &     trim(trname(n))//' MASS BY MERIDIONAL CIRC.'
        jlq_power(k) = 10.
C****   TMBAR-TM (CHANGE OF TRACER MASS BY MOIST CONVEC)(kg)
        k = k + 1
        jlnt_mc = k
        sname_jln(k,n) = 'tr_mc_'//trname(n)
        lname_jln(k,n) = 'CHANGE OF '//
     &     trim(trname(n))//' MASS BY MOIST CONVECTION'
        jlq_power(k) = 10.
C****   TMBAR-TM (CHANGE OF TRACER MASS BY Large-scale CONDENSE)  (kg)
        k = k + 1
        jlnt_lscond = k
        sname_jln(k,n) = 'tr_lscond'//trname(n)
        lname_jln(k,n) ='CHANGE OF '//
     &     trim(trname(n))//' MASS BY LARGE-SCALE CONDENSE'
        jlq_power(k) = 10.
C****   TMBAR-TM (CHANGE OF TRACER MASS BY DRY CONVEC)  (kg)
        k = k + 1
        jlnt_turb = k
        sname_jln(k,n) = 'tr_turb_'//trname(n)
        lname_jln(k,n) = 'CHANGE OF '//
     &     trim(trname(n))//' MASS BY TURBULENCE/DRY CONVECTION'
        jlq_power(k) = 10.
      end do

      if (k.gt. ktajl) then
        write (6,*) 
     &   'tjl_defs: Increase ktajl=',ktajl,' to at least ',k
        stop 'ktajl too small'
      end if

C**** Construct UNITS string for output
      do n=1,ntm
      do k=2,ktajl
      units_jln(k,n) = unit_string(ntm_power(n)+jlq_power(k),' kg/s')
      end do; end do

C****
C**** Tracer sources and sinks
C**** Defaults for jls (sources, sinks, etc.)
C**** This needs to be 'hand coded' depending on circumstances
      do k=1,ktajls  ! max number of sources and sinks
        jgrid_jls(k) = 1
        ia_jls(k) = ia_src
        scale_jls(k) = 1./(DTsrc+1.D-20)
      end do
      jls_source(:) = 0
        
      k = 0
      n = n_air  !no special sources

      n = n_SF6
        k = k + 1 
        jls_SF6_source = k
        sname_jls(k) = 'Layer_1_source_of_'//trname(n)
        lname_jls(k) = 'SF6 CFC-GRID SOURCE, LAYER 1'
        jls_source(k) = n
        jls_ltop(k) = 1
        jls_power(k) = -3.
        units_jls(k) = unit_string(jls_power(k),' kg/s')
      n = n_Rn222
        k = k + 1 
        jls_RN222_sink = k
        sname_jls(k) = 'Decay_of_'//trname(n)
        lname_jls(k) = 'LOSS OF RADON-222 BY DECAY'
        jls_source(k) = n
        jls_ltop(k) = lm
        jls_power(k) = -11.
        units_jls(k) = unit_string(jls_power(k),' kg/s')
        k = k + 1 
        jls_RN222_source = k
        sname_jls(k) = 'Ground_Source_of_'//trname(n)
        lname_jls(k) = 'RADON-222 SOURCE, LAYER 1'
        jls_source(k) = n
        jls_ltop(k) = 1
        jls_power(k) = -10.
        units_jls(k) = unit_string(jls_power(k),' kg/s')
! keep AIJ and AJL CO2 sources in same order !!
      n = n_CO2
        k = k + 1 
        jls_CO2_fossil_fuel = k
        sname_jls(k) = 'Fossil_fuel_source_'//trname(n)
        lname_jls(k) = 'CO2 Fossil fuel source (Marland)'
        jls_source(k) = n
        jls_ltop(k) = 1
        jls_power(k) = 3
        units_jls(k) = unit_string(jls_power(k),' kg/s')
        k = k + 1 
        jls_CO2_fertilization = k
        sname_jls(k) = 'fertilization_sink_'//trname(n)
        lname_jls(k) = 'CO2 fertilization sink (Friedlingstein)'
        jls_source(k) = n
        jls_ltop(k) = 1
        jls_power(k) = 3
        units_jls(k) = unit_string(jls_power(k),' kg/s')
        k = k + 1 
        jls_CO2_regrowth = k
        sname_jls(k) = 'Northern_forest_regrowth_'//trname(n)
        lname_jls(k) = 'CO2 Northern forest regrowth sink'
        jls_source(k) = n
        jls_ltop(k) = 1
        jls_power(k) = 3
        units_jls(k) = unit_string(jls_power(k),' kg/s')
        k = k + 1 
        jls_CO2_land_use = k
        sname_jls(k) = 'Land_Use_Modification_'//trname(n)
        lname_jls(k) = 'CO2 from Land use modification (Houton)'
        jls_source(k) = n
        jls_ltop(k) = 1
        jls_power(k) = 3
        units_jls(k) = unit_string(jls_power(k),' kg/s')
        k = k + 1 
        jls_CO2_Ecosystem = k
        sname_jls(k) = 'Ecosystem_exchange_'//trname(n)
        lname_jls(k) = 'CO2 Ecosystem exchange (Matthews)'
        jls_source(k) = n
        jls_ltop(k) = 1
        jls_power(k) = 3
        units_jls(k) = unit_string(jls_power(k),' kg/s')
        k = k + 1 
        jls_CO2_Ocean = k
        sname_jls(k) = 'Ocean_exchange_'//trname(n)
        lname_jls(k) = 'CO2 Ocean exchange'
        jls_source(k) = n
        jls_ltop(k) = 1
        jls_power(k) = 3
        units_jls(k) = unit_string(jls_power(k),' kg/s')

      if (k.gt. ktajls) then
        write (6,*) 
     &   'tjl_defs: Increase ktajls=',ktajls,' to at least ',k
        stop 'ktajls too small'
      end if

C**** CONTENTS OF TAIJLN(I,J,LM,N)  (SUM OVER TIME OF)
C****        TML (M*M * KG TRACER/KG AIR)
C**** Set defaults that are true for all tracers and layers
      ia_ijt    = ia_src
      ir_ijt(:) = ir_log2   !n
      ijtc_power(:) = ntm_power(:)+1   !n for concentration
      ijtm_power(:) = ntm_power(:)+4   !n for integrated mass
C**** Tracer concentrations (AIJLN)
      do n=1,ntm
      do l=1,lm
        write(sname_ijt(l,n),'(a,i2.2)') trim(TRNAME(n))//'_L_',l 
        write(lname_ijt(l,n),'(a,i2)')   trim(TRNAME(n))//' L ',l
        units_ijt(l,n) = unit_string(ijtc_power(n),cmr)
        scale_ijt(l,n) = MMR_to_VMR(n)*10.**(-ijtc_power(n))
      end do
      end do

C**** AIJKN
C****     1  TM (SUM OVER ALL LAYERS) (M*M * KG TRACER/KG AIR)
C****     2  TM1-TZM1 (SURFACE TRACER CONC.) (M*M * KG TRACER/KG AIR)
C****     3  TM (SUM OVER ALL LAYERS) (M*M * KG TRACER)
      do n=1,ntm
C**** Summation of mass over all layers
      k = 1        ! <<<<< Be sure to do this
      ijkn_mass = k
        write(sname_ijkn(k,n),'(a,i2)') trim(TRNAME(n))//'_Total_Mass'
        write(lname_ijkn(k,n),'(a,i2)') trim(TRNAME(n))//' Total Mass'
        units_ijkn(k,n) = unit_string(ijtm_power(n),' kg/m^2')
        scale_ijkn(k,n) = 10.**(-ijtm_power(n))
C**** Average concentration over layers
      k = k+1
      ijkn_conc = k
        write(sname_ijkn(k,n),'(a,i2)') trim(TRNAME(n))//'_Average'
        write(lname_ijkn(k,n),'(a,i2)') trim(TRNAME(n))//' Average'
        units_ijkn(k,n) = unit_string(ijtc_power(n),cmr)
        scale_ijkn(k,n) = MMR_to_VMR(n)*10.**(-ijtc_power(n))
C**** Surface concentration
      k = k+1
      ijkn_surf = k
        write(sname_ijkn(k,n),'(a,i2)') trim(TRNAME(n))//'_At_Surface'
        write(lname_ijkn(k,n),'(a,i2)') trim(TRNAME(n))//' At Surface'
        units_ijkn(k,n) = unit_string(ijtc_power(n),cmr)
        scale_ijkn(k,n) = MMR_to_VMR(n)*10.**(-ijtc_power(n))
      end do

      if (k .gt. ktaijk) then
        write (6,*) 
     &   'tjl_defs: Increase ktaijk=',ktaijk,' to at least ',k
        stop 'ktaijk too small'
      end if

C**** Tracer sources and sinks
C**** Defaults for ijts (sources, sinks, etc.)
C**** This needs to be 'hand coded' depending on circumstances
      k = 1
        ijts_SF6_source = k
        ijts_source(k) = n_SF6
        ia_ijts(k) = ia_src   
        lname_ijts(k) = 'SF6 Layer 1 SOURCE'
        sname_ijts(k) = 'SF6_CFC-GRID_SOURCE,_LAYER_1'
        ijts_power(k) = -15
        units_ijts(k) = unit_string(ijts_power(k),' kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k+1
        ijts_Rn222_source = k
        ijts_source(k) = n_Rn222
        ia_ijts(k) = ia_src   
        lname_ijts(k) = 'Rn222 L 1 SOURCE'
        sname_ijts(k) = 'Radon-222_SOURCE,_Layer_1'
        ijts_power(k) = -21.
        units_ijts(k) = unit_string(ijts_power(k),' kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
! keep AIJ and AJL CO2 sources in same order !!
      n = n_CO2
      k = k + 1 
        ijts_CO2_fossil_fuel = k
        ijts_source(k) = n
        ia_ijts(k) = ia_src   
        sname_ijts(k) = 'Fossil_fuel_source_'//trname(n)
        lname_ijts(k) = 'CO2 Fossil fuel src'
        ijts_power(k) = -11.
        units_ijts(k) = unit_string(ijts_power(k),' kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1 
        ijts_CO2_fertilization = k
        ijts_source(k) = n
        ia_ijts(k) = ia_src   
        sname_ijts(k) = 'fertilization_sink_'//trname(n)
        lname_ijts(k) = 'CO2 fertilization'
        ijts_power(k) = -11.
        units_ijts(k) = unit_string(ijts_power(k),' kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1 
        ijts_CO2_regrowth = k
        ijts_source(k) = n
        ia_ijts(k) = ia_src   
        sname_ijts(k) = 'Northern_forest_regrowth_'//trname(n)
        lname_ijts(k) = 'CO2 North forest regrowth'
        ijts_power(k) = -11.
        units_ijts(k) = unit_string(ijts_power(k),' kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1 
        ijts_CO2_land_use = k
        ijts_source(k) = n
        ia_ijts(k) = ia_src   
        sname_ijts(k) = 'Land_Use_Modification_'//trname(n)
        lname_ijts(k) = 'CO2 from Land use mods'
        ijts_power(k) = -11.
        units_ijts(k) = unit_string(ijts_power(k),' kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1 
        ijts_CO2_Ecosystem = k
        ijts_source(k) = n
        ia_ijts(k) = ia_src   
        sname_ijts(k) = 'Ecosystem_exchange_'//trname(n)
        lname_ijts(k) = 'CO2 Ecosystem exch'
        ijts_power(k) = -11.
        units_ijts(k) = unit_string(ijts_power(k),' kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k + 1 
        ijts_CO2_Ocean = k
        ijts_source(k) = n
        ia_ijts(k) = ia_src   
        sname_ijts(k) = 'Ocean_exchange_'//trname(n)
        lname_ijts(k) = 'CO2 Ocean exchange'
        ijts_power(k) = -11.
        units_ijts(k) = unit_string(ijts_power(k),' kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc

      if (k .gt. ktaijs) then
        write (6,*)'ijt_defs: Increase ktaijs=',ktaijs,' to at least ',k
        stop 'ktaijs too small'
      end if

C**** Initialize conservation diagnostics
C**** To add a new conservation diagnostic:
C****       Set up a QCON, and call SET_TCON to allocate array numbers,
C****       set up scales, titles, etc. 
C**** QCON denotes when the conservation diags should be accumulated
C**** 1:NPTS+1 ==> INST,  DYN,   COND,   RAD,   PREC,   LAND,  SURF,
C****            FILTER,STRDG/OCEAN, DAILY, OCEAN1, OCEAN2,
C**** First 12 are standard for all tracers and GCM
      QCON=(/ t,                                           !instant.
     *        T,  T,  F,  T,  T,  T,  T,  T,  F,  F,  F,   !2-12 (npts)
     *        F,  F,  F,  F,  F,  F,  F,  F,  F,  F,       !13-22
     *        F,  F,  F,  F,  F,  F,  F,  F,  F,  F  /)   !21-ktcon-1
      do n=1,ntm
        kt_power_inst(n)   = ntm_power(n)+2
        kt_power_change(n) = ntm_power(n)-4
        scale_inst(n)   = 10.**(-kt_power_inst(n))
        scale_change(n) = 10.**(-kt_power_change(n))
        inst_unit(n) = unit_string(kt_power_inst(n),  ' kg/m^2)')
         sum_unit(n) = unit_string(kt_power_change(n),' kg/s/m^2)')
      end do
      N = n_air
      CALL SET_TCON(QCON,TRNAME(N),inst_unit(n),
     *     sum_unit(n),scale_inst(n),scale_change(n), N,CONPTs)
      N = n_SF6
      qcon(13) = .true.; conpts(1) = 'L1 SOURCE'
      CALL SET_TCON(QCON,TRNAME(N),inst_unit(n),
     *     sum_unit(n),scale_inst(n),scale_change(n), N,CONPTs)
      N = n_Rn222
      qcon(13) = .true.; conpts(1) = 'L1 SOURCE'
      qcon(14) = .true.; conpts(2) = 'DECAY'
      CALL SET_TCON(QCON,TRNAME(N),inst_unit(n),
     *     sum_unit(n),scale_inst(n),scale_change(n), N,CONPTs)
      N = n_CO2
      qcon(13) = .true.; conpts(1) = 'FossilFuel'
      qcon(14) = .true.; conpts(2) = 'Fertilization'
      qcon(15) = .true.; conpts(3) = 'Forest Regrowth'
      qcon(16) = .true.; conpts(4) = 'Land Use'
      qcon(17) = .true.; conpts(5) = 'Ecosystem Exch'
      qcon(18) = .true.; conpts(6) = 'Ocean Exch'
      CALL SET_TCON(QCON,TRNAME(N),inst_unit(n),
     *     sum_unit(n),scale_inst(n),scale_change(n), N,CONPTs)

      return
      end subroutine init_tracer


      SUBROUTINE SET_TCON(QCON,NAME_CON,INST_UNIT,SUM_UNIT,INST_SC
     *     ,CHNG_SC, itr,CONPTs)
!@sum  SET_TCON assigns conservation diagnostic array indices
!@auth Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only: sday
      USE MODEL_COM, only: dtsrc,nfiltr
      USE TRACER_DIAG_COM, only: ktcon,title_tcon,scale_tcon,nsum_tcon
     *     ,nofmt,ia_tcon,name_tconsrv,lname_tconsrv,units_tconsrv
     *     ,ntcons
      USE DAGCOM, only: npts,ia_d5d,ia_d5s,ia_filt,ia_12hr,ia_src,CONPT
      IMPLICIT NONE
!@var QCON logical variable sets where conservation diags are saved
      LOGICAL, INTENT(IN),DIMENSION(ktcon-1) :: QCON
!@var INST_SC scale for instantaneous value
      REAL*8, INTENT(IN) :: INST_SC
!@var CHNG_SC scale for changes
      REAL*8, INTENT(IN) :: CHNG_SC
!@var NAME_CON name of conservation quantity
      CHARACTER*8, INTENT(IN) :: NAME_CON
!@var sname name of conservation quantity (no spaces)
      CHARACTER*8 :: sname
!@var INST_UNIT string for unit for instant. values
      CHARACTER*20, INTENT(IN) :: INST_UNIT
!@var SUM_UNIT string for unit for summed changes
      CHARACTER*20, INTENT(IN) :: SUM_UNIT
!@var ITR index for the tracer
      INTEGER, INTENT(IN) :: ITR
!@var CONPTS
      CHARACTER*16, DIMENSION(ntcons) :: CONPTS
      INTEGER NI,NM,NS,N,k

C**** remove spaces in NAME_CON for netcdf names
      sname=TRIM(NAME_CON)
      do k=1,len_trim(NAME_CON)
        if (sname(k:k).eq." ") sname(k:k)="_"
      end do
C****
      NI=1
      NOFMT(1,itr) = NI
      TITLE_TCON(NI,itr) = 
     *     "0INSTANT "//TRIM(NAME_CON)//" "//TRIM(INST_UNIT)
      SCALE_TCON(NI,itr) = INST_SC
      name_tconsrv(NI,itr) ="inst_"//sname
      lname_tconsrv(NI,itr) = "INSTANT "//TRIM(NAME_CON)
      units_tconsrv(NI,itr) = INST_UNIT
      IA_TCON(NI) = 12
      NM=NI
      DO N=2,ktcon-1
        IF (QCON(N)) THEN
          NM = NM + 1
          NOFMT(N,itr) = NM
          if (n.le.npts+1) then
            TITLE_TCON(NM,itr) = " CHANGE OF "//TRIM(NAME_CON)//" BY "//
     *         CONPT(N-1)
            name_tconsrv(NM,itr) =
     *           "chg_"//trim(sname)//"_"//TRIM(CONPT(N-1)(1:3))
          else
            TITLE_TCON(NM,itr) = " CHANGE OF "//TRIM(NAME_CON)//" BY "//
     *         CONPTs(N-npts-1)
            name_tconsrv(NM,itr) =
     *           "chg_"//trim(sname)//"_"//TRIM(CONPTs(N-npts-1)(1:3))
          end if
          lname_tconsrv(NM,itr) = TITLE_TCON(NM,itr)
          units_tconsrv(NM,itr) = SUM_UNIT
          SELECT CASE (N)
          CASE (2)
            SCALE_TCON(NM,itr) = CHNG_SC/DTSRC
            IA_TCON(NM) = ia_d5d
          CASE (3,4,5,6,7,9,11,12)
            SCALE_TCON(NM,itr) = CHNG_SC/DTSRC
            IA_TCON(NM) = ia_d5s
          CASE (8)
            SCALE_TCON(NM,itr) = CHNG_SC/(NFILTR*DTSRC)
            IA_TCON(NM) = ia_filt
          CASE (10)
            SCALE_TCON(NM,itr) = CHNG_SC*2./SDAY
            IA_TCON(NM) = ia_12hr
          CASE (13:)   ! special tracer sources
            SCALE_TCON(NM,itr) = CHNG_SC/DTSRC
            IA_TCON(NM) = ia_src
          END SELECT
        ELSE
          NOFMT(N,itr) = 0
        END IF
      END DO
      NS=NM+1
      IF (NS.gt.KTCON) THEN
        WRITE(6,*) "KTCON not large enough for extra conserv diags",
     *       KTCON,NI,NM,NS,NAME_CON
        STOP "Change KTCON in tracer diagnostic common block"
      END IF
      TITLE_TCON(NS,itr) = " SUM OF CHANGES "//TRIM(SUM_UNIT)
      name_Tconsrv(NS,itr) ="sum_chg_"//trim(sname)
      lname_Tconsrv(NS,itr) = " SUM OF CHANGES OF "//TRIM(NAME_CON)
      units_Tconsrv(NS,itr) = SUM_UNIT
      SCALE_TCON(NS,itr) = 1.
      IA_TCON(NS) = 12
      NSUM_TCON(NI,itr) = -1
      NSUM_TCON(NI+1:NS-1,itr) = NS
      NSUM_TCON(NS,itr) = 0
      RETURN
      END SUBROUTINE set_tcon


      SUBROUTINE tracer_IC(tname)
!@sum tracer_IC initializes tracers. Called once per day
!@auth Jean Lerner
      USE MODEL_COM, only: itime,lm
      USE GEOM, only: dxyp
      USE DYNAMICS, only: am  ! Air mass of each box (kg/m^2)
      USE TRACER_COM
      USE TRACER_DIAG_COM
      USE FILEMANAGER, only : openunit,closeunit
      IMPLICIT NONE
      INTEGER i,N,L,j,iu_data
      CHARACTER*8 tname
      CHARACTER*80 title
      REAL*4 CO2IC(IM,JM,LM)

      select case (tname)

      case default
        write(6,*) 'In TRACER_IC:',tname,' does not exist '
        stop
  
      case ('Air')
        n = n_air
        if (itime.eq.itime_tr0(n)) then
          CALL CALC_AM(lm)
          do l=1,lm
          do j=1,jm
            trm(:,j,l,n) = am(l,:,j)*dxyp(j)
          end do; enddo
          do i=2,im
            trm(i,1,:,n) =  trm(1,1,:,n)   !poles
            trm(i,jm,:,n) = trm(1,jm,:,n)   !poles
          enddo
          trmom(:,:,:,:,n) = 0.
        end if
      case ('SF6')
        n = n_SF6
        if (itime.eq.itime_tr0(n)) then
          trm(:,:,:,n) = 0.
          trmom(:,:,:,:,n) = 0.
        end if
      case ('Rn222')
        n = n_Rn222
        if (itime.eq.itime_tr0(n)) then
          CALL CALC_AM(lm)
          do l=1,lm
          do j=1,jm
            trm(:,j,l,n) = am(l,:,j)*dxyp(j)*1.d-22
          end do; end do
          do i=2,im
            trm(i,1,:,n) =  trm(1,1,:,n)   !poles
            trm(i,jm,:,n) = trm(1,jm,:,n)   !poles
          enddo
          trmom(:,:,:,:,n) = 0.
        end if
      case ('CO2')
        n = n_CO2
        if (itime.eq.itime_tr0(n)) then
          trm(:,:,:,n) = 0.
          trmom(:,:,:,:,n) = 0.
          CALL CALC_AM(lm)
          call openunit('CO2_IC',iu_data,.true.,.true.)
          read (iu_data) title,co2ic
          call closeunit(iu_data)
          write(6,*) title,' read from CO2_IC'
          do l=1,lm         !ppmv==>ppmm
          do j=1,jm
            trm(:,j,l,n) = co2ic(:,j,l)*am(l,:,j)*dxyp(j)*1.54d-6
          enddo; enddo
        end if
      END SELECT
      write(6,*) ' Tracer ',trname(n),' initialized at itime=',itime
      END SUBROUTINE tracer_IC


      MODULE CO2_SOURCES
      USE TRACER_COM
      real*8 co2_src(im,jm,6)
      END MODULE CO2_SOURCES


      subroutine read_CO2_sources(nt,iact)
!@sum Read in CO2 sources and sinks
!@auth Jean Lerner
C****
C**** There are two monthly sources and 4 annual sources
C**** Annual sources are read in at start and re-start of run only
C**** Monthly sources are interpolated each day
      USE MODEL_COM, only: itime,jday,DTsrc,JDperY,im,jm
      USE CONSTANT, only: hrday
      USE FILEMANAGER, only : openunit,closeunit, openunits,closeunits
      USE TRACER_COM, only: itime_tr0,trname
      use CO2_SOURCES, only: src=>co2_src
      implicit none
      character*80 title
      logical :: ifirst=.true.
      data jdlast /0/
      real*8 tlca(im,jm,2),tlcb(im,jm,2)  ! for 2 monthly sources
      integer*4 idofm(0:13)
      data idofm /-15,16,47,75,106,136,167,197,228,259,289,320,350,381/
      real*8 adj(6)
c     data adj/1.038,1.,1.,5.33,1.,1.749/          ! c
      data adj/3.81d0,3.67d0,3.67d0,19.54d0,3.67d0,6.42d0/     ! co2
      integer ann_units(4),mon_units(2)
      character*12 :: ann_files(4) = (/'CO2_FOS_FUEL','CO2_FERT    ',
     *   'CO2_REGROWTH','CO2_LAND_USE'/)
      logical :: ann_bins(4)=(/.true.,.true.,.true.,.true./)
      character*9 :: mon_files(2) = (/'CO2_VEG  ','CO2_OCEAN'/)
      logical :: mon_bins(2)=(/.true.,.true./)
      real*8 frac
      integer i,j,nt,iact,imon,iu,iumf,iuml,ioff,k,jdlast
      save ifirst,jdlast

      if (itime.lt.itime_tr0(nt)) return
 !    write(0,*) ' debugging traces iact,jday=',iact,jday,jdlast,ifirst
 !    write(6,*) ' debugging traces iact,jday=',iact,jday,jdlast,ifirst
 !    call sys_flush(6)
C****
C**** Annual Sources and sink
C**** Apply adjustment factors to bring sources into balance
C**** Annual sources are in KG C/M2/Y
C**** Sources are added once/hr; convert kg/year to kg
C****
      if (ifirst) then
        call openunits(ann_files,ann_units,ann_bins,4)
        k = 0
        do iu = ann_units(1),ann_units(4)
          k = k+1
          call readt (iu,0,src(1,1,k),im*jm,src(1,1,k),1)
          src(:,:,k) = src(:,:,k)*adj(k)/(hrday*JDperY)
        end do
        call closeunits(ann_units,4)
      endif
C****
C**** Monthly sources are interpolated to the current day
C****
      if (iact.eq.1 .and. .not.ifirst) return
      ifirst = .false.
      call openunits(mon_files,mon_units,mon_bins,2)
      iumf = mon_units(1)
      iuml = mon_units(2)
      if (jdlast.EQ.0) then ! NEED TO READ IN FIRST MONTH OF DATA
        imon=1          ! imon=January
        if (jday.le.16)  then ! JDAY in Jan 1-15, first month is Dec
          ioff = 0
          do iu=iumf,iuml
          ioff = ioff+1
          call readt(iu,0,tlca(1,1,ioff),im*jm,tlca(1,1,ioff),12)
          rewind iu
          end do
        else            ! JDAY is in Jan 16 to Dec 16, get first month
  120     imon=imon+1
          if (jday.gt.idofm(imon) .AND. imon.le.12) go to 120
          ioff = 0
          do iu=iumf,iuml
          ioff = ioff+1
          call readt(iu,0,tlca(1,1,ioff),im*jm,tlca(1,1,ioff),imon-1)
          if (imon.eq.13)  rewind iu
          end do
        end if
      else              ! Do we need to read in second month?
        if (jday.ne.jdlast+1) then ! Check that data is read in daily
          if (jday.ne.1 .OR. jdlast.ne.365) then
            write(6,*)
     *      'Incorrect values in TSRC: JDAY,JDLAST=',JDAY,JDLAST
            stop
          end if
          imon=imon-12  ! New year
          go to 130
        end if
        if (jday.le.idofm(imon)) go to 130
        imon=imon+1     ! read in new month of data
        tlca(:,:,:) = tlcb(:,:,:)
        if (imon.eq.13)  then
          do iu=iumf,iuml
            rewind iu
          end do
        end if
      end if
      ioff = 0
      do iu=iumf,iuml
        ioff = ioff+1
        call readt(iu,0,tlcb(1,1,ioff),im*jm,tlcb(1,1,ioff),1)
      end do
  130 jdlast=jday
c**** Interpolate two months of data to current day
      frac = float(idofm(imon)-jday)/(idofm(imon)-idofm(imon-1))
      do k=5,6
        src(:,:,k) = tlca(:,:,k-4)*frac + tlcb(:,:,k-4)*(1.-frac)
      end do
      write(6,*) trname(nt),'Sources interpolated to current day',frac
      call closeunits(mon_units,2)
C****
C**** Apply adjustment factors to bring sources into balance
C****
      do k=5,6
        src(:,:,k) = src(:,:,k)*adj(k)
      end do
C**** Monthly sources are in KG C/M2/S
C**** Sources are added once/hr; convert kg/sec to kg
      src(:,:,5:6) = src(:,:,5:6)*DTsrc
      return
      end subroutine read_CO2_sources



      SUBROUTINE tracer_source(tname)
!@sum tracer_source adds sources and sinks to tracers
!@auth Jean Lerner
      USE MODEL_COM, only: FEARTH,itime,dtsrc,JDperY,fland,itoice,ftype
      USE GEOM, only: dxyp
      USE QUSDEF
      USE DYNAMICS, only: am  ! Air mass of each box (kg/m^2)
      USE TRACER_COM
      USE TRACER_DIAG_COM
      USE CONSTANT, only: tf,sday,hrday
      USE PBLCOM, only: tsavg
      USE CO2_SOURCES, only: co2_src
      implicit none
      double precision, dimension(im,jm,lm) :: told
      integer :: i,j,k,l,ky,naij,najl,n
      double precision :: source,sarea,steppy,base,steppd,x,airm,anngas,
     *  steph,stepx,stepp,tsum
      character*8 tname

      select case (tname)

      case default
!     write(6,*) ' Sources for ',tname,' tracer are not in this routine'

C****
C**** Surface Sources of SF6 (Same grid as CFC)
C****
      case ('SF6')
      n = n_SF6
      if (itime.lt.itime_tr0(n)) return
C**** Source increases each year by .3pptv/year
C**** Distribute source over ice-free land
    ! STEPPY = NDYN*DT/(86400.*365.)
      steppy = DTsrc/(sday*JDperY)
C     Make sure index KY=1 in year that tracer turns on
    ! KY = 1 + (TAU-TAUTR0(N))/8760
      ky = 1 + (itime-itime_tr0(n))/(hrday*JDperY)
      WRITE (6,'(A,I2,A,I9)')' >> KY FOR SF6 IS',KY,' AT itime',itime
      base = (0.3d-12)*146./28.    !pptm
      x = base*ky*steppy
      airm = 9928.64 * 51018.7e10       !(kg/m**2 X m**2 = kg)
      anngas = x*airm/steppy
      told(:,:,1) = trm(:,:,1,n)
C**** Source over United States and Canada
      source = .37*anngas*steppy
      sarea  = 0.
      do j=31,35
      do i=12,22
        sarea = sarea + dxyp(j)*fearth(i,j)
      enddo; enddo
      do j=31,35
      do i=12,22
        trm(i,j,1,n) = trm(i,j,1,n) + source*dxyp(j)*fearth(i,j)/sarea
        trmom(mz, i,j,1,n) = trmom(mz ,i,j,1,n) 
     &                          - 1.5*source*dxyp(j)*fearth(i,j)/sarea
        trmom(mzz,i,j,1,n) = trmom(mzz,i,j,1,n) 
     &                           + .5*source*dxyp(j)*fearth(i,j)/sarea
      enddo; enddo
C**** Source over Europe and Russia
c     source = .37*anngas*steppy
      sarea  = 0.
      do j=33,39
      do i=35,45
        sarea = sarea + dxyp(j)*fearth(i,j)
      enddo; enddo
      do j=33,39
      do i=35,45
        trm(i,j,1,n) = trm(i,j,1,n) + source*dxyp(j)*fearth(i,j)/sarea
        trmom(mz, i,j,1,n) = trmom(mz ,i,j,1,n) 
     &                          - 1.5*source*dxyp(j)*fearth(i,j)/sarea
        trmom(mzz,i,j,1,n) = trmom(mzz,i,j,1,n) 
     &                           + .5*source*dxyp(j)*fearth(i,j)/sarea
      enddo; enddo
C**** Source over Far East
      source = .13*anngas*steppy
      sarea  = 0.
      do j=29,34
      do i=61,66
        sarea = sarea + dxyp(j)*fearth(i,j)
      enddo; enddo
      do j=29,34
      do i=61,66
        trm(i,j,1,n) = trm(i,j,1,n) + source*dxyp(j)*fearth(i,j)/sarea
        trmom(mz, i,j,1,n) = trmom(mz ,i,j,1,n) 
     &                          - 1.5*source*dxyp(j)*fearth(i,j)/sarea
        trmom(mzz,i,j,1,n) = trmom(mzz,i,j,1,n) 
     &                           + .5*source*dxyp(j)*fearth(i,j)/sarea
      enddo; enddo
C**** Source over Middle East
      source = .05*anngas*steppy
      sarea  = 0.
      do j=28,32
      do i=43,51
        sarea = sarea + dxyp(j)*fearth(i,j)
      enddo; enddo
      do j=28,32
      do i=43,51
        trm(i,j,1,n) = trm(i,j,1,n) + source*dxyp(j)*fearth(i,j)/sarea
        trmom(mz, i,j,1,n) = trmom(mz ,i,j,1,n) 
     &                          - 1.5*source*dxyp(j)*fearth(i,j)/sarea
        trmom(mzz,i,j,1,n) = trmom(mzz,i,j,1,n) 
     &                           + .5*source*dxyp(j)*fearth(i,j)/sarea
      enddo; enddo
C**** Source over South America
      source = .04*anngas*steppy
      i=27; j=18
      trm(i,j,1,n) = trm(i,j,1,n)             +     source*.5
      trmom(mz, i,j,1,n) = trmom(mz, i,j,1,n) - 1.5*source*.5
      trmom(mzz,i,j,1,n) = trmom(mzz,i,j,1,n) + 0.5*source*.5
      i=28; j=18
      trm(i,j,1,n) = trm(i,j,1,n)             +     source*.5
      trmom(mz, i,j,1,n) = trmom(mz, i,j,1,n) - 1.5*source*.5
      trmom(mzz,i,j,1,n) = trmom(mzz,i,j,1,n) + 0.5*source*.5
C**** Source over South Africa
      source = .02*anngas*steppy
      i=42; j=17
      trm(i,j,1,n) = trm(i,j,1,n)             +     source*.5
      trmom(mz, i,j,1,n) = trmom(mz, i,j,1,n) - 1.5*source*.5
      trmom(mzz,i,j,1,n) = trmom(mzz,i,j,1,n) + 0.5*source*.5
C**** Source over Australia and New Zealand
c     source = .02*anngas*steppy
      i=66; j=15
      trm(i,j,1,n) = trm(i,j,1,n)             +     source*.5
      trmom(mz, i,j,1,n) = trmom(mz, i,j,1,n) - 1.5*source*.5
      trmom(mzz,i,j,1,n) = trmom(mzz,i,j,1,n) + 0.5*source*.5
C**** diagnostics
      naij = ijts_SF6_source
      taijs(:,:,naij) = taijs(:,:,naij) + (trm(:,:,1,n)-told(:,:,1))
      najl = jls_SF6_source
      do j=1,jm
        tajls(j,1,najl) = tajls(j,1,najl)+sum(trm(:,j,1,n)-told(:,j,1))
      end do
      call DIAGTCA(13,n)

C****
C**** Sources and sinks for Radon-222
C****
      case ('Rn222')
      n = n_Rn222
      if (itime.lt.itime_tr0(n)) return
C**** ground source
    ! STEPPD = NDYN*DT/86400.
      steppd = DTsrc/sday
      told(:,:,1) = trm(:,:,1,n)
      do j=1,jm
      do i=1,im
C**** source from ice-free land
      if(tsavg(i,j).lt.tf)  then  !composit surface air temperature
        trm(i,j,1,n) = trm(i,j,1,n) + 1.0d-16*steppd*dxyp(j)*fearth(i,j)
        trmom(mz, i,j,1,n) = trmom(mz ,i,j,1,n) 
     &                          - 1.5*1.0d-16*steppd*dxyp(j)*fearth(i,j)
        trmom(mzz,i,j,1,n) = trmom(mzz,i,j,1,n) 
     &                          + 0.5*1.0d-16*steppd*dxyp(j)*fearth(i,j)
      else
        trm(i,j,1,n) = trm(i,j,1,n) + 3.2d-16*steppd*dxyp(j)*fearth(i,j)
        trmom(mz, i,j,1,n) = trmom(mz ,i,j,1,n) 
     &                          - 1.5*3.2d-16*steppd*dxyp(j)*fearth(i,j)
        trmom(mzz,i,j,1,n) = trmom(mzz,i,j,1,n) 
     &                          + 0.5*3.2d-16*steppd*dxyp(j)*fearth(i,j)
      end if
C**** source from ice-free ocean
      trm(i,j,1,n) = trm(i,j,1,n) +     1.6d-18*steppd*dxyp(j)*
     *               (1.-fland(i,j))*(1.-ftype(itoice,i,j))
      trmom(mz,i,j,1,n) = trmom(mz,i,j,1,n)
     &                            - 1.5*1.6d-18*steppd*dxyp(j)*
     *               (1.-fland(i,j))*(1.-ftype(itoice,i,j))
      trmom(mzz,i,j,1,n)= trmom(mzz,i,j,1,n)
     &                            + 0.5*1.6d-18*steppd*dxyp(j)*
     *               (1.-fland(i,j))*(1.-ftype(itoice,i,j))
      enddo  !i
      enddo  !j
      naij = ijts_Rn222_source
      taijs(:,:,naij) = taijs(:,:,naij) + (trm(:,:,1,n)-told(:,:,1))
      najl = jls_Rn222_source
      do j=1,jm
        tajls(j,1,najl) = tajls(j,1,najl)+sum(trm(:,j,1,n)-told(:,j,1))
      end do
      call DIAGTCA(13,n)

C**** radiative decay
      x = (1.-2.1d-6*DTsrc)
      told(:,:,:)  = trm(:,:,:,n)
      trm(:,:,:,n) = trm(:,:,:,n)*x
      trmom(:,:,:,:,n) = trmom(:,:,:,:,n)*x
      najl = jls_RN222_sink
      do l=1,lm
      do j=1,jm
        tajls(j,l,najl) = tajls(j,l,najl)+sum(trm(:,j,l,n)-told(:,j,l))
      enddo; enddo
      call DIAGTCA(14,n)

C****
C**** Sources and sinks for CO2
C****   Sources are either /sec or /year
C****
      case ('CO2')
      n = n_CO2
      if (itime.lt.itime_tr0(n)) return
      naij = ijts_CO2_fossil_fuel-1
      najl = jls_CO2_fossil_fuel-1
      do k=1,6
        do j=1,jm
          tsum = 0.
          do i=1,im
            x = co2_src(i,j,k)*dxyp(j)
            trm(i,j,1,n) = trm (i,j,1,n)+x
            trmom(mz, i,j,1,n) = trmom(mz, i,j,1,n)-x*1.5d0
            trmom(mzz,i,j,1,n) = trmom(mzz,i,j,1,n)+x*0.5d0
            tsum = tsum+x
            taijs(i,j,k+naij) = taijs(i,j,k+naij)+x
          end do
          tajls(j,1,k+najl) = tajls(j,1,k+najl)+tsum
        end do
        call DIAGTCA(npts+1+k,n)   ! conservation
      end do

      END SELECT
      END SUBROUTINE tracer_source


