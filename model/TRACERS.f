#include "rundeck_opts.h"

#ifdef TRACERS_ON
!@sum  TRACERS: tracer routines that are tracer-dependent
!@+    Routines included: 
!@+      Diagnostic specs: TRACER_DIAG_COM, INIT_TRACER, SET_TCON
!@+      Tracer definitions: TRACER_IC, set_tracer_source
!@+      Tracer independent routines: apply_tracer_source, tdecay
!@auth Jean Lerner/Gavin Schmidt
!@ver  1.0

      MODULE TRACER_DIAG_COM
!@sum Tracer diagnostic arrays
!@auth Jean Lerner
!ver   1.0
      USE TRACER_COM, only: ntm, ntsrcmax
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
!@var ijts_source tracer independent array for TAIJS surface src. diags
      INTEGER ijts_source(ntsrcmax,ntm)
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
!@var ijts_index: tracer index associated with a TAIJS diagnostic
      INTEGER, DIMENSION(ktaijs) :: ijts_index

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
!@var jls_source tracer independent array for TAJLS surface src. diags
      INTEGER jls_source(ntsrcmax,ntm)
!@var jls_decay tracer independent array for radioactive sinks
      INTEGER, DIMENSION(NTM) :: jls_decay
!@var jls_grav tracer independent array for grav. settling sink
      INTEGER, DIMENSION(NTM) :: jls_grav
!@var jls_index: Number of source/sink tracer JL diagnostics/tracer 
      integer, dimension(ktajls) :: jls_index
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
      INTEGER, DIMENSION(ktcon,ntm) ::  IA_TCON
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
!@var itcon_surf Index array for surface source/sink conservation diags
      INTEGER, DIMENSION(NTSRCMAX,NTM) :: itcon_surf
!@var itcon_decay Index array for decay conservation diags
      INTEGER, DIMENSION(NTM) :: itcon_decay
!@var itcon_grav Index array for gravitational settling conserv. diags
      INTEGER, DIMENSION(NTM) :: itcon_grav

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
      use DAGCOM, only: ia_src,ia_12hr
      USE MODEL_COM, only: dtsrc
      use BDIJ, only: ir_log2
      USE TRACER_COM
      USE TRACER_DIAG_COM
      USE CONSTANT, only: mair,sday
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
        scale_jls(k) = 1./DTsrc
      end do
      jls_index(:) = 0
        
      k = 0
      n = n_air  !no special sources

      n = n_SF6
        k = k + 1 
        jls_source(1,n)=k
        sname_jls(k) = 'Layer_1_source_of_'//trname(n)
        lname_jls(k) = 'SF6 CFC-GRID SOURCE, LAYER 1'
        jls_index(k) = n
        jls_ltop(k) = 1
        jls_power(k) = -3.
        units_jls(k) = unit_string(jls_power(k),' kg/s')
      n = n_Rn222
        k = k + 1 
        jls_decay(n) = k   ! special array for all radioactive sinks
        sname_jls(k) = 'Decay_of_'//trname(n)
        lname_jls(k) = 'LOSS OF RADON-222 BY DECAY'
        jls_index(k) = n
        jls_ltop(k) = lm
        jls_power(k) = -11.
        units_jls(k) = unit_string(jls_power(k),' kg/s')
        k = k + 1 
        jls_source(1,n)=k
        sname_jls(k) = 'Ground_Source_of_'//trname(n)
        lname_jls(k) = 'RADON-222 SOURCE, LAYER 1'
        jls_index(k) = n
        jls_ltop(k) = 1
        jls_power(k) = -10.
        units_jls(k) = unit_string(jls_power(k),' kg/s')
! keep AIJ and AJL CO2 sources in same order !!
      n = n_CO2
        k = k + 1 
        jls_source(1,n)=k
        sname_jls(k) = 'Fossil_fuel_source_'//trname(n)
        lname_jls(k) = 'CO2 Fossil fuel source (Marland)'
        jls_index(k) = n
        jls_ltop(k) = 1
        jls_power(k) = 3
        units_jls(k) = unit_string(jls_power(k),' kg/s')
        k = k + 1 
        jls_source(2,n)=k
        sname_jls(k) = 'fertilization_sink_'//trname(n)
        lname_jls(k) = 'CO2 fertilization sink (Friedlingstein)'
        jls_index(k) = n
        jls_ltop(k) = 1
        jls_power(k) = 3
        units_jls(k) = unit_string(jls_power(k),' kg/s')
        k = k + 1 
        jls_source(3,n)=k
        sname_jls(k) = 'Northern_forest_regrowth_'//trname(n)
        lname_jls(k) = 'CO2 Northern forest regrowth sink'
        jls_index(k) = n
        jls_ltop(k) = 1
        jls_power(k) = 3
        units_jls(k) = unit_string(jls_power(k),' kg/s')
        k = k + 1 
        jls_source(4,n)=k
        sname_jls(k) = 'Land_Use_Modification_'//trname(n)
        lname_jls(k) = 'CO2 from Land use modification (Houton)'
        jls_index(k) = n
        jls_ltop(k) = 1
        jls_power(k) = 3
        units_jls(k) = unit_string(jls_power(k),' kg/s')
        k = k + 1 
        jls_source(5,n)=k
        sname_jls(k) = 'Ecosystem_exchange_'//trname(n)
        lname_jls(k) = 'CO2 Ecosystem exchange (Matthews)'
        jls_index(k) = n
        jls_ltop(k) = 1
        jls_power(k) = 3
        units_jls(k) = unit_string(jls_power(k),' kg/s')
        k = k + 1 
        jls_source(6,n)=k
        sname_jls(k) = 'Ocean_exchange_'//trname(n)
        lname_jls(k) = 'CO2 Ocean exchange'
        jls_index(k) = n
        jls_ltop(k) = 1
        jls_power(k) = 3
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
     &   'tij_defs: Increase ktaijk=',ktaijk,' to at least ',k
        stop 'ktaijk too small'
      end if

C**** Tracer sources and sinks
C**** Defaults for ijts (sources, sinks, etc.)
C**** This needs to be 'hand coded' depending on circumstances
      k = 1
        ijts_source(1,n_SF6) = k
        ijts_index(k) = n_SF6
        ia_ijts(k) = ia_src   
        lname_ijts(k) = 'SF6 Layer 1 SOURCE'
        sname_ijts(k) = 'SF6_CFC-GRID_SOURCE,_LAYER_1'
        ijts_power(k) = -15
        units_ijts(k) = unit_string(ijts_power(k),' kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
      k = k+1
        ijts_source(1,n_Rn222) = k
        ijts_index(k) = n_Rn222
        ia_ijts(k) = ia_src   
        lname_ijts(k) = 'Rn222 L 1 SOURCE'
        sname_ijts(k) = 'Radon-222_SOURCE,_Layer_1'
        ijts_power(k) = -21.
        units_ijts(k) = unit_string(ijts_power(k),' kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc
! keep AIJ and AJL CO2 sources in same order !!
      n = n_CO2
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
        scale_inst(n)   = 10d0**(-kt_power_inst(n))
        scale_change(n) = 10d0**(-kt_power_change(n))
        inst_unit(n) = unit_string(kt_power_inst(n),  ' kg/m^2)')
         sum_unit(n) = unit_string(kt_power_change(n),' kg/s/m^2)')
      end do
      N = n_air
      CALL SET_TCON(QCON,TRNAME(N),inst_unit(n),
     *     sum_unit(n),scale_inst(n),scale_change(n), N,CONPTs)
      N = n_SF6
      itcon_surf(1,N) = 13
      qcon(itcon_surf(1,N)) = .true.; conpts(1) = 'L1 SOURCE'
      CALL SET_TCON(QCON,TRNAME(N),inst_unit(n),
     *     sum_unit(n),scale_inst(n),scale_change(n), N,CONPTs)
      N = n_Rn222
      itcon_surf(1,N) = 13
      qcon(itcon_surf(1,N)) = .true.; conpts(1) = 'L1 SOURCE'
      itcon_decay(n) = 14
      qcon(itcon_decay(n)) = .true.; conpts(2) = 'DECAY'
      CALL SET_TCON(QCON,TRNAME(N),inst_unit(n),
     *     sum_unit(n),scale_inst(n),scale_change(n), N,CONPTs)
      N = n_CO2
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
      CALL SET_TCON(QCON,TRNAME(N),inst_unit(n),
     *     sum_unit(n),scale_inst(n),scale_change(n), N,CONPTs)

C**** Here are some more examples of conservation diag configuration
c      n=n_dust
c      itcon_grav(n) = xx
c      qcon(itcon_grav(n)) = .true.; conpts(2) = 'SETTLING'
c      CALL SET_TCON(QCON,TRNAME(N),inst_unit(n),
c     *     sum_unit(n),scale_inst(n),scale_change(n), N,CONPTs)


C**** print out total tracer diagnostic array size
      WRITE (6,'(A14,2I8)') "KTACC=",KTACC

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
      IA_TCON(NI,itr) = 12
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
            IA_TCON(NM,itr) = ia_d5d
          CASE (3,4,5,6,7,9,11,12)
            SCALE_TCON(NM,itr) = CHNG_SC/DTSRC
            IA_TCON(NM,itr) = ia_d5s
          CASE (8)
            SCALE_TCON(NM,itr) = CHNG_SC/(NFILTR*DTSRC)
            IA_TCON(NM,itr) = ia_filt
          CASE (10)
            SCALE_TCON(NM,itr) = CHNG_SC*2./SDAY
            IA_TCON(NM,itr) = ia_12hr
          CASE (13:)   ! special tracer sources
            SCALE_TCON(NM,itr) = CHNG_SC/DTSRC
            IA_TCON(NM,itr) = ia_src
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
      IA_TCON(NS,itr) = 12
      NSUM_TCON(NI,itr) = -1
      NSUM_TCON(NI+1:NS-1,itr) = NS
      NSUM_TCON(NS,itr) = 0
      RETURN
      END SUBROUTINE set_tcon


      SUBROUTINE tracer_IC
!@sum tracer_IC initializes tracers when they are first switched on
!@auth Jean Lerner
      USE MODEL_COM, only: itime,im,jm,lm
      USE GEOM, only: dxyp
      USE DYNAMICS, only: am  ! Air mass of each box (kg/m^2)
      USE PBLCOM, only : npbl,trabl
      USE TRACER_COM, only : ntm,trm,trmom,itime_tr0,trname,needtrs
      USE FILEMANAGER, only : openunit,closeunit
      IMPLICIT NONE
      INTEGER i,n,l,j,iu_data,ipbl,it
      CHARACTER*80 title
      REAL*4 co2ic(im,jm,lm)

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

      end select

C**** Initialise pbl profile if necessary
      if (needtrs(n)) then
        do it=1,4
          do ipbl=1,npbl
            trabl(ipbl,n,:,:,it) = trm(:,:,1,n)
          end do
        end do
      end if

      write(6,*) ' Tracer ',trname(n),' initialized at itime=',itime
      end if
      end do
C****
      end subroutine tracer_IC


      MODULE CO2_SOURCES
      USE TRACER_COM
!@var co2_src C02 surface sources and sinks (kg/s)
      real*8 co2_src(im,jm,6)
      END MODULE CO2_SOURCES

      subroutine daily_tracer(iact)
!@sum daily_tracer is called once a day for tracers
!@auth Jean Lerner
C**** Note this routine must always exist (but can be a dummy routine)
      USE TRACER_COM, only : ntm,trname
      IMPLICIT NONE
      INTEGER n,iact

C**** Initialize tracers here to allow for tracers that 'turn on'
C****  at any time
      call tracer_IC

C**** Tracer specific call for CO2
      do n=1,ntm
        if (trname(n).eq."CO2") call read_CO2_sources(n,iact)
      end do

      return
C****
      end subroutine daily_tracer


      subroutine read_CO2_sources(nt,iact)
!@sum reads in CO2 sources and sinks
!@auth Jean Lerner
C****
C**** There are two monthly sources and 4 annual sources
C**** Annual sources are read in at start and re-start of run only
C**** Monthly sources are interpolated each day
      USE MODEL_COM, only: itime,jday,DTsrc,JDperY,im,jm,idofm=>JDmidOfM
      USE CONSTANT, only: sday
      USE FILEMANAGER, only : openunit,closeunit, openunits,closeunits
      USE TRACER_COM, only: itime_tr0,trname
      use CO2_SOURCES, only: src=>co2_src
      implicit none
      character*80 title
      logical :: ifirst=.true.
      data jdlast /0/
      real*8 tlca(im,jm,2),tlcb(im,jm,2)  ! for 2 monthly sources
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
C**** Sources need to be kg/m^2 s; convert /year to /s
C****
      if (ifirst) then
        call openunits(ann_files,ann_units,ann_bins,4)
        k = 0
        do iu = ann_units(1),ann_units(4)
          k = k+1
          call readt (iu,0,src(1,1,k),im*jm,src(1,1,k),1)
          src(:,:,k) = src(:,:,k)*adj(k)/(sday*JDperY)
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
C**** Monthly sources are in KG C/M2/S => src in kg/m^2 s
C****
      do k=5,6
        src(:,:,k) = src(:,:,k)*adj(k)
      end do
      return
      end subroutine read_CO2_sources

      SUBROUTINE set_tracer_source
!@sum tracer_source calculates non-interactive sources for tracers
!@auth Jean Lerner/Gavin Schmidt
      USE MODEL_COM, only: FEARTH,itime,JDperY,fland,psf,pmtop
      USE GEOM, only: dxyp,areag
      USE QUSDEF
      USE DYNAMICS, only: am  ! Air mass of each box (kg/m^2)
      USE TRACER_COM
      USE FLUXES, only : trsource,tot_trsource
      USE SEAICE_COM, only : rsi
      USE CONSTANT, only: tf,sday,hrday,bygrav,mair
      USE PBLCOM, only: tsavg
      USE CO2_SOURCES, only: co2_src
      implicit none
      integer :: i,j,ns,l,ky,n
      double precision :: source,sarea,steppy,base,steppd,x,airm,anngas,
     *  steph,stepx,stepp

C**** All sources are saved as kg/s
      do n=1,ntm

      if (itime.lt.itime_tr0(n)) cycle
      select case (trname(n))

      case default
!     write(6,*) ' Sources for ',trname(n),' tracer are not in this routine'
C****
C**** Surface Sources of SF6 (Same grid as CFC)
C****
      case ('SF6')
        trsource(:,:,:,n)=0
C**** Source increases each year by .3pptv/year
C**** Distribute source over ice-free land
        steppy = 1./(sday*JDperY)
C     Make sure index KY=1 in year that tracer turns on
        ky = 1 + (itime-itime_tr0(n))/(hrday*JDperY)
        WRITE (6,'(A,I2,A,I9)')' >> KY FOR SF6 IS',KY,' AT itime',itime
        base = (0.3d-12)*tr_mm(n)/mair !pptm
        x = base*ky*steppy
        airm = (psf-pmtop)*100.*bygrav * AREAG !(kg/m**2 X m**2 = kg)
        anngas = x*airm/steppy
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
            if(tsavg(i,j).lt.tf)  then !composite surface air temperature
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
            do i=1,im
              trsource(i,j,ns,n) = co2_src(i,j,ns)*dxyp(j)
            end do
          end do
        end do

      end select

C**** tot_trsource is required for PBL calculations
      tot_trsource(:,:,n) = 0.
      do ns=1,ntsurfsrc(n)
        tot_trsource(:,:,n) = tot_trsource(:,:,n)+trsource(:,:,ns,n)
      end do

      end do
C****
      END SUBROUTINE set_tracer_source

      SUBROUTINE apply_tracer_source(dtsurf)
!@sum apply_tracer_source adds non-interactive surface sources to tracers
!@auth Jean Lerner/Gavin Schmidt
      USE MODEL_COM, only : jm
      USE TRACER_COM, only : ntm,trm,trmom,ntsurfsrc
      USE QUSDEF, only : mz,mzz
      USE FLUXES, only : trsource,tot_trsource
      USE TRACER_DIAG_COM, only : taijs,tajls,ijts_source,jls_source
     *     ,itcon_surf
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: dtsurf
      INTEGER n,ns,naij,najl,j

C**** This is tracer independent coding designed to work for all
C**** surface sources that do not vary as a function of the surface
C**** variables (non-interactive). With no daily varying sources setting 
C**** these sources can be done once a day, otherwise it should be
C**** every source time step. 

      do n=1,ntm
        do ns=1,ntsurfsrc(n)
          trm(:,:,1,n) = trm(:,:,1,n) + trsource(:,:,ns,n)*dtsurf
C**** diagnostics
          naij = ijts_source(ns,n)
          taijs(:,:,naij) = taijs(:,:,naij) + trsource(:,:,ns,n)*dtsurf
          najl = jls_source(ns,n)
          do j=1,jm
            tajls(j,1,najl) = tajls(j,1,najl)+sum(trsource(:,j,ns,n))
     *           *dtsurf
          end do
          call DIAGTCA(itcon_surf(ns,n),n)
        end do
C**** modify vertical moments        
        trmom( mz,:,:,1,n) = trmom( mz,:,:,1,n)-1.5*tot_trsource(:,:,n)
     *       *dtsurf
        trmom(mzz,:,:,1,n) = trmom(mzz,:,:,1,n)+0.5*tot_trsource(:,:,n)
     *       *dtsurf
      end do
C****
      RETURN
      END SUBROUTINE apply_tracer_source

      SUBROUTINE TDECAY
!@sum TDECAY decays radioactive tracers every source time step
!@auth Gavin Schmidt/Jean Lerner
      USE MODEL_COM, only : im,jm,lm,itime,dtsrc
      USE TRACER_COM, only : ntm,trm,trmom,trdecy,itime_tr0
#ifdef TRACERS_WATER
     *     ,trwm
#endif
      USE TRACER_DIAG_COM, only : tajls,jls_decay,itcon_decay
      IMPLICIT NONE
      real*8, save, dimension(ntm) :: expdec = 1.
      real*8, dimension(im,jm,lm) :: told
      integer, save :: ifirst=1
      integer n,najl,j,l

      if (ifirst.eq.1) then               
        do n=1,ntm
          if (trdecy(n).gt.0.0) expdec(n)=exp(-trdecy(n)*dtsrc)
        end do
        ifirst = 0
      end if

      do n=1,ntm
        if (trdecy(n).gt.0. .and. itime.ge.itime_tr0(n)) then
C**** Atmospheric decay
          told(:,:,:)=trm(:,:,:,n)
#ifdef TRACERS_WATER
     *               +trwm(:,:,:,n)
          trwm(:,:,:,n)=expdec(n)*trwm(:,:,:,n)
#endif
          trm(:,:,:,n)=expdec(n)*trm(:,:,:,n)
          trmom(:,:,:,:,n)=expdec(n)*trmom(:,:,:,:,n)
          najl = jls_decay(n)
          do l=1,lm
          do j=1,jm
          tajls(j,l,najl)=tajls(j,l,najl)+sum(trm(:,j,l,n)
#ifdef TRACERS_WATER
     *               +trwm(:,j,l,n)
#endif
     *               -told(:,j,l))
          enddo 
          enddo
          call DIAGTCA(itcon_decay(n),n)
        end if
      end do
C****
      return
      end

      subroutine checktr(subr)

      SUBROUTINE TRGRAV
!@sum TRGRAV gravitationally settles particular tracers 
!@auth Gavin Schmidt/Reha Cakmur
      USE CONSTANT, only : visc_air,grav
      USE MODEL_COM, only : im,jm,lm,itime,dtsrc,zatmo
      USE TRACER_COM, only : ntm,trm,trmom,itime_tr0,trradius,trpdens
      USE TRACER_DIAG_COM, only : tajls,jls_grav,itcon_grav
      USE FLUXES, only : trgrdep
      USE DYNAMICS, only : gz
      IMPLICIT NONE
      real*8, save, dimension(ntm) :: stokevdt = 0.
      integer, save :: ifirst=1
      real*8, dimension(im,jm,lm) :: told
      real*8 fgrfluxd,fgrfluxu
      integer n,najl,i,j,l

      if (ifirst.eq.1) then               
C**** Calculate settling velocity based on Stokes' Law using particle
C**** density and effective radius
        do n=1,ntm
          if (trradius(n).gt.0.0) stokevdt(n)=dtsrc*2.*grav*trpdens(n)
     *         *trradius(n)**2/(9.*visc_air)
        end do
        ifirst = 0
      end if

      do n=1,ntm
        if (trradius(n).gt.0. .and. itime.ge.itime_tr0(n)) then
C**** Gravitional settling 
          do l=1,lm
          do j=1,jm
          do i=1,imaxj(j)
            told(i,j,l)=trm(i,j,l,n)
C**** Calculate height differences using geopotential
            if (l.eq.1) then   ! layer 1 calc
              fgrfluxd=stokevdt(n)*grav/(gz(i,j,l)-zatmo(i,j))
              trgrdep(i,j)=fgrfluxd*trm(i,j,l,n)
            else               ! above layer 1
              fgrfluxd=stokevdt(n)*grav/(gz(i,j,l)-gz(i,j,l-1))
            end if
            if (l.lt.lm) then  ! below top layer
              fgrfluxu=stokevdt(n)*grav/(gz(i,j,l+1)-gz(i,j,l))
            else               ! top layer
              fgrfluxu=0.
            end if
            trm(i,j,l,n)=trm(i,j,l  ,n)*(1.-fgrfluxd)
     *                 + trm(i,j,l+1,n)*    fgrfluxu
            trmom(mz ,i,j,l,n)=trmom(mz ,i,j,l,n)*(1.-fgrfluxd)
            trmom(mzz,i,j,l,n)=trmom(mzz,i,j,l,n)*(1.-fgrfluxd)
            trmom(mzx,i,j,l,n)=trmom(mzx,i,j,l,n)*(1.-fgrfluxd)
            trmom(myz,i,j,l,n)=trmom(myz,i,j,l,n)*(1.-fgrfluxd)
          end do
          end do
          end do

          najl = jls_grav(n)
          do l=1,lm
          do j=1,jm
          tajls(j,l,najl)=tajls(j,l,najl)+sum(trm(:,j,l,n)-told(:,j,l))
          enddo 
          enddo
          call DIAGTCA(itcon_grav(n),n)
        end if
      end do
C****
      return
      end

!@sum  CHECKTR Checks whether tracer variables are reasonable
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : ls1,im,jm,lm
      USE DYNAMICS, only : am
      USE GEOM, only : dxyp,imaxj
      USE TRACER_COM
      IMPLICIT NONE
      LOGICAL QCHECKO
      INTEGER I,J,L,N,imax,jmax,lmax
      REAL*8 relerr, errmax
!@var SUBR identifies where CHECK was called from
      CHARACTER*6, INTENT(IN) :: SUBR

      do n=1,ntm
        CALL CHECK3(trm(1,1,1,n),IM,JM,LM,SUBR,trname(n)(1:3))
        CALL CHECK3(trmom(1,1,1,1,n),NMOM,IM,JM*LM,SUBR,
     *       'X'//trname(n)(1:2))
#ifdef TRACERS_WATER
        CALL CHECK3(trwm(1,1,1,n),IM,JM,LM,SUBR,'W'//trname(n)(1:2))
#endif

C**** check whether air mass is conserved
        
        if (trname(n).eq.'Air') then
          errmax = 0. ; lmax=0 ; imax=0 ; jmax=0 
          do l=1,lm
          do j=1,jm
          do i=1,imaxj(j)
            relerr=abs(trm(i,j,l,n)-am(l,i,j)*dxyp(j))/
     *           (am(l,i,j)*dxyp(j))
            if (relerr.gt.errmax) then
              lmax=l ; imax=i ; jmax=j ; errmax=relerr
            end if
          end do
          end do
          end do
          print*,"Relative error in air mass after ",subr,":",imax
     *         ,jmax,lmax,errmax,trm(imax,jmax,lmax,n),am(lmax,imax
     *         ,jmax)*dxyp(jmax)
        end if
      end do
      
      return
      end subroutine checktr

#endif
