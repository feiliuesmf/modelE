#include "rundeck_opts.h"

      MODULE TRDIAG_COM
!@sum Tracer diagnostic arrays
!@+    Mostly tracer independent, but this may depend on applications
!@auth Jean Lerner
!ver   1.0
      USE MODEL_COM, only: im,jm,lm
      USE DIAG_COM, only: npts !npts are conservation quantities
     &     ,sname_strlen,units_strlen,lname_strlen,jm_budg
#if (defined TRACERS_ON) || (defined TRACERS_OCEAN)
      USE TRACER_COM, only: ntm
#ifdef TRACERS_AEROSOLS_SOA
     &                     ,nsoa
#endif  /* TRACERS_AEROSOLS_SOA */
     *     , ntsurfsrcmax, nt3Dsrcmax
#ifdef TRACERS_AMP
     *     ,ntmAMP
      USE AERO_CONFIG,only: nbins
#endif
#endif   /* TRACERS_ON or OCEAN */
#ifdef NEW_IO
      use cdl_mod
#endif
      IMPLICIT NONE
      SAVE

C**** TAIJS  <<<< KTAIJS and IJTS_xx are Tracer-Dependent >>>>
C**** TAJLS  <<<< KTAJLS and JLS_xx are Tracer-Dependent >>>>

#ifdef TRACERS_ON
!@dbparam to_volume_MixRat: For printout of tracer concentration
!@+   to_volume_MixRat=1: printout is in Volume Mixing Ratio
!@+   to_volume_MixRat=0: printout is in Mass Mixing Ratio
      INTEGER, DIMENSION(NTM) :: to_volume_MixRat=0
!@dbparam to_conc: For printout of 3D tracer concentration in kg/m3
!@+   to_conc=0: printout is as defined by to_volume_MixRat
!@+   to_conc=1: printout is in kg/m3
      INTEGER, DIMENSION(NTM) :: to_conc=0
#endif  /* TRACERS_ON */
#if (defined TRACERS_WATER) || (defined TRACERS_OCEAN)
!@dbparam to_per_mil For printout of tracer concentration in permil
      INTEGER, DIMENSION(NTM) :: to_per_mil = 0
#endif
#ifdef TRACERS_ON
!@var MMR_to_VMR: converts tracer mass mixing ratio to volume mr
      REAL*8, DIMENSION(NTM) :: MMR_to_VMR

!!! WARNING: if new diagnostics are added, keep io_trdiag up-to-date !!!
C**** TAIJLN
!@var TAIJLN 3D tracer diagnostics (all tracers)
      real*8, allocatable, dimension(:,:,:,:) :: taijln
      real*8, allocatable, dimension(:,:,:,:) :: taijln_loc
!@var SNAME_IJT, UNITS_IJT: Names and units of lat-sigma tracer IJ diags
      character(len=sname_strlen), dimension(ntm) :: sname_ijt
      character(len=units_strlen), dimension(ntm) :: units_ijt
!@var LNAME_IJT: descriptions of tracer IJ diagnostics
      character(len=lname_strlen), dimension(ntm) ::
     &     lname_ijt = 'unused'
!@var SCALE_IJT: printout scaling factor for tracer IJ diagnostics
      REAL*8, dimension(ntm) :: scale_ijt
!@var IJTC_POWER: power of 10 used for tracer IJ concentration diags
      integer, dimension(ntm) :: ijtc_power
!@var IJTM_POWER: power of 10 used for tracer IJ mass unit diags
      integer, dimension(ntm) :: ijtm_power

C**** TAIJN
!@param KTAIJ number of 2D diags describing surface and column load along 
!@+   with wet and dry deposition
!@+   please just increase this if needed - don't bother with pp options
      integer, parameter :: ktaij=22

!@var IJT_XX names for taijn diagnostics
      integer tij_conc,tij_surf,tij_surfbv,tij_mass
!@var IJT_XX names for water-based taijn diagnostics
      integer tij_rvr,tij_seaice,tij_prec,tij_evap,tij_grnd,tij_lk1
     *     ,tij_lk2,tij_soil,tij_snow,tij_uflx,tij_vflx,tij_icocflx
     *     ,tij_tusi,tij_tvsi,tij_rvro,tij_icb
      integer tij_drydep,tij_gsdep ! TRACERS_DRYDEP

!@var TAIJN lat/lon tracer diagnostics (all tracers)
      real*8, allocatable, dimension(:,:,:,:)      :: taijn
      real*8, allocatable, dimension(:,:,:,:) :: taijn_loc
!@var SCALE_TIJ: printout scaling factor for tracer IJK diagnostics
      REAL*8, dimension(ktaij,ntm) :: scale_tij
!@var SNAME_TIJ,UNITS_TIJ: Names and units of lat-sigma tracer diags
      character(len=sname_strlen), dimension(ktaij,ntm) :: sname_tij
      character(len=units_strlen), dimension(ktaij,ntm) :: units_tij
!@var DNAME_TIJ, DENOM_TIJ: Short names, indices of tij denominators.
!@+   Currently, denom is specified along with the standard metadata,
!@+   and the dnames are looked up afterward.  Exception: denominators
!@+   that are not part of the taijn array must be specified using dname.
!@+   For now the denom index indicates which _tracer_ index is to be
!@+   used as a denominator.  Ratios of 2 accumulation types are
!@+   not yet needed.
      character(len=sname_strlen), dimension(ktaij,ntm) :: dname_tij=''
      integer, dimension(ktaij,ntm) :: denom_tij=0
!@var LNAME_TIJ: descriptions of tracer IJK diags
      character(len=lname_strlen), dimension(ktaij,ntm) ::
     &     lname_tij = 'unused'

C**** TAIJS  <<<< KTAIJS and IJTS_xx are Tracer-Dependent >>>>
!@var ijs_XXX index for diags not specific to a certain tracer
      INTEGER :: ijs_ai,ijs_isoprene,ijs_NO2_1030,ijs_NO2_1030c,
     &ijs_NO2_1330,ijs_NO2_1330c,ijts_Sdrydep,
     &ijts_clrsky=0,ijts_pocean=0

!@param KTAIJS number of special lat/lon tracer diagnostics
!@+   please just increase this if needed - don't bother with pp options
      INTEGER,PARAMETER :: ktaijs=1800

!@param MaxSubCl Maximum number of sub classes of tracers for rad. diagnostics
      INTEGER,PARAMETER :: MaxSubCl=4
!@param MaxDMc Maximum number of special wet depo diags for MC clouds
!@param MaxDLs Maximum number of special wet depo diags for LS clouds
      INTEGER,PARAMETER :: MaxDMc=6,MaxDLs=6
!@param MaxSpec Maximum number special diagnostics not associated with specific
!@+ tracer
      INTEGER,PARAMETER :: MaxSpec=3
!@dbparam diag_rad switches on/off comprehensive radiative diags for tracers
      INTEGER :: diag_rad=0 ! =off (default)
!@var TAIJS  lat/lon special tracer diagnostics; sources, sinks, etc.
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: TAIJS
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: TAIJS_loc
!@var ijts_source tracer independent array for TAIJS surface src. diags
      INTEGER ijts_source(ntsurfsrcmax,ntm)
!@var ijts_isrc tracer independent array for TAIJS interactive srf. src.
      INTEGER ijts_isrc(ntsurfsrcmax,ntm)
!@var ijts_aq tracer independent array for TAIJS aqueous change
      INTEGER ijts_aq(ntm)
!@var ijts_alb BC impact on snow/ice albedo, grain size,sw and lw radiation
      INTEGER ijts_alb(2)
!@var ijts_tau tracer independent array for TAIJS hydrated opt. thick.
      INTEGER ijts_tau(2,ntm)
!@var ijts_tausub index for TAIJS opt. thick. for tracer sub classes
      INTEGER ijts_tausub(2,Ntm,MaxSubCl)
!@var ijts_sqex index for TAIJS total extinction for 6 radiation bands
      INTEGER ijts_sqex(2,6,Ntm)
!@var ijts_sqexsub index for TAIJS total extinction for 6 radiation bands for
!@+   tracer sub classes
      INTEGER ijts_sqexsub(2,6,Ntm,MaxSubCl)
!@var ijts_sqsc index for TAIJS scattering extinction for 6 radiation bands
      INTEGER ijts_sqsc(2,6,Ntm)
!@var ijts_sqscsub index for TAIJS scattering extinction for 6 radiation
!@+   bands for tracer sub classes
      INTEGER ijts_sqscsub(2,6,Ntm,MaxSubCl)
!@var ijts_sqcb index for TAIJS sct asymmetry factor for 6 radiation bands
      INTEGER ijts_sqcb(2,6,Ntm)
!@var ijts_sqcbsub index for TAIJS sct asymmetry factor for 6 radiation bands
!@+   for tracer sub classes
      INTEGER ijts_sqcbsub(2,6,Ntm,MaxSubCl)
!@var ijts_fc tracer independent array for TAIJS SW/LW rad. forcings
      INTEGER ijts_fc(8,ntm)
#ifdef AUXILIARY_OX_RADF
!@var ijts_auxfc auxiliary Ox array for TAIJS SW/LW rad. forcings
      INTEGER ijts_auxfc(4)
#endif
!@var ijts_fcsub index for TAIJS SW/LW rad. forc. for tracer sub classes
      INTEGER ijts_fcsub(8,Ntm,MaxSubCl)
!@var ijts_spec index for TAIJS diags. not associated with single tracer
      INTEGER :: ijts_spec(MaxSpec)
!@var ijts_gasex index for TAIJS associated with atm/oc gas exchange
      INTEGER :: ijts_gasex(3,Ntm)
#ifdef TRACERS_AMP
!@var ijts_AMPpdf special diagnostic for not-transported tracers
      INTEGER ijts_AMPpdf(1,nbins)
#endif
!@var ijts_AMPe tracer independent array for emissions
      INTEGER ijts_AMPe(Ntm)
!@var ijts_AMPp tracer independent array for AMP processes
      INTEGER ijts_AMPp(7,Ntm)
!@var ijts_trdpmc indices of taijs special wet depo diags for MC clouds
      INTEGER :: ijts_trdpmc(MaxDMc,Ntm)
!@var ijts_trdpls indices of taijs special wet depo diags for LS clouds
      INTEGER :: ijts_trdpls(MaxDLs,Ntm)
!@var ijts_wet tracer independent array for TAIJS wet depo diagnostics
      INTEGER :: ijts_wet(Ntm)
!@var ijts_3Dsource tracer independent array for TAIJS 3D src. diags
      INTEGER ijts_3Dsource(nt3Dsrcmax,ntm)
!@var SNAME_IJTS, UNITS_IJTS: Names & units of lat-sigma tracer diags
      character(len=sname_strlen), dimension(ktaijs) :: sname_ijts
      character(len=units_strlen), dimension(ktaijs) :: units_ijts
!@var LNAME_IJTS: descriptions of tracer IJTS diags
      character(len=lname_strlen), dimension(ktaijs) ::
     &     lname_ijts = 'unused'
!@var DNAME_IJTS, DENOM_IJTS: Short names, indices of ijts denominators.
!@+   Currently, dname is specified along with the standard metadata and
!@+   the denom indices are looked up afterward.
      character(len=sname_strlen), dimension(ktaijs) :: dname_ijts=''
      integer, dimension(ktaijs) :: denom_ijts=0
!@var SCALE_IJTS: printout scaling factor for tracer IJTS diagnostics
      REAL*8, dimension(ktaijs) :: scale_ijts
!@var IA_IJTS: idacc-number for tracer source/sink IJ diags
      integer ia_ijts(ktaijs)
!@var ijts_power: power of 10 used for tracer IJ source/sink diags
      INTEGER, DIMENSION(ktaijs) :: ijts_power
!@var ijts_HasArea: does accumulation need to be divided by grid area?
      LOGICAL, DIMENSION(ktaijs) :: ijts_HasArea

C**** TAIJLS 3D special tracer diagnostics

!@param ktaijl number of TAIJLS tracer diagnostics;
      INTEGER, PARAMETER :: ktaijl=50
#ifdef ACCMIP_LIKE_DIAGS 
     &                            + 12
#endif
#ifdef SOA_DIAGS
     &                            + 12
#ifdef TRACERS_TERP
     &                            + 1
#endif  /* TRACERS_TERP */
     &                            + 16*nsoa
#endif  /* SOA_DIAGS */
!@var TAIJLS  3D tracer diagnostics (tracer dependent)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: TAIJLS
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: TAIJLS_loc
!@var SNAME_IJLT: Names of 3D tracer IJL diagnostics
      character(len=sname_strlen), dimension(ktaijl) :: sname_ijlt
!@var LNAME_IJLT,UNITS_IJLT: descriptions/units of 3D tracer diagnostics
      character(len=lname_strlen), dimension(ktaijl) ::
     &     lname_ijlt = 'unused'
      character(len=units_strlen), dimension(ktaijl) :: units_ijlt
!@var SCALE_IJLT: printout scaling factor for 3D tracer diagnostics
      REAL*8, dimension(ktaijl) :: scale_ijlt
!@var IR_IJLT: range index of IJL diagnostics
      integer, dimension(ktaijl) :: ir_ijlt
!@var IA_IJLT: accumulation index for IJL diagnostics
      integer, dimension(ktaijl) :: ia_ijlt
!@var ijlt_power: power of 10 used for tracer IJL 3D diags
      INTEGER, DIMENSION(ktaijs) :: ijlt_power
!@var ijlt_XXX diag names associated with 3D tracer special diags
      INTEGER :: ijlt_OH,ijlt_NO3,ijlt_HO2,ijlt_COp,ijlt_COd,
     & ijlt_Oxp,ijlt_Oxd,ijlt_CH4d,ijlt_OxpHO2,ijlt_OxpCH3O2,ijlt_OxpRO2
     & ,ijlt_OxlOH,ijlt_OxlHO2,ijlt_OxlALK,ijlt_phO1D,ijlt_pO1D,ijlt_pOH
     & ,ijlt_NOxLgt,ijlt_NOvmr,ijlt_NO2vmr,ijlt_JO1D,ijlt_JNO2
     & ,ijlt_JH2O2,ijlt_prodSO4aq,ijlt_prodSO4gs
#ifdef SOA_DIAGS
!@var ijlt_soa_changeL_isoprene gas-phase changeL of isoprene SOA (ug/m3)
#ifdef TRACERS_TERP
!@var ijlt_soa_changeL_terpenes gas-phase changeL of terpenes (ug/m3)
#endif  /* TRACERS_TERP */
!@var ijlt_soa_voc2nox VOC/NOx ratio (ppbC/ppb)
!@var ijlt_soa_pcp Total non-volatile SOA-absorbing mass (ug/m3)
!@var ijlt_soa_aerotot PCP plus SOA (g+a) (ug/m3)
!@var ijlt_soa_aerotot_gas Gas-phase semivolatile potential SOA (g+a) (ug/m3)
!@var ijlt_soa_xmf_isop Molar fraction of isoprene SOA
!@var ijlt_soa_xmf_apin Molar fraction of a-pinene SOA
!@var ijlt_soa_zcoef_isop Activity coefficient for isoprene SOA
!@var ijlt_soa_zcoef_apin Activity coefficient for a-pinene SOA
!@var ijlt_soa_meanmw Mean organic aerosol molecular weight (g/mol)
!@var ijlt_soa_iternum Total iterations for SOA calculations (count)
!@var ijlt_soa_M0 Final M0 value
!@var ijlt_soa_y0_ug_g gas-phase y0_ug (ug/m3)
!@var ijlt_soa_y0_ug_a aerosol-phase y0_ug (ug/m3)
!@var ijlt_soa_y_ug_g gas-phase y_ug (ug/m3)
!@var ijlt_soa_y_ug_a aerosol-phase y_ug (ug/m3)
!@var ijlt_soa_changeL_g_before gas-phase changeL before SOA (ug/m3)
!@var ijlt_soa_changeL_a_before aerosol-phase changeL before SOA (ug/m3)
!@var ijlt_soa_changeL_g_after gas-phase changeL after SOA (ug/m3)
!@var ijlt_soa_changeL_a_after aerosol-phase changeL after SOA (ug/m3)
!@var ijlt_soa_apartmass Effective apartmass
!@var ijlt_soa_kpart Partitioning coefficient (m3/ug)
!@var ijlt_soa_kp Final partitioning coefficient (m3/ug)
!@var ijlt_soa_soamass Final soamass value
!@var ijlt_soa_partfact Final partfact value
!@var ijlt_soa_evap Previous time-step soa that survived evaporation (ug/m3)
!@var ijlt_soa_cond Previous time-step soa precursor that condensed (ug/m3)
!@var ijlt_soa_chem Condensed soa precursor that was produced chemically at this time-step (ug/m3)
      integer :: ijlt_soa_changeL_isoprene
#ifdef TRACERS_TERP
      integer :: ijlt_soa_changeL_terpenes
#endif  /* TRACERS_TERP */
      integer :: ijlt_soa_voc2nox, ijlt_soa_pcp, ijlt_soa_aerotot,
     $     ijlt_soa_aerotot_gas, ijlt_soa_xmf_isop, ijlt_soa_xmf_apin,
     $     ijlt_soa_zcoef_isop, ijlt_soa_zcoef_apin, ijlt_soa_meanmw,
     $     ijlt_soa_iternum,  ijlt_soa_m0
      integer, dimension(nsoa) :: ijlt_soa_y0_ug_g,  ijlt_soa_y0_ug_a,
     $     ijlt_soa_y_ug_g, ijlt_soa_y_ug_a, ijlt_soa_changeL_g_before,
     $     ijlt_soa_changeL_a_before, ijlt_soa_changeL_g_after,
     $     ijlt_soa_changeL_a_after, ijlt_soa_apartmass, ijlt_soa_kpart,
     $     ijlt_soa_kp, ijlt_soa_soamass, ijlt_soa_partfact,
     $     ijlt_soa_evap, ijlt_soa_cond, ijlt_soa_chem
#endif  /* SOA_DIAGS */
#ifdef TRACERS_AMP
!@var ijlt_AMPext special diagnostic for not-transported tracers
!@var ijlt_AMPm tracer independent array for AMP modes
      INTEGER :: ijlt_AMPext(6),ijlt_AMPm(2,ntm)
#endif 
!@var ijlt_3Dtau 3D tracer independent array for hydrated opt. thick.
      INTEGER ijlt_3Dtau(ntm)

C**** TAJLN
!@param ktajl,ktajlx number of TAJL tracer diagnostics;
!@+          ktajlx includes composites
      INTEGER, PARAMETER :: ktajl=10, ktajlx=ktajl+2
!@var TAJLN  vertical tracer diagnostics (all tracers)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:)     :: TAJLN
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: TAJLN_loc
!@var jlnt_xx Names for TAJLN diagnostics
      INTEGER jlnt_conc,jlnt_mass,jlnt_nt_tot,jlnt_nt_mm,jlnt_vt_tot,
     &  jlnt_vt_mm,jlnt_mc,jlnt_turb,jlnt_lscond, jlnt_bebe,
     &  jlnt_bepb,jlnt_cldh2o

!@var SNAME_JLN: Names of lat-sigma tracer JL diagnostics
      character(len=sname_strlen), dimension(ktajlx,ntm) :: sname_jln
!@var LNAME_JLN,UNITS_JLN: descriptions/units of tracer JL diagnostics
      character(len=lname_strlen), dimension(ktajlx,ntm) ::
     &     lname_jln = 'unused'
      character(len=units_strlen), dimension(ktajlx,ntm) :: units_jln
!@var SCALE_JLQ: printout scaling factor for tracer JL diagnostics
      REAL*8, dimension(ktajlx) :: scale_jlq
!@var IA_JLQ,JGRID_JLQ: idacc-numbers,gridtypes for tracer JL diags
      integer, dimension(ktajlx) :: ia_jlq,jgrid_jlq
!@var JLQ_POWER: power of 10 used for tracer JL diagnostics
!@+        It is associated with a specific physical process
      integer jlq_power(ktajlx)
!@var scale_jln: Scale for jl maps
      REAL*8, DIMENSION(ntm) :: scale_jln

C**** TAJLS  <<<< KTAJLS and JLS_xx are Tracer-Dependent >>>>
!@param ktajls number of source/sink TAJLS tracer diagnostics;
!@+   please just increase this if needed - don't bother with pp options
      INTEGER,PARAMETER :: ktajls=1017

!@var jls_XXX index for non-tracer specific or special diags
      INTEGER jls_OHconk,jls_HO2con,jls_NO3
     *     ,jls_phot,jls_incloud(2,ntm), jls_OHcon,jls_H2Omr
     *     ,jls_N2O5sulf,jls_day,jls_COd,jls_COp,jls_Oxd,jls_Oxp
     *     ,jls_ClOcon,jls_H2Ocon,jls_H2Ochem,jls_OxdT,jls_OxpT

!@var TAJLS  JL special tracer diagnostics for sources, sinks, etc
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: TAJLS
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: TAJLS_loc
!@var jls_source tracer independent array for TAJLS surface src. diags
      INTEGER jls_source(ntsurfsrcmax,ntm)
!@var jls_isrc tracer indep. array for TAJLS interactive surface src. diags
      INTEGER jls_isrc(ntsurfsrcmax,ntm)
!@var jls_3Dsource tracer independent array for TAJLS 3D source diags
      INTEGER jls_3Dsource(nt3Dsrcmax,ntm)
!@var jls_decay tracer independent array for radioactive sinks
      INTEGER, DIMENSION(NTM) :: jls_decay
!@var jls_grav tracer independent array for grav. settling sink
      INTEGER, DIMENSION(NTM) :: jls_grav
!@var jls_prec tracer independent array for precipitation/wet dep
      INTEGER, DIMENSION(2,NTM) :: jls_prec
!@var jls_trdpmc indices of tajls special wet depo diags for MC clouds
      INTEGER :: jls_trdpmc(MaxDMc,Ntm)
!@var jls_trdpls indices of tajls special wet depo diags for LS clouds
      INTEGER :: jls_trdpls(MaxDLs,Ntm)
!@var jls_wet tracer independent array for wet deposition (for old dust)
      INTEGER, DIMENSION(NTM) :: jls_wet

!@var jls_spec index for TAJLS for special diagnostics not associated with
!@var jls_spec single tracer
      INTEGER :: jls_spec(MaxSpec)
!@var jwt_jls: Weighting index for jls diags 1=simple average, 2=by area
      integer, dimension(ktajls) :: jwt_jls
!@var SNAME_JLS: Names of lat-sigma tracer JL sources/sinks
      character(len=sname_strlen), dimension(ktajls) :: sname_jls
!@var LNAME_JLS,UNITS_JLS: descriptions/units of tracer JLS diags
      character(len=lname_strlen), dimension(ktajls) ::
     &     lname_jls = 'unused'
      character(len=units_strlen), dimension(ktajls) :: units_jls
!@var SCALE_JLS: printout scaling factor for tracer JLS diagnostics
      REAL*8, dimension(ktajls) :: scale_jls
!@var IA_JLS,JGRID_JLS: idacc-numbers,gridtypes for tracer JL diags
      integer, dimension(ktajls) :: ia_jls,jgrid_jls
!@var JLS_POWER: power of 10 used for tracer JLS diagnostics
      integer jls_power(ktajls)
!@var jls_ltop: Top layer for this diagnostic
      integer jls_ltop(ktajls)
#endif  /* TRACERS_ON */

#if (defined TRACERS_ON) || (defined TRACERS_OCEAN)
C**** TCONSRV
!@param NTCONS Maximum Number of special tracer conservation points
      INTEGER, PARAMETER :: ntcons=20
#ifdef TRACERS_AMP
     &                             +3
#endif
!@param KTCON total number of conservation diagnostics for tracers
      INTEGER, PARAMETER :: KTCON=npts+ntcons+2
!@param ntmxcon total number of conservation quantities
#ifdef TRACERS_SPECIAL_Shindell
C**** include some extra troposphere only ones
      INTEGER, PARAMETER :: ntmxcon = ntm+3
#else 
#ifdef TRACERS_OCEAN
      INTEGER, PARAMETER :: ntmxcon = ntm+20 ! arbitrary increase, may change
#else
      INTEGER, PARAMETER :: ntmxcon = ntm
#endif
#endif
!@var TCONSRV conservation diagnostics for tracers
      REAL*8, DIMENSION(JM_BUDG,ktcon,ntmxcon) :: TCONSRV,TCONSRV_loc 
!@var SCALE_TCON scales for tracer conservation diagnostics
      REAL*8, DIMENSION(ktcon,ntmxcon) :: SCALE_TCON
!@var TITLE_TCON titles for tracer conservation diagnostics
      CHARACTER*38, DIMENSION(ktcon,ntmxcon) :: TITLE_TCON
!@var IA_TCON IDACC numbers for tracer conservation diagnostics
      INTEGER, DIMENSION(ktcon,ntmxcon) ::  IA_TCON
!@var NSUM_TCON indices for summation of conservation diagnostics
      INTEGER, DIMENSION(ktcon,ntmxcon) :: NSUM_TCON
!@var NOFMT indices for TCONSRV array
      INTEGER, DIMENSION(ktcon,ntmxcon) :: NOFMT
!@var CONPTS names of special processes for tracer conservation diags
      CHARACTER*16, DIMENSION(ntcons) :: CONPTS
!@var kt_power_inst,kt_power_change: Exponents for tracer conservation
      INTEGER, DIMENSION(ntm) :: kt_power_inst,kt_power_change
!@var name_tconsrv,lname_tconsrv,units_tconsrv: for tracer conservation
      character(len=sname_strlen), dimension(ktcon,ntmxcon) ::
     &     name_tconsrv='unused'
      character(len=units_strlen), dimension(ktcon,ntmxcon) ::
     &     units_tconsrv
      character(len=lname_strlen), dimension(ktcon,ntmxcon) ::
     &     lname_tconsrv
!@var SCALE_INST,SCALE_CHANGE: Scale factors for tracer conservation
      REAL*8, dimension(ntmxcon) :: SCALE_INST,SCALE_CHANGE
!@var itcon_surf Index array for surface source/sink conservation diags
      INTEGER, DIMENSION(ntsurfsrcmax,ntmxcon) :: itcon_surf
!@var itcon_3Dsrc Index array for 3D source/sink conservation diags
      INTEGER, DIMENSION(nt3Dsrcmax,ntmxcon) :: itcon_3Dsrc
!@var itcon_decay Index array for decay conservation diags
      INTEGER, DIMENSION(ntmxcon) :: itcon_decay
!@var itcon_mc Index array for moist convection conserv. diags
      INTEGER, DIMENSION(ntmxcon) :: itcon_mc
!@var itcon_ss Index array for large-scale condensation conserv. diags
      INTEGER, DIMENSION(ntmxcon) :: itcon_ss
!@var itcon_amp Index array for microphysical processes diags
      INTEGER, DIMENSION(7,ntmxcon) :: itcon_amp
!@var itcon_amp Index array for microphysical processes diags
      INTEGER, DIMENSION(ntmxcon) :: itcon_ampe
!@var itcon_amp Index array for microphysical processes diags
      INTEGER, DIMENSION(2,ntmxcon) :: itcon_ampm
!@var itcon_dd Index array for dry deposition conserv. diags
      INTEGER, DIMENSION(ntmxcon,2) :: itcon_dd
!@var itcon_wt Index array for dust/mineral dust deposition conserv. diags
      INTEGER, DIMENSION(ntmxcon) :: itcon_wt
!@var natmtrcons, nocntrcons number of atmospheric/ocean tcon diags
      INTEGER :: natmtrcons=0, nocntrcons=0
#endif  /* TRACERS_ON  or  TRACERS_OCEAN */

!@var PDSIGJL temporary storage for mean pressures for jl diags
      REAL*8, DIMENSION(JM,LM)            :: PDSIGJL

#ifdef TRACERS_WATER
!@var TRP_acc, TRE_acc accumulation arrays for some SUBDD diags
!!    REAL*8 TRP_acc(ntm,IM,JM), TRE_acc(ntm,IM,JM)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:), public :: TRP_acc,TRE_acc
#endif
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_DUST)
!@var PM2p5_acc, PM10_acc accumulation arrays for some SUBDD diags
!@+ s prefix means SFC and l1 means L=1 accumulations. c prefix means
!@+ concentration units (kg/m3) rather than the mass mixing ratio
      REAL*8, ALLOCATABLE, DIMENSION(:,:), public ::  ! (IM,JM)
     &sPM2p5_acc,sPM10_acc,l1PM2p5_acc,l1PM10_acc,csPM2p5_acc,csPM10_acc
#endif

!@var trcsurf global array of tracer mixing ratio at surface [kg/kg]
      REAL*8,ALLOCATABLE,DIMENSION(:,:,:) :: trcsurf
!@var trcSurfByVol global array of tracer concentration at surface [kg/m^3]
      real(kind=8),allocatable,dimension(:,:,:) :: trcSurfByVol

!@var trcSurfMixR_acc global array of tracers to accumulate mxixing ratio at
!@var                 surface for subdd diagnostics [kg/kg]
!@var trcSurfByVol_acc global array of tracers to accumulate concentration at
!@var                  surface for subdd diagnostics [kg/m^3]
      real(kind=8),allocatable,dimension(:,:,:) :: trcSurfMixR_acc
     &     ,trcSurfByVol_acc

#ifdef NEW_IO
#ifdef TRACERS_ON
! This section declares arrays used to write tracer
! diagnostics accumulations in a format suitable for offline
! postprocessing.  The size of these arrays cannot be known
! precisely a priori due to the large number of outputs not
! declared in advance.  Thus, metadata arrays are declared with
! a larger-than-needed size, and acc arrays are re-allocated
! to the correct size after the size is determined.
! The "n" and "s" instances of diagnostics classes are
! merged as follows:
!@var taijl_out combines taijln and taijls.  Distributed.
!@var taij_out  combines taijn  and taijs.   Distributed.
!@var tajl_out  combines tajln  and tajls.
!@var tconsrv_out combines the ktcon/ntmxcon dims of tconsrv
!@+               and omits its unused elements.
! All contents of txxx_out arrays are PER UNIT AREA.
! Denominator information is stored in denom_xxx.
! Metadata for scaled outputs is consolidated in cdl_xxx.
! 
      integer, parameter :: ktaij_ = (ktaij*ntm+ktaijs
#ifdef TRACERS_SPECIAL_O18
     &     + 4                  ! include dexcess + D17O diags
#endif
#ifdef TRACERS_DRYDEP
     &     + ntm                ! include dry dep % diags
#endif
     &     )*3/2  ! make 50% larger for denoms and extra specials
      integer :: ktaij_out ! actual number of qtys in taij_out
      real*8, dimension(:,:,:), allocatable :: taij_out
      integer, dimension(ktaij_) :: ir_taij,ia_taij,denom_taij
      character(len=lname_strlen), dimension(ktaij_) :: lname_taij
      character(len=sname_strlen), dimension(ktaij_) :: sname_taij
      character(len=units_strlen), dimension(ktaij_) :: units_taij
      real*8, dimension(ktaij_) :: scale_taij
      type(cdl_type) :: cdl_taij,cdl_taij_latlon
      real*8, dimension(:,:,:), allocatable :: hemis_taij

      integer, parameter :: ktaijl_ = (ntm+ktaijl
#ifdef TRACERS_SPECIAL_O18
     *     + 2        ! include dexcess + D17O diags
#endif
     &     )*3/2  ! make 50% larger for denoms and extra specials
      integer :: ktaijl_out ! actual number of qtys in taijl_out
      real*8, dimension(:,:,:,:), allocatable :: taijl_out
      integer, dimension(ktaijl_) :: ir_taijl,ia_taijl,denom_taijl
      character(len=lname_strlen), dimension(ktaijl_) :: lname_taijl
      character(len=sname_strlen), dimension(ktaijl_) :: sname_taijl
      character(len=units_strlen), dimension(ktaijl_) :: units_taijl
      real*8, dimension(ktaijl_) :: scale_taijl
      type(cdl_type) :: cdl_taijl,cdl_taijl_latlon

      integer, parameter :: ktajl_ = (ktajlx*ntm+ktajls
#ifdef TRACERS_SPECIAL_O18
     &     + 4                  ! include dexcess + D17O diags
#endif
#ifdef TRACERS_COSMO
     &     + 2                  ! beryllium ratios
#endif
     &     )*3/2  ! make 50% larger for denoms and extra specials
      integer :: ktajl_out ! actual number of qtys in tajl_out
      real*8, dimension(:,:,:), allocatable :: tajl_out
      integer, dimension(ktajl_) :: pow_tajl,ia_tajl,denom_tajl,
     &     jgrid_tajl,lgrid_tajl,ltop_tajl
      character(len=lname_strlen), dimension(ktajl_) :: lname_tajl
      character(len=sname_strlen), dimension(ktajl_) :: sname_tajl
      character(len=units_strlen), dimension(ktajl_) :: units_tajl
      real*8, dimension(ktajl_) :: scale_tajl
      type(cdl_type) :: cdl_tajl
      real*8, dimension(:,:,:), allocatable :: hemis_tajl,vmean_tajl
#endif /* TRACERS_ON */

#if (defined TRACERS_ON) || (defined TRACERS_OCEAN)
      integer :: ktcon_out ! actual number of qtys in tconsrv_out
      real*8, dimension(:,:), allocatable :: tconsrv_out,hemis_tconsrv
      real*8, dimension(:), allocatable :: scale_tcon_out
      integer, dimension(:), allocatable ::  ia_tcon_out
      character(len=sname_strlen), dimension(:), allocatable ::
     &     sname_tconsrv_out
      type(cdl_type) :: cdl_tconsrv
#endif
#endif  /* NEW IO */

      END MODULE TRDIAG_COM

#ifdef TRACERS_ON
      SUBROUTINE SET_TCON(QCON,NAME_CON,QSUM,INST_UNIT,SUM_UNIT
     *     ,INST_SC,CHNG_SC, itr,CONPTs)
!@sum  SET_TCON assigns conservation diagnostic array indices
!@auth Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only: sday
      USE MODEL_COM, only: dtsrc,nfiltr
      USE DIAG_COM, only: npts,ia_d5d,ia_d5s,ia_filt,ia_12hr,ia_src
     *     ,conpt0
      USE TRDIAG_COM, only: ktcon,title_tcon,scale_tcon,nsum_tcon
     *     ,nofmt,ia_tcon,name_tconsrv,lname_tconsrv,units_tconsrv
     *     ,ntcons
      IMPLICIT NONE
!@var QCON denotes at which points conservation diags are saved
      LOGICAL, INTENT(IN),DIMENSION(ktcon-1) :: QCON
!@var QSUM sets whether each diag is included in final sum
!@+   should be zero for diags set using DIAGTCB (i.e. using difference)
      LOGICAL, INTENT(IN),DIMENSION(ktcon-1) :: QSUM
      LOGICAL, DIMENSION(ktcon-1) :: QSUM_CON   ! local version
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
!@var CONPTS conservation diag points for special tracer diags
      CHARACTER*16, INTENT(IN), DIMENSION(ntcons) :: CONPTS
!@var CONPTS_sname like CONPTS but without spaces
      CHARACTER*16, DIMENSION(ntcons) :: CONPTS_sname
!@var CONPT0_sname like CONPT0 but without spaces
      CHARACTER*10, DIMENSION(npts) :: CONPT0_sname
      CHARACTER*11 CHGSTR
      INTEGER NI,NM,NS,N,k
      CHARACTER*40 clean_str

C**** make nice netcdf names
      sname=trim(clean_str(name_con))
      do n=1,ntcons
         conpts_sname(n) = trim(clean_str(conpts(n)))
      enddo
      do n=1,npts
         conpt0_sname(n) = trim(clean_str(conpt0(n)))
      enddo
C****
      NI=1
      NOFMT(1,itr) = NI
      TITLE_TCON(NI,itr) =
     *     "0INSTANT "//TRIM(NAME_CON)//" "//TRIM(INST_UNIT)
      SCALE_TCON(NI,itr) = INST_SC
      name_tconsrv(NI,itr) ="inst_"//sname
      lname_tconsrv(NI,itr) = "INSTANT "//TRIM(NAME_CON)
      units_tconsrv(NI,itr) = INST_UNIT
      NSUM_TCON(NI,itr) = -1
      IA_TCON(NI,itr) = 12
      NM=NI
      DO N=2,ktcon-1
        IF (QCON(N)) THEN
          NM = NM + 1
          NOFMT(N,itr) = NM
          QSUM_CON(NM)=.FALSE.
          IF (QSUM(N)) QSUM_CON(NM)=.TRUE.
          CHGSTR=" CHANGE OF "
          if (n.le.npts+1) then
            TITLE_TCON(NM,itr) = CHGSTR//TRIM(NAME_CON)//" BY "//
     *         CONPT0(N-1)
            name_tconsrv(NM,itr) =
     *           "chg_"//trim(sname)//"_"//TRIM(CONPT0_sname(N-1))
          else
            IF (.not. QSUM(N)) CHGSTR="     DELTA "
            TITLE_TCON(NM,itr) = CHGSTR//TRIM(NAME_CON)//" BY "//
     *           CONPTs(N-npts-1)
            name_tconsrv(NM,itr) =
     *           "chg_"//trim(sname)//"_"//TRIM(CONPTs_sname(N-npts-1))
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
        call stop_model(
     &       "Change KTCON in tracer diagnostic common block",255)
      END IF
      DO NM=NI+1,NS-1
        NSUM_TCON(NM,itr) = -1
        IF (QSUM_CON(NM)) NSUM_TCON(NM,itr) = NS
      END DO
      TITLE_TCON(NS,itr) = " SUM OF CHANGES "//TRIM(SUM_UNIT)
      name_Tconsrv(NS,itr) ="sum_chg_"//trim(sname)
      lname_Tconsrv(NS,itr) = " SUM OF CHANGES OF "//TRIM(NAME_CON)
      units_Tconsrv(NS,itr) = SUM_UNIT
      SCALE_TCON(NS,itr) = 1.
      IA_TCON(NS,itr) = 12
      NSUM_TCON(NS,itr) = 0
      RETURN
      END SUBROUTINE set_tcon

      SUBROUTINE io_trdiag(kunit,it,iaction,ioerr)
!@sum  io_trdiag reads and writes tracer diagnostics arrays to file
!@auth Jean Lerner
!@ver  1.0
      USE MODEL_COM, only: im,jm,lm
      USE MODEL_COM, only: ioread,iowrite,iowrite_mon,iowrite_single
     *     ,irerun,ioread_single,lhead
      USE DIAG_COM, only : jm_budg
      USE DOMAIN_DECOMP_1D, only : AM_I_ROOT, ESMF_BCAST, grid, GET
      USE TRACER_COM, only: ntm
      USE TRDIAG_COM, only: taijln_loc, taijln, taijls_loc, taijls,
     *     taijn_loc,  taijn, taijs_loc,  taijs, tajln_loc,  tajln,
     *     tajls_loc,  tajls, tconsrv_loc, tconsrv, pdsigjl, ktaij,
     *     ktaijs, ktajlx, ktajls, ktcon, ktaijl, ntmxcon

      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var HEADER Character string label for individual records
      CHARACTER*80 :: HEADER, MODULE_HEADER = "TRdiag01"
!@var it input/ouput value of hour
      INTEGER, INTENT(INOUT) :: it
      INTEGER it_check ! =it if all diag. TA..4 were kept up-to-date
!@param KTACC total number of tracer diagnostic words
      INTEGER, PARAMETER ::
     *     ktacc=IM*JM*LM*NTM + IM*JM*LM*ktaijl + IM*JM*ktaij*NTM + IM
     *     *JM*ktaijs +JM_BUDG*LM*ktajlx*NTM + JM_BUDG*LM*ktajls +
     *     JM_BUDG*ktcon*ntmxcon
!@var TA..4(...) dummy arrays for reading diagnostics files
      real*4, allocatable, dimension(:,:,:,:) :: taijln4
      real*4, allocatable, dimension(:,:,:,:) :: taijls4
      real*4, allocatable, dimension(:,:,:,:) :: taijn4
      REAL*4, allocatable, DIMENSION(:,:,:)   :: TAIJS4
      REAL*4, allocatable, DIMENSION(:,:,:,:) :: TAJLN4
      REAL*4, allocatable, DIMENSION(:,:,:)   :: TAJLS4
      REAL*4, allocatable, DIMENSION(:,:,:)   :: TCONSRV4
      INTEGER :: status

      INTEGER :: J_0H, J_1H, I_0H, I_1H

      CALL GET( grid,  J_STRT_HALO = J_0H,  J_STOP_HALO = J_1H )
      I_0H=GRID%I_STRT_HALO
      I_1H=GRID%I_STOP_HALO

      write (MODULE_HEADER(lhead+1:80),'(a,i8,a)')
     *   'R8 TACC(',ktacc,'),it'

      SELECT CASE (IACTION)
      CASE (IOWRITE,IOWRITE_SINGLE)  
C***  PACK distributed arrays into global ones in preparation for output

        call gather_trdiag

        SELECT CASE (IACTION)
          CASE (IOWRITE,IOWRITE_MON) ! output to standard restart file
          IF (AM_I_ROOT())  WRITE (kunit,err=10) MODULE_HEADER,
     *                         TAIJLN,TAIJLS,TAIJN,TAIJS,
     *                         TAJLN ,TAJLS,TCONSRV,it
          CASE (IOWRITE_SINGLE)    ! output to acc file
            MODULE_HEADER(LHEAD+1:LHEAD+2) = 'R4'
            IF (AM_I_ROOT()) WRITE (kunit,err=10) MODULE_HEADER,
     *       REAL(TAIJLN,KIND=4),REAL(TAIJLS,KIND=4),
     *       REAL(TAIJN,KIND=4),
     *       REAL(TAIJS,KIND=4) ,REAL(TAJLN,KIND=4),
     *       REAL(TAJLS,KIND=4) ,REAL(TCONSRV,KIND=4),it
        END SELECT
      CASE (IOREAD:)          ! input from restart file
        SELECT CASE (IACTION)
        CASE (ioread_single)    ! accumulate diagnostic files
          if(am_i_root()) then
            ALLOCATE (taijln4(IM,JM,LM,ntm), stat=status )
            ALLOCATE (taijls4(IM,JM,LM,ktaijl), stat=status )
            ALLOCATE (taijn4(IM,JM,ktaij,ntm), stat=status )
            ALLOCATE (TAIJS4(IM,JM,ktaijs), stat=status )
            ALLOCATE (TAJLN4(JM_BUDG,LM,ktajlx,ntm), stat=status )
            ALLOCATE (TAJLS4(JM_BUDG,LM,ktajls), stat=status )
            ALLOCATE (TCONSRV4(JM_BUDG,ktcon,ntmxcon), stat=status )

            READ (kunit,err=10) HEADER,
     *           TAIJLN4,TAIJLS4,TAIJN4,TAIJS4,TAJLN4,TAJLS4,TCONSRV4
     *           ,it_check
            if (it.ne.it_check) then
              PRINT*,"io_trdiag: compare TAIJLN,TAIJLN4, ... dimensions"
              go to 10  ! or should this be just a warning ??
            end if
            IF (HEADER(1:LHEAD).NE.MODULE_HEADER(1:LHEAD)) THEN
              PRINT*,"Discrepancy in module version ",HEADER
     *             ,MODULE_HEADER
              GO TO 10
            END IF

C**** Accumulate diagnostics on global arrays (convert to real*8)
            TAIJLN = TAIJLN + TAIJLN4
            TAIJLS = TAIJLS + TAIJLS4
            TAIJN  = TAIJN  + TAIJN4
            TAIJS  = TAIJS  + TAIJS4
            TAJLN  = TAJLN  + TAJLN4
            TAJLS  = TAJLS  + TAJLS4
            TCONSRV= TCONSRV+ TCONSRV4
            
            DEALLOCATE ( taijln4 )
            DEALLOCATE ( taijls4 )
            DEALLOCATE ( taijn4 )
            DEALLOCATE ( TAIJS4 )
            DEALLOCATE ( TAJLN4 )
            DEALLOCATE ( TAJLS4 )
            DEALLOCATE ( TCONSRV4 )

          end if
C*** Unpack read global data into (real*8) local distributed arrays

          call scatter_trdiag


!         IF (HEADER(1:LHEAD).NE.MODULE_HEADER(1:LHEAD)) THEN
!           PRINT*,"Discrepancy in module version ",HEADER
!    *           ,MODULE_HEADER
!           GO TO 10
!         END IF
        CASE (ioread)  ! restarts
          if ( AM_I_ROOT() ) then
            READ (kunit,err=10) HEADER,
     *                           TAIJLN,TAIJLS,TAIJN,TAIJS,
     *                           TAJLN,TAJLS,TCONSRV,it
            IF (HEADER(1:LHEAD).NE.MODULE_HEADER(1:LHEAD)) THEN
              PRINT*,"Discrepancy in module version ",HEADER
     *             ,MODULE_HEADER
              GO TO 10
            END IF
          end if
C*** Unpack read global data into local distributed arrays

          call scatter_trdiag

          CALL ESMF_BCAST( grid, it )
        END SELECT
      END SELECT

      RETURN
 10   IOERR=1
      RETURN
      END SUBROUTINE io_trdiag
#endif

#ifdef NEW_IO
#ifdef TRACERS_ON /* only declare NEW_IO routines when needed */
      subroutine def_rsf_trdiag(fid,r4_on_disk)
!@sum  def_rsf_trdiag defines tracer diag array structure in restart+acc files
!@auth M. Kelley
!@ver  beta
      use trdiag_com, only :
     &     TAIJLN=>TAIJLN_loc,
     &     TAIJLS=>TAIJLS_loc,
     &     TAIJN=>TAIJN_loc,
     &     TAIJS=>TAIJS_loc,
     &     TAJLN,
     &     TAJLS,
     &     TAIJL=>TAIJL_out,
     &     TAIJ=>TAIJ_out,
     &     TAJL=>TAJL_out
      use domain_decomp_atm, only : grid
      use pario, only : defvar
      implicit none
      integer fid           !@var fid file id
      logical :: r4_on_disk !@var r4_on_disk if true, real*8 stored as real*4
      if(r4_on_disk) then ! acc file
        call defvar(grid,fid,taijl,
     &       'taijl(dist_im,dist_jm,lm,ktaijl)',r4_on_disk=.true.)
        call defvar(grid,fid,taij,
     &       'taij(dist_im,dist_jm,ktaij)',r4_on_disk=.true.)
        call defvar(grid,fid,tajl,
     &       'tajl(jm_budg,lm,ktajl)',r4_on_disk=.true.)
      else
        call defvar(grid,fid,taijln,'taijln(dist_im,dist_jm,lm,ntm)')
        call defvar(grid,fid,taijls,'taijls(dist_im,dist_jm,lm,ktaijl)')
        call defvar(grid,fid,taijs,'taijs(dist_im,dist_jm,ktaijs)')
        call defvar(grid,fid,taijn,'taijn(dist_im,dist_jm,ktaij,ntm)')
        call defvar(grid,fid,tajln,'tajln(jm_budg,lm,ktajlx,ntm)')
        call defvar(grid,fid,tajls,'tajls(jm_budg,lm,ktajls)')
      endif

      call def_rsf_tcons(fid,r4_on_disk)

      return
      end subroutine def_rsf_trdiag

      subroutine new_io_trdiag(fid,iaction)
!@sum  new_io_trdiag read/write tracer acc arrays from/to restart+acc files
!@auth M. Kelley
!@ver  beta new_ prefix avoids name clash with the default version
      use model_com, only : ioread,iowrite,iowrite_single
      use trdiag_com, only :
     &     TAIJLN=>TAIJLN_loc,
     &     TAIJLS=>TAIJLS_loc,
     &     TAIJN=>TAIJN_loc,
     &     TAIJS=>TAIJS_loc,
     &     TAJLN,
     &     TAJLS,
     &     TAIJL=>TAIJL_out,
     &     TAIJ=>TAIJ_out,
     &     TAJL=>TAJL_out
      use domain_decomp_atm, only : grid
      use pario, only : write_dist_data,read_dist_data,
     &     write_data,read_data
      implicit none
      integer fid   !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
      select case (iaction)
      case (iowrite_single)     ! output to acc file
        call write_dist_data(grid,fid,'taijl',taijl)
        call write_dist_data(grid,fid,'taij',taij)
        call write_data(grid,fid,'tajl',tajl)
      case (iowrite)            ! output to restart file
        call gather_zonal_trdiag
        call write_dist_data(grid,fid,'taijln',taijln)
        call write_dist_data(grid,fid,'taijls',taijls)
        call write_dist_data(grid,fid,'taijs',taijs)
        call write_dist_data(grid,fid,'taijn',taijn)
        call write_data(grid,fid,'tajln',tajln)
        call write_data(grid,fid,'tajls',tajls)
      case (ioread)            ! input from restart file
        call read_dist_data(grid,fid,'taijln',taijln)
        call read_dist_data(grid,fid,'taijls',taijls)
        call read_dist_data(grid,fid,'taijs',taijs)
        call read_dist_data(grid,fid,'taijn',taijn)
        call read_data(grid,fid,'tajln',tajln)
        call read_data(grid,fid,'tajls',tajls)
        call scatter_zonal_trdiag
      end select

      call new_io_tcons(fid,iaction)

      return
      end subroutine new_io_trdiag

      subroutine def_meta_trdiag(fid)
!@sum  def_meta_trdiag defines tracer metadata in acc files
!@auth M. Kelley
!@ver  beta
      use trdiag_com
      use pario, only : defvar,write_attr
      use domain_decomp_atm, only : grid
      use cdl_mod, only : defvar_cdl
      implicit none
      integer :: fid         !@var fid file id

      call write_attr(grid,fid,'taij','reduction','sum')
      call write_attr(grid,fid,'taij','split_dim',3)
      call defvar(grid,fid,ia_taij(1:ktaij_out),
     &     'ia_taij(ktaij)')
      call defvar(grid,fid,denom_taij(1:ktaij_out),
     &     'denom_taij(ktaij)')
      call defvar(grid,fid,scale_taij(1:ktaij_out),
     &     'scale_taij(ktaij)')
      call defvar(grid,fid,sname_taij(1:ktaij_out),
     &     'sname_taij(sname_strlen,ktaij)')
      call defvar(grid,fid,hemis_taij,'hemis_taij(one,shnhgm,ktaij)',
     &     r4_on_disk=.true.)
      call write_attr(grid,fid,'hemis_taij','reduction','sum')
      call defvar_cdl(grid,fid,cdl_taij,
     &     'cdl_taij(cdl_strlen,kcdl_taij)')
#ifdef CUBED_SPHERE
      call defvar_cdl(grid,fid,cdl_taij_latlon,
     &     'cdl_taij_latlon(cdl_strlen,kcdl_taij_latlon)')
#endif

      call write_attr(grid,fid,'taijl','reduction','sum')
      call write_attr(grid,fid,'taijl','split_dim',4)
      call defvar(grid,fid,ia_taijl(1:ktaijl_out),
     &     'ia_taijl(ktaijl)')
      call defvar(grid,fid,denom_taijl(1:ktaijl_out),
     &     'denom_taijl(ktaijl)')
      call defvar(grid,fid,scale_taijl(1:ktaijl_out),
     &     'scale_taijl(ktaijl)')
      call defvar(grid,fid,sname_taijl(1:ktaijl_out),
     &     'sname_taijl(sname_strlen,ktaijl)')
      call defvar_cdl(grid,fid,cdl_taijl,
     &     'cdl_taijl(cdl_strlen,kcdl_taijl)')
#ifdef CUBED_SPHERE
      call defvar_cdl(grid,fid,cdl_taijl_latlon,
     &     'cdl_taijl_latlon(cdl_strlen,kcdl_taijl_latlon)')
#endif

      call write_attr(grid,fid,'tajl','reduction','sum')
      call write_attr(grid,fid,'tajl','split_dim',3)
      call defvar(grid,fid,ia_tajl(1:ktajl_out),
     &     'ia_tajl(ktajl)')
      call defvar(grid,fid,denom_tajl(1:ktajl_out),
     &     'denom_tajl(ktajl)')
      call defvar(grid,fid,scale_tajl(1:ktajl_out),
     &     'scale_tajl(ktajl)')
      call defvar(grid,fid,sname_tajl(1:ktajl_out),
     &     'sname_tajl(sname_strlen,ktajl)')
      call defvar_cdl(grid,fid,cdl_tajl,
     &     'cdl_tajl(cdl_strlen,kcdl_tajl)')
      call defvar(grid,fid,hemis_tajl,'hemis_tajl(shnhgm,lm,ktajl)',
     &     r4_on_disk=.true.)
      call write_attr(grid,fid,'hemis_tajl','reduction','sum')
      call defvar(grid,fid,vmean_tajl,
     &     'vmean_tajl(jm_budg_plus3,one,ktajl)',r4_on_disk=.true.)
      call write_attr(grid,fid,'vmean_tajl','reduction','sum')

      call def_meta_tcons(fid)

      return
      end subroutine def_meta_trdiag

      subroutine write_meta_trdiag(fid)
!@sum  write_meta_trdiag write tracer accumulation metadata to file
!@auth M. Kelley
      use trdiag_com
      use pario, only : write_dist_data,write_data
      use domain_decomp_atm, only : grid
      use cdl_mod, only : write_cdl
      implicit none
      integer :: fid         !@var fid file id

      call write_data(grid,fid,'hemis_taij',hemis_taij)
      call write_data(grid,fid,'ia_taij',ia_taij(1:ktaij_out))
      call write_data(grid,fid,'denom_taij',denom_taij(1:ktaij_out))
      call write_data(grid,fid,'scale_taij',scale_taij(1:ktaij_out))
      call write_data(grid,fid,'sname_taij',sname_taij(1:ktaij_out))
      call write_cdl(grid,fid,'cdl_taij',cdl_taij)
#ifdef CUBED_SPHERE
      call write_cdl(grid,fid,'cdl_taij_latlon',cdl_taij_latlon)
#endif

      call write_data(grid,fid,'ia_taijl',ia_taijl(1:ktaijl_out))
      call write_data(grid,fid,'denom_taijl',denom_taijl(1:ktaijl_out))
      call write_data(grid,fid,'scale_taijl',scale_taijl(1:ktaijl_out))
      call write_data(grid,fid,'sname_taijl',sname_taijl(1:ktaijl_out))
      call write_cdl(grid,fid,'cdl_taijl',cdl_taijl)
#ifdef CUBED_SPHERE
      call write_cdl(grid,fid,'cdl_taijl_latlon',cdl_taijl_latlon)
#endif

      call write_data(grid,fid,'hemis_tajl',hemis_tajl)
      call write_data(grid,fid,'vmean_tajl',vmean_tajl)
      call write_data(grid,fid,'ia_tajl',ia_tajl(1:ktajl_out))
      call write_data(grid,fid,'denom_tajl',denom_tajl(1:ktajl_out))
      call write_data(grid,fid,'scale_tajl',scale_tajl(1:ktajl_out))
      call write_data(grid,fid,'sname_tajl',sname_tajl(1:ktajl_out))
      call write_cdl(grid,fid,'cdl_tajl',cdl_tajl)

      call write_meta_tcons(fid)

      return
      end subroutine write_meta_trdiag

#endif /* TRACERS_ON */

#if (defined TRACERS_ON) || (defined TRACERS_OCEAN)

      subroutine def_rsf_tcons(fid,r4_on_disk)
!@sum  def_rsf_tcons defines tracer diag array structure in restart+acc files
!@auth M. Kelley
!@ver  beta
      use trdiag_com, only :
     &     TCONSRV,
     &     TCONSRV_out
      use domain_decomp_atm, only : grid
      use pario, only : defvar
      implicit none
      integer fid           !@var fid file id
      logical :: r4_on_disk !@var r4_on_disk if true, real*8 stored as real*4
      if(r4_on_disk) then ! acc file
        call defvar(grid,fid,tconsrv_out,
     &       'tconsrv(jm_budg,ktcon)',r4_on_disk=.true.)
      else
        call defvar(grid,fid,tconsrv,
     &       'tconsrv(jm_budg,ktcon,ntmxcon)',r4_on_disk=r4_on_disk)
      endif
      return
      end subroutine def_rsf_tcons

      subroutine new_io_tcons(fid,iaction)
!@sum  new_io_tcons read/write tconsrv arrays from/to restart+acc files
!@auth M. Kelley
!@ver  beta new_ prefix avoids name clash with the default version
      use model_com, only : ioread,iowrite,iowrite_single
      use trdiag_com, only :
     &     TCONSRV,
     &     TCONSRV_out
      use domain_decomp_atm, only : grid
      use pario, only : write_dist_data,read_dist_data,
     &     write_data,read_data
      implicit none
      integer fid   !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
      select case (iaction)
      case (iowrite_single)     ! output to acc file
        call write_data(grid,fid,'tconsrv',tconsrv_out)
      case (iowrite)            ! output to restart file
        call gather_zonal_tcons
        call write_data(grid,fid,'tconsrv',tconsrv)
      case (ioread)            ! input from restart file
        call read_data(grid,fid,'tconsrv',tconsrv)
        call scatter_zonal_tcons
      end select
      return
      end subroutine new_io_tcons

      subroutine def_meta_tcons(fid)
!@sum  def_meta_trdiag defines tconsrv metadata in acc files
!@auth M. Kelley
!@ver  beta
      use trdiag_com
      use pario, only : defvar,write_attr
      use domain_decomp_atm, only : grid
      use cdl_mod, only : defvar_cdl
      implicit none
      integer :: fid         !@var fid file id

      call write_attr(grid,fid,'tconsrv','reduction','sum')
      call write_attr(grid,fid,'tconsrv','split_dim',2)
      call defvar(grid,fid,hemis_tconsrv,'hemis_tconsrv(shnhgm,ktcon)',
     &     r4_on_disk=.true.)
      call write_attr(grid,fid,'hemis_tconsrv','reduction','sum')
      call defvar(grid,fid,ia_tcon_out,'ia_tconsrv(ktcon)')
      call defvar(grid,fid,scale_tcon_out,'scale_tconsrv(ktcon)')
      call defvar(grid,fid,sname_tconsrv_out,
     &     'sname_tconsrv(sname_strlen,ktcon)')
      call defvar_cdl(grid,fid,cdl_tconsrv,
     &     'cdl_tconsrv(cdl_strlen,kcdl_tconsrv)')

      return
      end subroutine def_meta_tcons

      subroutine write_meta_tcons(fid)
!@sum  write_meta_tcons write tconsrv accumulation metadata to file
!@auth M. Kelley
      use trdiag_com
      use pario, only : write_dist_data,write_data
      use domain_decomp_atm, only : grid
      use cdl_mod, only : write_cdl
      implicit none
      integer :: fid         !@var fid file id

      call write_data(grid,fid,'hemis_tconsrv',hemis_tconsrv)
      call write_data(grid,fid,'ia_tconsrv',ia_tcon_out)
      call write_data(grid,fid,'scale_tconsrv',scale_tcon_out)
      call write_data(grid,fid,'sname_tconsrv',sname_tconsrv_out)
      call write_cdl(grid,fid,'cdl_tconsrv',cdl_tconsrv)

      return
      end subroutine write_meta_tcons

#endif /* TRACERS_ON or TRACERS_OCEAN */
#endif /* NEW_IO */


      SUBROUTINE ALLOC_TRDIAG_COM
      USE DIAG_COM, only : jm_budg
      USE TRDIAG_COM
      USE DOMAIN_DECOMP_ATM, only : GET, AM_I_ROOT, GRID
      use diag_zonal, only : get_alloc_bounds
      implicit none
      INTEGER :: J_0H,J_1H, I_0H,I_1H
      INTEGER :: status
      integer :: j_0budg,j_1budg
      integer :: img, jmg
      
      CALL GET(grid, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)
      I_0H=GRID%I_STRT_HALO
      I_1H=GRID%I_STOP_HALO

      call get_alloc_bounds(grid,
     &     j_strt_budg=j_0budg,j_stop_budg=j_1budg)

#ifdef TRACERS_WATER
      ALLOCATE ( TRP_acc(ntm,I_0H:I_1H,J_0H:J_1H),stat=status)
      ALLOCATE ( TRE_acc(ntm,I_0H:I_1H,J_0H:J_1H),stat=status)
#endif
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_DUST)
      ALLOCATE (  sPM2p5_acc(I_0H:I_1H,J_0H:J_1H),stat=status)
      ALLOCATE (   sPM10_acc(I_0H:I_1H,J_0H:J_1H),stat=status)
      ALLOCATE ( csPM2p5_acc(I_0H:I_1H,J_0H:J_1H),stat=status)
      ALLOCATE (  csPM10_acc(I_0H:I_1H,J_0H:J_1H),stat=status)
      ALLOCATE ( l1PM2p5_acc(I_0H:I_1H,J_0H:J_1H),stat=status)
      ALLOCATE (  l1PM10_acc(I_0H:I_1H,J_0H:J_1H),stat=status)

      sPM2p5_acc=0.d0; sPM10_acc=0.d0; l1PM2p5_acc=0.d0; l1PM10_acc=0.d0
      csPM2p5_acc=0.d0; csPM10_acc=0.d0
#endif
#ifdef TRACERS_ON 
      ALLOCATE(trcsurf(I_0H:I_1H,J_0H:J_1H,Ntm),stat=status)
      ALLOCATE(trcSurfByVol(I_0H:I_1H,J_0H:J_1H,Ntm),stat=status)
      ALLOCATE(trcSurfMixR_acc(I_0H:I_1H,J_0H:J_1H,Ntm),stat=status)
      ALLOCATE(trcSurfByVol_acc(I_0H:I_1H,J_0H:J_1H,Ntm),stat=status)
      ALLOCATE ( TAIJLN_loc(I_0H:I_1H,J_0H:J_1H,LM,ntm), stat=status )
      ALLOCATE ( TAIJLS_loc(I_0H:I_1H,J_0H:J_1H,LM,ktaijl), stat=status)
      ALLOCATE ( TAIJN_loc( I_0H:I_1H,J_0H:J_1H,ktaij,ntm),stat=status )
      ALLOCATE ( TAIJS_loc( I_0H:I_1H,J_0H:J_1H,ktaijs   ),stat=status )
      ALLOCATE ( TAJLN_loc(  J_0BUDG:J_1BUDG,LM,ktajlx,ntm),stat=status)
      ALLOCATE ( TAJLS_loc(  J_0BUDG:J_1BUDG,LM,ktajls    ),stat=status)

      if(am_i_root()) then
         img = IM
         jmg = JM
      else
         img = 1
         jmg = 1
      end if
      ALLOCATE ( TAIJLN(img,jmg,LM,ntm), stat=status )
      ALLOCATE ( TAIJLS(img,jmg,LM,ktaijl), stat=status )
      ALLOCATE ( TAIJN( img,jmg,ktaij,ntm), stat=status )
      ALLOCATE ( TAIJS( img,jmg,ktaijs   ), stat=status )

      ALLOCATE ( TAJLN(JM_BUDG,LM,ktajlx,ntm), stat=status )
      ALLOCATE ( TAJLS(JM_BUDG,LM,ktajls    ), stat=status )

#ifdef NEW_IO
! ktaij_, ktaijl_, ktajl_ are larger than necessary.  These arrays
! will be reallocated to the proper sizes later.
      ALLOCATE ( TAIJ_out( I_0H:I_1H,J_0H:J_1H,ktaij_),stat=status )
      ALLOCATE ( TAIJL_out( I_0H:I_1H,J_0H:J_1H,LM,ktaijl_),stat=status)
      if(am_i_root()) then
        ALLOCATE ( TAJL_out(JM_BUDG,LM,ktajl_), stat=status )
      endif
#endif

#endif  /* TRACERS_ON */
      RETURN
      END SUBROUTINE ALLOC_TRDIAG_COM

#if (defined TRACERS_ON) || (defined TRACERS_OCEAN)
C**** reset/gather/scatter routines for tconsrv special case
C**** can be called from both ocean tracer or atm tracer code
      subroutine reset_tcons
      USE TRDIAG_COM, only: TCONSRV_loc, TCONSRV
      USE DOMAIN_DECOMP_ATM, only : am_i_root
      implicit none
      TCONSRV_loc=0.
      if(am_i_root()) TCONSRV=0.
      return
      end subroutine reset_tcons

      subroutine gather_zonal_tcons
      USE DIAG_COM, only : ia_inst,jm_budg
      USE TRDIAG_COM, only: TCONSRV_loc, TCONSRV, ia_tcon,ntmxcon,ktcon
      USE DOMAIN_DECOMP_ATM, only : GRID,am_i_root,sumxpe
      implicit none
      real*8, dimension(:,:,:), allocatable :: tconsrv_sv
      integer :: n,k

c TCONSRV is a mixture of accumulations and instantaneous values, hence the
c complicated logic
! set inst qtys to zero in the "global" array before summing over PEs
      if(am_i_root()) then
        allocate(tconsrv_sv(jm_budg,ktcon,ntmxcon)); tconsrv_sv=tconsrv
        do n=1,ntmxcon
          do k=1,ktcon
            if(ia_tcon(k,n).eq.ia_inst) tconsrv(:,k,n)=0.
          enddo
        enddo
      endif

      call sumxpe (TCONSRV_loc,TCONSRV, increment=.true. )

      do n=1,ntmxcon ! keep instantaneous values
        do k=1,ktcon
          if(ia_tcon(k,n).ne.ia_inst) tconsrv_loc(:,k,n)=0. 
        enddo
      enddo

      if(am_i_root()) then
        do n=1,ntmxcon
          do k=1,ktcon
            if(ia_tcon(k,n).eq.ia_inst .and. all(tconsrv(:,k,n)==0.))
     &         tconsrv(:,k,n) = tconsrv_sv(:,k,n)
          enddo
        enddo
        deallocate(tconsrv_sv)
      endif

      return
      end subroutine gather_zonal_tcons

      subroutine scatter_zonal_tcons
      USE TRDIAG_COM, only: TCONSRV_loc, TCONSRV
      USE DOMAIN_DECOMP_ATM, only : GRID
      USE DIAG_ZONAL, only : unpack_lc
      implicit none
c      call unpack_lc   (grid, TCONSRV, TCONSRV_loc)
      tconsrv_loc = 0
      return
      end subroutine scatter_zonal_tcons
#endif

#ifdef TRACERS_ON
      subroutine reset_trdiag
      USE TRDIAG_COM, only: TAIJLN_loc, TAIJLS_loc, TAIJN_loc, 
     *     TAIJS_loc, TAJLN_loc, TAJLS_loc, TAJLN, TAJLS
      USE DOMAIN_DECOMP_ATM, only : am_i_root
      implicit none

       TAJLN_loc=0. ; TAJLS_loc=0. 
      TAIJLN_loc=0. ; TAIJLS_loc=0. ; TAIJN_loc=0. ; TAIJS_loc=0.
      if(am_i_root()) then
        TAJLN=0. ; TAJLS=0. 
      endif
      call reset_tcons

      return
      end subroutine reset_trdiag

      subroutine gather_trdiag
      USE TRDIAG_COM, only : TAIJLN, TAIJLN_loc, TAIJLS, TAIJLS_loc,
     *     TAIJN, TAIJN_loc,TAIJS, TAIJS_loc
      USE DOMAIN_DECOMP_1D, ONLY : GRID, PACK_DATA
      implicit none

      CALL PACK_DATA (GRID, TAIJLN_loc, TAIJLN)
      CALL PACK_DATA (GRID, TAIJLS_loc, TAIJLS)
      CALL PACK_DATA (GRID, TAIJN_loc , TAIJN)
      CALL PACK_DATA (GRID, TAIJS_loc , TAIJS)
      call gather_zonal_trdiag
      call gather_zonal_tcons

      return
      end subroutine gather_trdiag

      subroutine gather_zonal_trdiag
      USE TRDIAG_COM, only : TAJLN , TAJLN_loc, TAJLS, TAJLS_loc
      USE DOMAIN_DECOMP_ATM, ONLY : GRID
      USE DIAG_ZONAL, only : pack_lc
      implicit none
      call pack_lc   (grid, TAJLN_loc,  TAJLN )
      call pack_lc   (grid, TAJLS_loc,  TAJLS )
      return
      end subroutine gather_zonal_trdiag

      subroutine scatter_trdiag
      USE TRDIAG_COM, only : TAIJLN, TAIJLN_loc, TAIJLS, TAIJLS_loc,
     *     TAIJN, TAIJN_loc,TAIJS, TAIJS_loc
      USE DOMAIN_DECOMP_1D, ONLY : GRID, UNPACK_DATA
      implicit none
      CALL UNPACK_DATA (GRID, TAIJLN, TAIJLN_loc)
      CALL UNPACK_DATA (GRID, TAIJLS, TAIJLS_loc)
      CALL UNPACK_DATA (GRID, TAIJN , TAIJN_loc)
      CALL UNPACK_DATA (GRID, TAIJS , TAIJS_loc)
      call scatter_zonal_trdiag
      call scatter_zonal_tcons
      return
      end subroutine scatter_trdiag

      subroutine scatter_zonal_trdiag
      USE TRDIAG_COM, only : TAJLN , TAJLN_loc, TAJLS, TAJLS_loc
      USE DOMAIN_DECOMP_ATM, ONLY : GRID
      USE DIAG_ZONAL, only : unpack_lc
      implicit none
      call unpack_lc (grid, TAJLN,  TAJLN_loc )
      call unpack_lc (grid, TAJLS,  TAJLS_loc )
      return
      end subroutine scatter_zonal_trdiag

C**** routines for accumulating zonal mean diags (lat/lon grid)

      SUBROUTINE INC_TAJLS(I,J,L,TJL_INDEX,ACC)
!@sum inc_tajl adds ACC located at atmospheric gridpoint I,J,L
!@+   to the latitude-height zonal sum TAJLS(TJL_INDEX).
!@auth M. Kelley
      USE TRDIAG_COM, only : tajls=>tajls_loc
      USE DIAG_COM, only : wtbudg
      USE GEOM, only : j_budg
      IMPLICIT NONE
!@var I,J,L atm gridpoint indices for the accumulation
      INTEGER, INTENT(IN) :: I,J,L
!@var JL_INDEX index of the diagnostic being accumulated
      INTEGER, INTENT(IN) :: TJL_INDEX
!@var ACC increment of the diagnostic being accumulated
      REAL*8, INTENT(IN) :: ACC
C**** accumulate I,J value on the budget grid using j_budg to assign
C**** each point to a zonal mean (not bitwise reproducible for MPI).
      TAJLS(J_BUDG(I,J),L,TJL_INDEX) = TAJLS(J_BUDG(I,J),L,TJL_INDEX) +
     *     ACC!*wtbudg(I,J) ! cannot use wtbudg!=1 b/c kg units for tracers

      RETURN
      END SUBROUTINE INC_TAJLS

      SUBROUTINE INC_TAJLS2(I,J,L,TJL_INDEX,ACC)
!@sum inc_tajl adds ACC located at atmospheric gridpoint I,J,L
!@+   to the latitude-height zonal sum TAJLS(TJL_INDEX).
!@auth M. Kelley
      USE TRDIAG_COM, only : tajls=>tajls_loc
      USE DIAG_COM, only : wtbudg2
      USE GEOM, only : j_budg
      IMPLICIT NONE
!@var I,J,L atm gridpoint indices for the accumulation
      INTEGER, INTENT(IN) :: I,J,L
!@var JL_INDEX index of the diagnostic being accumulated
      INTEGER, INTENT(IN) :: TJL_INDEX
!@var ACC increment of the diagnostic being accumulated
      REAL*8, INTENT(IN) :: ACC
C**** accumulate I,J value on the budget grid using j_budg to assign
C**** each point to a zonal mean (not bitwise reproducible for MPI).
      TAJLS(J_BUDG(I,J),L,TJL_INDEX) = TAJLS(J_BUDG(I,J),L,TJL_INDEX) +
     *     ACC*wtbudg2(I,J)

      RETURN
      END SUBROUTINE INC_TAJLS2

      SUBROUTINE INC_TAJLS_COLUMN(I,J,L1,L2,NL,TJL_INDEX,ACC)
!@sum inc_tajl_column adds ACC(L1:L2) located at atmospheric gridpoint I,J
!@+   to the latitude-height zonal sums TAJLS(:,L1:L2,TJL_INDEX).
!@auth M. Kelley
      USE TRDIAG_COM, only : tajls=>tajls_loc
      USE DIAG_COM, only : wtbudg
      USE GEOM, only : j_budg
      IMPLICIT NONE
!@var I,J,L1,L2 atm gridpoint indices for the accumulation
      INTEGER, INTENT(IN) :: I,J,L1,L2,NL
!@var JL_INDEX index of the diagnostic being accumulated
      INTEGER, INTENT(IN) :: TJL_INDEX
!@var ACC increment of the diagnostic being accumulated
      REAL*8, INTENT(IN) :: ACC(NL)
      INTEGER :: JB,L
C**** accumulate I,J value on the budget grid using j_budg to assign
C**** each point to a zonal mean (not bitwise reproducible for MPI).
      JB = J_BUDG(I,J)
      DO L=L1,L2
        TAJLS(JB,L,TJL_INDEX) = TAJLS(JB,L,TJL_INDEX) +
     *     ACC(L)!*wtbudg(I,J) ! cannot use wtbudg!=1 b/c kg units for tracers
      ENDDO
      RETURN
      END SUBROUTINE INC_TAJLS_COLUMN

      SUBROUTINE INC_TAJLS2_COLUMN(I,J,L1,L2,NL,TJL_INDEX,ACC)
!@sum inc_tajl_column adds ACC(L1:L2) located at atmospheric gridpoint I,J
!@+   to the latitude-height zonal sums TAJLS(:,L1:L2,TJL_INDEX).
!@auth M. Kelley
      USE TRDIAG_COM, only : tajls=>tajls_loc
      USE DIAG_COM, only : wtbudg2
      USE GEOM, only : j_budg
      IMPLICIT NONE
!@var I,J,L1,L2 atm gridpoint indices for the accumulation
      INTEGER, INTENT(IN) :: I,J,L1,L2,NL
!@var JL_INDEX index of the diagnostic being accumulated
      INTEGER, INTENT(IN) :: TJL_INDEX
!@var ACC increment of the diagnostic being accumulated
      REAL*8, INTENT(IN) :: ACC(NL)
      INTEGER :: JB,L
C**** accumulate I,J value on the budget grid using j_budg to assign
C**** each point to a zonal mean (not bitwise reproducible for MPI).
      JB = J_BUDG(I,J)
      DO L=L1,L2
        TAJLS(JB,L,TJL_INDEX) = TAJLS(JB,L,TJL_INDEX) +
     *     ACC(L)*wtbudg2(I,J)
      ENDDO
      RETURN
      END SUBROUTINE INC_TAJLS2_COLUMN

      SUBROUTINE INC_TAJLN(I,J,L,TJL_INDEX,N,ACC)
!@sum inc_tajln adds ACC located at atmospheric gridpoint I,J,L
!@+   and tracer n to the latitude-height zonal sum TAJLN(TJL_INDEX,N).
!@auth M. Kelley
      USE TRDIAG_COM, only : tajln=>tajln_loc
      USE DIAG_COM, only : wtbudg
      USE GEOM, only : j_budg
      IMPLICIT NONE
!@var I,J,L atm gridpoint indices, N tracer no. for the accumulation
      INTEGER, INTENT(IN) :: I,J,L,N
!@var TJL_INDEX index of the diagnostic being accumulated
      INTEGER, INTENT(IN) :: TJL_INDEX
!@var ACC increment of the diagnostic being accumulated
      REAL*8, INTENT(IN) :: ACC
      TAJLN(J_BUDG(I,J),L,TJL_INDEX,N)=TAJLN(J_BUDG(I,J),L,TJL_INDEX,N)
     *     + ACC!*wtbudg(I,J) ! cannot use wtbudg!=1 b/c kg units for tracers

      RETURN
      END SUBROUTINE INC_TAJLN

      SUBROUTINE INC_TAJLN_COLUMN(I,J,L1,L2,NL,TJL_INDEX,N,ACC)
!@sum inc_tajln adds ACC located at atmospheric gridpoint I,J,L
!@+   and tracer n to the latitude-height zonal sum TAJLN(TJL_INDEX,N).
!@auth M. Kelley
      USE TRDIAG_COM, only : tajln=>tajln_loc
      USE DIAG_COM, only : wtbudg
      USE GEOM, only : j_budg
      IMPLICIT NONE
!@var I,J,L atm gridpoint indices, N tracer no. for the accumulation
      INTEGER, INTENT(IN) :: I,J,L1,L2,NL,N
!@var TJL_INDEX index of the diagnostic being accumulated
      INTEGER, INTENT(IN) :: TJL_INDEX
!@var ACC increment of the diagnostic being accumulated
      REAL*8, INTENT(IN) :: ACC(NL)
      INTEGER :: JB,L
      JB = J_BUDG(I,J)
      DO L=L1,L2
        TAJLN(JB,L,TJL_INDEX,N)=TAJLN(JB,L,TJL_INDEX,N) +
     *    ACC(L) !*wtbudg(I,J) ! cannot use wtbudg!=1 b/c kg units for tracers
      ENDDO
      RETURN
      END SUBROUTINE INC_TAJLN_COLUMN
#endif   /* TRACERS_ON */
