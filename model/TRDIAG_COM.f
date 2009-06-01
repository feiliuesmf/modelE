#include "rundeck_opts.h"

      MODULE TRDIAG_COM
!@sum Tracer diagnostic arrays
!@+    Mostly tracer independent, but this may depend on applications
!@auth Jean Lerner
!ver   1.0
      USE MODEL_COM, only: im,jm,lm
      USE DIAG_COM, only: npts !npts are conservation quantities
     &     ,sname_strlen,units_strlen,lname_strlen
#if (defined TRACERS_ON) || (defined TRACERS_OCEAN)
      USE TRACER_COM, only: ntm
     *     , ntsurfsrcmax, nt3Dsrcmax
#ifdef TRACERS_AMP
     *     ,ntmAMP
      USE AERO_CONFIG,only: nbins
#endif
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
#endif
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
      !real*4, dimension(IM,JM,LM,ntm)         :: taijln4
      !real*8, allocatable, dimension(:,:,:,:) :: taijln4_loc
!@var SNAME_IJT, UNITS_IJT: Names and units of lat-sigma tracer IJ diags
      character(len=sname_strlen), dimension(lm,ntm) :: sname_ijt
      character(len=units_strlen), dimension(lm,ntm) :: units_ijt
!@var LNAME_IJT: descriptions of tracer IJ diagnostics
      character(len=lname_strlen), dimension(lm,ntm) ::
     &     lname_ijt = 'unused'
!@var SCALE_IJT: printout scaling factor for tracer IJ diagnostics
      REAL*8, dimension(lm,ntm) :: scale_ijt
!@var IR_IJT: range index of IJ diagnostics
      integer, dimension(ntm) :: ir_ijt
!@var IA_IJT: idacc-number for tracer IJ diags (same for all tracers)
      integer ia_ijt
!@var IJTC_POWER: power of 10 used for tracer IJ concentration diags
      integer, dimension(ntm) :: ijtc_power
!@var IJTM_POWER: power of 10 used for tracer IJ mass unit diags
      integer, dimension(ntm) :: ijtm_power
!@var IJT_XX Names for TAIJLN diagnostics (Not currently used)

C**** TAIJN
!@param KTAIJ number of 2D diags describing surface and column load along with wet and dry deposition
!@+   please just increase this if needed - don't bother with pp options
      integer, parameter :: ktaij=18

!@var IJT_XX names for taijn diagnostics
      integer tij_conc,tij_surf,tij_surfbv,tij_mass
!@var IJT_XX names for water-based taijn diagnostics
      integer tij_rvr,tij_seaice,tij_prec,tij_evap,tij_grnd,tij_lk1
     *     ,tij_lk2,tij_soil,tij_snow,tij_uflx,tij_vflx,tij_icocflx   
      integer tij_kw,tij_alpha  ! TRACERS_GASEXCH_ocean
      integer tij_drydep,tij_gsdep ! TRACERS_DRYDEP

!@var TAIJN lat/lon tracer diagnostics (all tracers)
      real*8, allocatable, dimension(:,:,:,:)      :: taijn
      real*8, allocatable, dimension(:,:,:,:) :: taijn_loc
      !real*4, dimension(IM,JM,ktaij,ntm)      :: taijn4
      !real*8, allocatable, dimension(:,:,:,:) :: taijn4_loc
!@var SCALE_TIJ: printout scaling factor for tracer IJK diagnostics
      REAL*8, dimension(ktaij,ntm) :: scale_tij
!@var SNAME_TIJ,UNITS_TIJ: Names and units of lat-sigma tracer diags
      character(len=sname_strlen), dimension(ktaij,ntm) :: sname_tij
      character(len=units_strlen), dimension(ktaij,ntm) :: units_tij
!@var LNAME_TIJ: descriptions of tracer IJK diags
      character(len=lname_strlen), dimension(ktaij,ntm) ::
     &     lname_tij = 'unused'

C**** TAIJS  <<<< KTAIJS and IJTS_xx are Tracer-Dependent >>>>
!@var ijs_XXX index for diags not specific to a certain tracer
      INTEGER :: ijs_ai,ijs_isoprene,ijs_NO2_col,ijs_NO2_count
      INTEGER, DIMENSION(LM) :: ijs_OH,ijs_NO3,ijs_HO2,ijs_JH2O2
#ifdef HTAP_LIKE_DIAGS
     & ,ijs_COp,ijs_COd,ijs_Oxp,ijs_Oxd,ijs_CH4d
#endif

!@param KTAIJS number of special lat/lon tracer diagnostics
!@+   please just increase this if needed - don't bother with pp options
      INTEGER,PARAMETER :: ktaijs=2223

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
!@param nDustEmij index of dust emission in ijts_isrc
      INTEGER,PARAMETER :: nDustEmij=1
!@param nDustTurbij index of dust dry turbulent deposition in ijts_isrc
      INTEGER,PARAMETER :: nDustTurbij=3  ! not used?
!@param nDustEv1ij index of number of dust events below threshold wind
!@param nDustEv1ij in ijts_spec
!@param nDustEv2ij index of number of dust events above threshold wind
!@param nDustEv2ij in ijts_spec
!@param nDustWthij index of threshold velocity in ijts_spec
      INTEGER,PARAMETER :: nDustEv1ij=1,nDustEv2ij=2,nDustWthij=3
#endif
#ifdef TRACERS_DUST
!@param nDustEm2ij index of dust emission according to cubic scheme
!@param nDustEm2ij in ijts_isrc
      INTEGER,PARAMETER :: nDustEm2ij=2
#endif

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
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:)       :: TAIJS
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: TAIJS_loc
      !REAL*4, DIMENSION(IM,JM,ktaijs)       :: TAIJS4
      !REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: TAIJS4_loc
!@var ijts_source tracer independent array for TAIJS surface src. diags
      INTEGER ijts_source(ntsurfsrcmax,ntm)
!@var ijts_isrc tracer independent array for TAIJS interactive srf. src.
      INTEGER ijts_isrc(ntsurfsrcmax,ntm)
!@var ijts_aq tracer independent array for TAIJS aqueous change
      INTEGER ijts_aq(ntm)
!@var ijts_alb BC impact on snow/ice albedo, grain size,sw and lw radiation
      INTEGER ijts_alb(2,ntm)
!@var ijts_tau tracer independent array for TAIJS hydrated opt. thick.
      INTEGER ijts_tau(2,ntm)
!@var ijts_3Dtau 3D tracer independent array for TAIJS hydrated opt. thick.
      INTEGER ijts_3Dtau(lm,ntm)
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
      INTEGER ijts_fc(6,ntm)
!@var ijts_fcsub index for TAIJS SW/LW rad. forc. for tracer sub classes
      INTEGER ijts_fcsub(6,Ntm,MaxSubCl)
!@var ijts_spec index for TAIJS for special diags. not associated with single
!@var ijts_spec tracer
      INTEGER :: ijts_spec(MaxSpec)
#ifdef TRACERS_AMP
!@var ijts_AMPpdf special diagnostic for not-transported tracers
      INTEGER ijts_AMPpdf(1,nbins)
#endif
!@var ijts_AMPext special diagnostic for not-transported tracers
      INTEGER ijts_AMPext(lm,6)
!@var ijts_AMPm tracer independent array for AMP modes
      INTEGER ijts_AMPm(lm,2,Ntm)
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
!@var SCALE_IJTS: printout scaling factor for tracer IJTS diagnostics
      REAL*8, dimension(ktaijs) :: scale_ijts
!@var IA_IJTS: idacc-number for tracer source/sink IJ diags
      integer ia_ijts(ktaijs)
!@var ijts_power: power of 10 used for tracer IJ source/sink diags
      INTEGER, DIMENSION(ktaijs) :: ijts_power
!@var ijts_index: tracer index associated with a TAIJS diagnostic
      INTEGER, DIMENSION(ktaijs) :: ijts_index

#if (defined TRACERS_DUST) && (defined TRACERS_DRYDEP)
!@var rts_save saves rts as global field for tracer diagnostics
      REAL*8 :: rts_save(Im,Jm)
#endif

C**** TAJLN
!@param ktajl,ktajlx number of TAJL tracer diagnostics;
!@+          ktajlx includes composites
      INTEGER, PARAMETER :: ktajl=10
     &                     ,ktajlx=ktajl+2
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
      integer ktajlp !@var helps keeps track of
!@var scale_jln: Scale for jl maps
      REAL*8, DIMENSION(ntm) :: scale_jln

C**** TAJLS  <<<< KTAJLS and JLS_xx are Tracer-Dependent >>>>
!@param ktajls number of source/sink TAJLS tracer diagnostics;
!@+   please just increase this if needed - don't bother with pp options
      INTEGER,PARAMETER :: ktajls=400

!@var jls_XXX index for non-tracer specific or special diags
      INTEGER jls_OHconk,jls_HO2con,jls_NO3
     *     ,jls_phot,jls_incloud(2,ntm), jls_OHcon,jls_H2Omr
     *     ,jls_N2O5sulf,jls_day,jls_COd,jls_COp,jls_Oxd,jls_Oxp
     *     ,jls_ClOcon,jls_H2Ocon,jls_H2Ochem,jls_OxdT,jls_OxpT

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
!@param nDustEmjl index of dust emission in jls_source
      INTEGER,PARAMETER :: nDustEmjl=1
!@param nDustTurbjl index of dust dry turbulent deposition in jls_source
      INTEGER,PARAMETER :: nDustTurbjl=3
!@param nDustEv1jl index of number of dust events below threshold wind
!@param nDustEv1jl in jls_spec
!@param nDustEv2jl index of number of dust events above threshold wind
!@param nDustEv2jl in jls_spec
!@param nDustWthjl index of threshold velocity in ijts_spec
      INTEGER,PARAMETER :: nDustEv1jl=1,nDustEv2jl=2,nDustWthjl=3
#endif
#ifdef TRACERS_DUST
!@param nDustEm2jl index of dust emission according to cubic scheme
!@param nDustEm2jl in jls_source
      INTEGER,PARAMETER :: nDustEm2jl=2
#endif
!@var TAJLS  JL special tracer diagnostics for sources, sinks, etc
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: TAJLS
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: TAJLS_loc
!@var jls_source tracer independent array for TAJLS surface src. diags
      INTEGER jls_source(ntsurfsrcmax,ntm)
!@var jls_isrc tracer independent array for TAJLS interactive surface src. diags
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

C**** TCONSRV
!@param NTCONS Maximum Number of special tracer conservation points
      INTEGER, PARAMETER :: ntcons=20
!@param KTCON total number of conservation diagnostics for tracers
      INTEGER, PARAMETER :: KTCON=npts+ntcons+2
!@param ntmxcon total number of conservation quantities
#ifdef TRACERS_SPECIAL_Shindell
C**** include some extra troposphere only ones
      INTEGER, PARAMETER :: ntmxcon = ntm+3
#else
      INTEGER, PARAMETER :: ntmxcon = ntm
#endif
!@var TCONSRV conservation diagnostics for tracers
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: TCONSRV
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: TCONSRV_loc
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
     &     name_tconsrv
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
      INTEGER,DIMENSION(ntmxcon) :: itcon_wt
#endif
!@var PDSIGJL temporary storage for mean pressures for jl diags
      REAL*8, DIMENSION(JM,LM)            :: PDSIGJL

#ifdef TRACERS_WATER
!@var TRP_acc, TRE_acc accumulation arrays for some SUBDD diags
!!    REAL*8 TRP_acc(ntm,IM,JM), TRE_acc(ntm,IM,JM)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:), public :: TRP_acc,TRE_acc
#endif

#ifdef NEW_IO
c declarations that facilitate summation of acc-files when using
c the new i/o system
      target :: 
     &     TAIJLN_loc,TAIJN_loc,TAIJS_loc,
     &     TAJLN_loc,TAJLS_loc,TCONSRV_loc
      real*8, allocatable, dimension(:,:,:,:), target ::
     &     taijln_fromdisk,taijn_fromdisk,tajln_fromdisk
      real*8, allocatable, dimension(:,:,:), target ::
     &     taijs_fromdisk,tajls_fromdisk,tconsrv_fromdisk
      real*8, dimension(:,:,:,:), pointer ::
     &     taijln_ioptr,taijn_ioptr,tajln_ioptr
      real*8, dimension(:,:,:), pointer ::
     &     taijs_ioptr,tajls_ioptr,tconsrv_ioptr
#endif

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

C**** remove spaces in NAME_CON for netcdf names
      sname=TRIM(NAME_CON)
      do k=1,len_trim(NAME_CON)
        if (sname(k:k).eq." ") sname(k:k)="_"
      end do

C**** remove spaces, invalid characters in CONPTS, CONPT0 for netcdf names
      conpts_sname = conpts
      do n=1,ntcons
      do k=1,len_trim(conpts_sname(n))
         if (conpts_sname(n)(k:k).eq." ") conpts_sname(n)(k:k)="_"
         if (conpts_sname(n)(k:k).eq."+") conpts_sname(n)(k:k)="_"
      enddo
      enddo
      conpt0_sname = conpt0
      do n=1,npts
      do k=1,len_trim(conpt0_sname(n))
         if (conpt0_sname(n)(k:k).eq." ") conpt0_sname(n)(k:k)="_"
         if (conpt0_sname(n)(k:k).eq."+") conpt0_sname(n)(k:k)="_"
      enddo
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

#endif
#ifdef TRACERS_ON
      SUBROUTINE io_trdiag(kunit,it,iaction,ioerr)
!@sum  io_trdiag reads and writes tracer diagnostics arrays to file
!@auth Jean Lerner
!@ver  1.0
      USE MODEL_COM, only: im,jm,lm
      USE MODEL_COM, only: ioread,iowrite,iowrite_mon,iowrite_single
     *     ,irerun,ioread_single,lhead
      USE DIAG_COM, only : jm_budg
      USE DOMAIN_DECOMP_1D, only : AM_I_ROOT
      USE DOMAIN_DECOMP_1D, only : ESMF_BCAST
      USE DOMAIN_DECOMP_1D, only : grid, GET
      USE TRACER_COM, only: ntm

      USE TRDIAG_COM, only : taijln_loc, taijln
      USE TRDIAG_COM, only : taijn_loc,  taijn
      USE TRDIAG_COM, only : taijs_loc,  taijs
      USE TRDIAG_COM, only : tajln_loc,  tajln
      USE TRDIAG_COM, only : tajls_loc,  tajls
      USE TRDIAG_COM, only : tconsrv_loc,tconsrv

      USE TRDIAG_COM, only : pdsigjl

      USE TRDIAG_COM, only : ktaij, ktaijs, ktajlx, ktajls, ktcon
      USE TRDIAG_COM, only : ntmxcon

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
     *  ktacc=IM*JM*LM*NTM + IM*JM*ktaij*NTM + IM*JM*ktaijs +
     *     JM_BUDG*LM*ktajlx*NTM + JM_BUDG*LM*ktajls + JM_BUDG*ntmxcon
     *     *ktcon
!@var TA..4(...) dummy arrays for reading diagnostics files
      real*4, allocatable, dimension(:,:,:,:) :: taijln4
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
     *                         TAIJLN,TAIJN,TAIJS,
     *                         TAJLN ,TAJLS,TCONSRV,it
          CASE (IOWRITE_SINGLE)    ! output to acc file
            MODULE_HEADER(LHEAD+1:LHEAD+2) = 'R4'
            IF (AM_I_ROOT()) WRITE (kunit,err=10) MODULE_HEADER,
     *       REAL(TAIJLN,KIND=4),REAL(TAIJN,KIND=4),
     *       REAL(TAIJS,KIND=4) ,REAL(TAJLN,KIND=4),
     *       REAL(TAJLS,KIND=4) ,REAL(TCONSRV,KIND=4),it
        END SELECT
      CASE (IOREAD:)          ! input from restart file
        SELECT CASE (IACTION)
        CASE (ioread_single)    ! accumulate diagnostic files
          if(am_i_root()) then
            ALLOCATE (taijln4(IM,JM,LM,ntm), stat=status )
            ALLOCATE (taijn4(IM,JM,ktaij,ntm), stat=status )
            ALLOCATE (TAIJS4(IM,JM,ktaijs), stat=status )
            ALLOCATE (TAJLN4(JM_BUDG,LM,ktajlx,ntm), stat=status )
            ALLOCATE (TAJLS4(JM_BUDG,LM,ktajls), stat=status )
            ALLOCATE (TCONSRV4(JM_BUDG,ktcon,ntmxcon), stat=status )

            READ (kunit,err=10) HEADER,
     *         TAIJLN4,TAIJN4,TAIJS4,TAJLN4,TAJLS4,TCONSRV4,it_check
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
            TAIJN  = TAIJN  + TAIJN4
            TAIJS  = TAIJS  + TAIJS4
            TAJLN  = TAJLN  + TAJLN4
            TAJLS  = TAJLS  + TAJLS4
            TCONSRV= TCONSRV+ TCONSRV4
            
            DEALLOCATE ( taijln4 )
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
     *                           TAIJLN,TAIJN,TAIJS,
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
     &     TAIJLN=>TAIJLN_ioptr,
     &     TAIJN=>TAIJN_ioptr,
     &     TAIJS=>TAIJS_ioptr,
     &     TAJLN=>TAJLN_ioptr,
     &     TAJLS=>TAJLS_ioptr,
     &     TCONSRV=>TCONSRV_ioptr
      use domain_decomp_atm, only : grid
      use pario, only : defvar
      implicit none
      integer fid           !@var fid file id
      logical :: r4_on_disk !@var r4_on_disk if true, real*8 stored as real*4
      call defvar(grid,fid,taijln,
     &     'taijln(dist_im,dist_jm,lm,ntm)',r4_on_disk=r4_on_disk)
      call defvar(grid,fid,taijs,
     &     'taijs(dist_im,dist_jm,ktaijs)',r4_on_disk=r4_on_disk)
      call defvar(grid,fid,taijn,
     &     'taijn(dist_im,dist_jm,ktaij,ntm)',r4_on_disk=r4_on_disk)
      call defvar(grid,fid,tajln,
     &     'tajln(jm_budg,lm,ktajlx,ntm)',r4_on_disk=r4_on_disk)
      call defvar(grid,fid,tajls,
     &     'tajls(jm_budg,lm,ktajls)',r4_on_disk=r4_on_disk)
      call defvar(grid,fid,tconsrv,
     &     'tconsrv(jm_budg,ktcon,ntmxcon)',r4_on_disk=r4_on_disk)
      return
      end subroutine def_rsf_trdiag

      subroutine new_io_trdiag(fid,iaction)
!@sum  new_io_trdiag read/write tracer acc arrays from/to restart+acc files
!@auth M. Kelley
!@ver  beta new_ prefix avoids name clash with the default version
      use model_com, only : ioread,iowrite
c in the postprocessing case where arrays are read from disk and summed,
c these i/o pointers point to temporary arrays.  Otherwise, they point to
c the instances of the arrays used during normal operation. 
      use trdiag_com, only :
     &     TAIJLN=>TAIJLN_ioptr,
     &     TAIJN=>TAIJN_ioptr,
     &     TAIJS=>TAIJS_ioptr,
     &     TAJLN=>TAJLN_ioptr,
     &     TAJLS=>TAJLS_ioptr,
     &     TCONSRV=>TCONSRV_ioptr
      use domain_decomp_atm, only : grid
      use pario, only : write_dist_data,read_dist_data
      implicit none
      integer fid   !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
      select case (iaction)
      case (iowrite)            ! output to restart or acc file
        call gather_zonal_trdiag
        call write_dist_data(grid,fid,'taijln',taijln)
        call write_dist_data(grid,fid,'taijs',taijs)
        call write_dist_data(grid,fid,'taijn',taijn)
        call write_data(grid,fid,'tajln',tajln)
        call write_data(grid,fid,'tajls',tajls)
        call write_data(grid,fid,'tconsrv',tconsrv)
      case (ioread)            ! input from restart or acc file
        call read_dist_data(grid,fid,'taijln',taijln)
        call read_dist_data(grid,fid,'taijs',taijs)
        call read_dist_data(grid,fid,'taijn',taijn)
        call read_data(grid,fid,'tajln',tajln)
        call read_data(grid,fid,'tajls',tajls)
        call read_data(grid,fid,'tconsrv',tconsrv)
        call scatter_zonal_trdiag
      end select
      return
      end subroutine new_io_trdiag

      subroutine set_ioptrs_tracacc_default
c point i/o pointers for diagnostic accumlations to the
c instances of the arrays used during normal operation. 
      use trdiag_com
      implicit none
      taijln_ioptr  => taijln_loc
      taijs_ioptr   => taijs_loc
      taijn_ioptr   => taijn_loc
      tajls_ioptr   => tajls!_loc
      tajln_ioptr   => tajln!_loc
      tconsrv_ioptr => tconsrv!_loc
      return
      end subroutine set_ioptrs_tracacc_default

      subroutine set_ioptrs_tracacc_sumfiles
c point i/o pointers for diagnostic accumlations to temporary
c arrays that hold data read from disk
      use trdiag_com
      use domain_decomp_atm, only : grid,am_i_root
      implicit none
      integer :: i_0h,i_1h,j_0h,j_1h
c
c allocate temporary arrays (local size)
c
      if(allocated(taijln_fromdisk)) return
      I_0H = grid%I_STRT_HALO
      I_1H = grid%I_STOP_HALO
      J_0H = grid%J_STRT_HALO
      J_1H = grid%J_STOP_HALO

      allocate(taijln_fromdisk(i_0h:i_1h,j_0h:j_1h,lm,ntm))
      taijln_ioptr => taijln_fromdisk
      allocate(taijs_fromdisk(i_0h:i_1h,j_0h:j_1h,ktaijs))
      taijs_ioptr => taijs_fromdisk
      allocate(taijn_fromdisk(i_0h:i_1h,j_0h:j_1h,ktaij,ntm))
      taijn_ioptr => taijn_fromdisk
      if(am_i_root()) then
        allocate(tajls_fromdisk(jm_budg,lm,ktajls))
        tajls_ioptr => tajls_fromdisk
        allocate(tajln_fromdisk(jm_budg,lm,ktajlx,ntm))
        tajln_ioptr => tajln_fromdisk
        allocate(tconsrv_fromdisk(jm_budg,ktcon,ntmxcon))
        tconsrv_ioptr => tconsrv_fromdisk
      endif
      return
      end subroutine set_ioptrs_tracacc_sumfiles

      subroutine sumfiles_tracacc
c increment diagnostic accumlations with the data that was
c read from disk and stored in the _fromdisk arrays.
      use trdiag_com
      use domain_decomp_atm, only : am_i_root
      implicit none
      taijln_loc  = taijln_loc  + taijln_fromdisk
      taijs_loc   = taijs_loc   + taijs_fromdisk
      taijn_loc   = taijn_loc   + taijn_fromdisk
      if(am_i_root()) then
        tajls   = tajls   + tajls_fromdisk
        tajln   = tajln   + tajln_fromdisk
        tconsrv = tconsrv + tconsrv_fromdisk
      endif
      return
      end subroutine sumfiles_tracacc
#endif /* TRACERS_ON */
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
      
      CALL GET(grid, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)
      I_0H=GRID%I_STRT_HALO
      I_1H=GRID%I_STOP_HALO

      call get_alloc_bounds(grid,
     &     j_strt_budg=j_0budg,j_stop_budg=j_1budg)

#ifdef TRACERS_WATER
      ALLOCATE ( TRP_acc(ntm,I_0H:I_1H,J_0H:J_1H),stat=status)
      ALLOCATE ( TRE_acc(ntm,I_0H:I_1H,J_0H:J_1H),stat=status)
#endif
#ifdef TRACERS_ON 
      ALLOCATE ( TAIJLN_loc(I_0H:I_1H,J_0H:J_1H,LM,ntm), stat=status )
      ALLOCATE ( TAIJN_loc( I_0H:I_1H,J_0H:J_1H,ktaij,ntm),stat=status )
      ALLOCATE ( TAIJS_loc( I_0H:I_1H,J_0H:J_1H,ktaijs   ),stat=status )
      ALLOCATE ( TAJLN_loc(  J_0BUDG:J_1BUDG,LM,ktajlx,ntm),stat=status)
      ALLOCATE ( TAJLS_loc(  J_0BUDG:J_1BUDG,LM,ktajls    ),stat=status)
      ALLOCATE ( TCONSRV_loc(J_0BUDG:J_1BUDG,ktcon,ntmxcon),stat=status)

      if(am_i_root()) then
        ALLOCATE ( TAIJLN(IM,JM,LM,ntm), stat=status )
        ALLOCATE ( TAIJN( IM,JM,ktaij,ntm), stat=status )
        ALLOCATE ( TAIJS( IM,JM,ktaijs   ), stat=status )
        ALLOCATE ( TAJLN(JM_BUDG,LM,ktajlx,ntm), stat=status )
        ALLOCATE ( TAJLS(JM_BUDG,LM,ktajls    ), stat=status )
        ALLOCATE ( TCONSRV(JM_BUDG,ktcon,ntmxcon), stat=status )
      endif
#endif
      RETURN
      END SUBROUTINE ALLOC_TRDIAG_COM

#ifdef TRACERS_ON
      subroutine reset_trdiag
      USE TRDIAG_COM, only: TAIJLN_loc, TAIJN_loc, 
     *     TAIJS_loc, TAJLN_loc, TAJLS_loc, TCONSRV_loc,
     *     TAJLN, TAJLS, TCONSRV
      USE DOMAIN_DECOMP_ATM, only : am_i_root
      implicit none

       TAJLN_loc=0. ; TAJLS_loc=0. ; TCONSRV_loc=0.
      TAIJLN_loc=0. ; TAIJN_loc=0. ; TAIJS_loc=0.
      if(am_i_root()) then
        TAJLN=0. ; TAJLS=0. ; TCONSRV=0.
      endif

      return
      end subroutine reset_trdiag

      subroutine gather_trdiag
      USE TRDIAG_COM, only : TAIJLN, TAIJLN_loc, TAIJN, TAIJN_loc,
     *     TAIJS, TAIJS_loc
      USE DOMAIN_DECOMP_1D, ONLY : GRID, PACK_DATA
      implicit none

      CALL PACK_DATA (GRID, TAIJLN_loc, TAIJLN)
      CALL PACK_DATA (GRID, TAIJN_loc , TAIJN)
      CALL PACK_DATA (GRID, TAIJS_loc , TAIJS)
      call gather_zonal_trdiag

      return
      end subroutine gather_trdiag

      subroutine gather_zonal_trdiag
      USE TRDIAG_COM, only : TAJLN , TAJLN_loc, TAJLS, TAJLS_loc,
     *     TCONSRV, TCONSRV_loc
      USE DOMAIN_DECOMP_ATM, ONLY : GRID
      USE DIAG_ZONAL, only : pack_lc
      implicit none
      call pack_lc   (grid, TAJLN_loc,  TAJLN )
      call pack_lc   (grid, TAJLS_loc,  TAJLS )
      call pack_lc   (grid, TCONSRV_loc,TCONSRV )
      return
      end subroutine gather_zonal_trdiag

      subroutine scatter_trdiag
      USE TRDIAG_COM, only : TAIJLN, TAIJLN_loc, TAIJN, TAIJN_loc,
     *     TAIJS, TAIJS_loc
      USE DOMAIN_DECOMP_1D, ONLY : GRID, UNPACK_DATA
      implicit none
      CALL UNPACK_DATA (GRID, TAIJLN, TAIJLN_loc)
      CALL UNPACK_DATA (GRID, TAIJN , TAIJN_loc)
      CALL UNPACK_DATA (GRID, TAIJS , TAIJS_loc)
      call scatter_zonal_trdiag
      return
      end subroutine scatter_trdiag

      subroutine scatter_zonal_trdiag
      USE TRDIAG_COM, only : TAJLN , TAJLN_loc, TAJLS, TAJLS_loc,
     *     TCONSRV, TCONSRV_loc
      USE DOMAIN_DECOMP_ATM, ONLY : GRID
      USE DIAG_ZONAL, only : unpack_lc
      implicit none
      call unpack_lc (grid, TAJLN,  TAJLN_loc )
      call unpack_lc (grid, TAJLS,  TAJLS_loc )
      call unpack_lc (grid, TCONSRV,TCONSRV_loc )
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

      SUBROUTINE INC_TAJLN(I,J,L,TJL_INDEX,N,ACC)
!@sum inc_tajln adds ACC located at atmospheric gridpoint I,J,L
!@+   and tracer n to the latitude-height zonal sum TAJLN(TJL_INDEX,N).
!@auth M. Kelley
      USE TRDIAG_COM, only : tajln=>tajln_loc
      USE DIAG_COM, only : wtbudg
      USE GEOM, only : j_budg
      IMPLICIT NONE
!@var I,J,L atm gridpoint indices, N tracer # for the accumulation
      INTEGER, INTENT(IN) :: I,J,L,N
!@var TJL_INDEX index of the diagnostic being accumulated
      INTEGER, INTENT(IN) :: TJL_INDEX
!@var ACC increment of the diagnostic being accumulated
      REAL*8, INTENT(IN) :: ACC
      TAJLN(J_BUDG(I,J),L,TJL_INDEX,N)=TAJLN(J_BUDG(I,J),L,TJL_INDEX,N)
     *     + ACC!*wtbudg(I,J) ! cannot use wtbudg!=1 b/c kg units for tracers

      RETURN
      END SUBROUTINE INC_TAJLN
#endif
