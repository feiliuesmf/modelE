#include "rundeck_opts.h"
      

!@sum  TRACER_COM: Exists alone to minimize the number of dependencies
!@+    This version for simple trace gases/chemistry/isotopes
!@auth Jean Lerner/Gavin Schmidt

      MODULE TRACER_COM
!@sum  TRACER_COM tracer variables
!@auth Jean Lerner
C
      USE QUSDEF, only: nmom
      USE RESOLUTION, only: im,jm,lm
      use OldTracer_mod, only: trName
      use OldTracer_mod, only: tr_mm
      use OldTracer_mod, only: ntm_power
      use OldTracer_mod, only: t_qlimit
      use OldTracer_mod, only: needtrs
      use OldTracer_mod, only: trdecay

      use OldTracer_mod, only: itime_tr0
      use OldTracer_mod, only: mass2vol
      use OldTracer_mod, only: vol2mass

      use OldTracer_mod, only: dodrydep
      use OldTracer_mod, only: F0
      use OldTracer_mod, only: HSTAR

      use OldTracer_mod, only: do_fire
      use OldTracer_mod, only: nBBsources
      use OldTracer_mod, only: emisPerFireByVegType
      use OldTracer_mod, only: trpdens
      use OldTracer_mod, only: trradius

      use OldTracer_mod, only: tr_wd_TYPE
      use OldTracer_mod, only: tr_RKD
      use OldTracer_mod, only: tr_DHD
      use OldTracer_mod, only: fq_aer
      use OldTracer_mod, only: rc_washt
      use OldTracer_mod, only: isDust

      use OldTracer_mod, only: nGAS, nPART, nWATER

      use OldTracer_mod, only: tr_H2ObyCH4
      use OldTracer_mod, only: dowetdep
      use OldTracer_mod, only: trw0
      use OldTracer_mod, only: ntrocn
      use OldTracer_mod, only: conc_from_fw
      use OldTracer_mod, only: trglac
      use OldTracer_mod, only: ntisurfsrc

      use OldTracer_mod, only: trli0
      use OldTracer_mod, only: trsi0
#ifdef TRACERS_AEROSOLS_VBS
      use TRACERS_VBS, only: vbs_bins
#endif

      use TracerBundle_mod
      use TracerSource_mod, only: N_MAX_SECT
c     
      IMPLICIT NONE
      SAVE

      type (TracerBundle_type) :: tracers

!@dbparam COUPLED_CHEM: if 0 => uncoupled, if 1 => coupled
      integer :: COUPLED_CHEM = 0

C**** Each tracer has a variable name and a unique index
!@param NTM number of tracers

!@var ntm_O18: Number of TRACERS_SPECIAL_O18 tracers.
#ifdef TRACERS_SPECIAL_O18
      integer, parameter :: ntm_o18=2
#else
      integer, parameter :: ntm_o18=0
#endif  /* TRACERS_SPECIAL_O18 */
!@var ntm_gasexch: Number of TRACERS_GASEXCH_ocean tracers.
#if defined(TRACERS_GASEXCH_ocean) || defined(TRACERS_GASEXCH_land)
      integer, parameter :: ntm_gasexch=1
#else
      integer, parameter :: ntm_gasexch=0
#endif  /* TRACERS_GASEXCH_ocean */
!@var ntm_lerner: Number of TRACERS_SPECIAL_Lerner tracers.
#ifdef TRACERS_SPECIAL_Lerner
      integer, parameter :: ntm_lerner=9
#else
      integer, parameter :: ntm_lerner=0
#endif  /* TRACERS_SPECIAL_Lerner */
!@var ntm_water: Number of TRACERS_WATER tracers.
#if (defined TRACERS_WATER) 
      integer, parameter :: ntm_water=1
#else
      integer, parameter :: ntm_water=0
#endif  /* TRACERS_WATER */
!@var ntm_shindell_trop: Number of TRACERS_SPECIAL_Shindell tracers.
!@var ntm_shindell_strat: Number of Shindell strat chem tracers.
#ifdef TRACERS_SPECIAL_Shindell
      integer, parameter :: ntm_shindell_trop=15
      integer, parameter :: ntm_shindell_strat=10
#else
      integer, parameter :: ntm_shindell_trop=0
      integer, parameter :: ntm_shindell_strat=0
#endif  /* TRACERS_SPECIAL_Shindell */
!@var ntm_terp: Number of TRACERS_TERP tracers.
#ifdef TRACERS_TERP
      integer, parameter :: ntm_terp=1
#else
      integer, parameter :: ntm_terp=0
#endif  /* TRACERS_TERP */
!@var ntm_shindell_extra: Number of SHINDELL_STRAT_EXTRA tracers.
#ifdef SHINDELL_STRAT_EXTRA
#ifdef ACCMIP_LIKE_DIAGS
      integer, parameter :: ntm_shindell_extra=3
#else
      integer, parameter :: ntm_shindell_extra=1
#endif /* ACCMIP_LIKE_DIAGS*/
#else
      integer, parameter :: ntm_shindell_extra=0
#endif  /* SHINDELL_STRAT_EXTRA */
!@var ntm_koch: Number of TRACERS_AEROSOLS_Koch tracers.
!@var ntm_vbs: Number of TRACERS_AEROSOLS_VBS tracers.
#ifdef TRACERS_AEROSOLS_Koch
#ifdef SULF_ONLY_AEROSOLS
#ifdef TRACERS_SPECIAL_Shindell
      integer, parameter :: ntm_koch=4
#else
      integer, parameter :: ntm_koch=5
#endif  /* TRACERS_SPECIAL_Shindell */
      integer, parameter :: ntm_vbs=0
#elif (defined TRACERS_AEROSOLS_VBS)
#ifdef TRACERS_SPECIAL_Shindell
      integer, parameter :: ntm_koch=9
#else
      integer, parameter :: ntm_koch=10
#endif  /* TRACERS_SPECIAL_Shindell */
      integer, parameter :: ntm_vbs=2*vbs_bins
#else
#ifdef TRACERS_SPECIAL_Shindell
      integer, parameter :: ntm_koch=12
#else
      integer, parameter :: ntm_koch=13
#endif  /* TRACERS_SPECIAL_Shindell */
      integer, parameter :: ntm_vbs=0
#endif  /* SULF_ONLY_AEROSOLS or TRACERS_AEROSOLS_VBS */
#else
      integer, parameter :: ntm_koch=0
      integer, parameter :: ntm_vbs=0
#endif  /* TRACERS_AEROSOLS_Koch */
!@var ntm_ococean: Number of TRACERS_AEROSOLS_OCEAN tracers.
#ifdef TRACERS_AEROSOLS_OCEAN
      integer, parameter :: ntm_ococean=1
#else
      integer, parameter :: ntm_ococean=0
#endif  /* TRACERS_AEROSOLS_OCEAN */

!@var ntm_dust: Number of dust aerosol tracers.
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP)||\
    (defined TRACERS_TOMAS) 
#if (defined TRACERS_DUST) || (defined TRACERS_AMP)||\
    (defined TRACERS_TOMAS) 
#ifdef TRACERS_DUST_Silt4
      integer, parameter :: ntm_dust=5
#else
      integer, parameter :: ntm_dust=4
#endif  /* TRACERS_DUST_Silt4 */
#else /* TRACERS_MINERALS || TRACERS_QUARZHEM */
!@var ntm_minerals: Number of TRACERS_MINERALS tracers.
#ifdef TRACERS_MINERALS
      integer, parameter :: ntm_minerals=20
#else
      integer, parameter :: ntm_minerals = 0
#endif  /* TRACERS_MINERALS */
!@var ntm_quarzhem: Number of TRACERS_QUARZHEM tracers.
#ifdef TRACERS_QUARZHEM
      integer, parameter :: ntm_quarzhem=3
#else
      integer, parameter :: ntm_quarzhem = 0
#endif  /* TRACERS_QUARZHEM */
      integer, parameter :: ntm_dust = ntm_minerals + ntm_quarzhem
#endif /* TRACERS_MINERALS || TRACERS_QUARZHEM */
#else
      integer, parameter :: ntm_dust=0
#endif

!@var ntm_het: Number of TRACERS_HETCHEM tracers.
#ifdef TRACERS_HETCHEM
#ifdef TRACERS_NITRATE
      integer, parameter :: ntm_het=6
#else
      integer, parameter :: ntm_het=3
#endif  /* TRACERS_NITRATE */
#else
      integer, parameter :: ntm_het=0
#endif  /* TRACERS_HETCHEM */
!@var ntm_nitrate: Number of TRACERS_NITRATE tracers.
#ifdef TRACERS_NITRATE
      integer, parameter :: ntm_nitrate=3
#else
      integer, parameter :: ntm_nitrate=0
#endif  /* TRACERS_NITRATE */
!@param ntm_soa: Number of SOA tracers.
#ifdef TRACERS_AEROSOLS_SOA
!@+            Index g means gas phase, index p means particulate
!@+            (aerosol) phase. g should ALWAYS be EXACTLY BEFORE
!@+            the corresponding p and all SOA species should be
!@+            one after the other. Do not forget to declare them
!@+            to the chemistry species as well.
#ifdef TRACERS_TERP
      integer, parameter :: ntm_soa=8
#else
      integer, parameter :: ntm_soa=4
#endif  /* TRACERS_TERP */
!@param nsoa the total number of aerosol-phase SOA related species (=ntm_soa/2)
      integer, parameter :: nsoa=ntm_soa/2
#else
      integer, parameter :: ntm_soa=0
#endif  /* TRACERS_AEROSOLS_SOA */
!@var ntm_nitrate: Number of TRACERS_COSMO tracers.
!@var ntm_nitrate: Number of TRACERS_RADON tracers.
#ifdef TRACERS_COSMO
#ifdef TRACERS_RADON
      integer, parameter :: ntm_cosmo=4
#else
      integer, parameter :: ntm_cosmo=2
#endif  /* TRACERS_RADON */
#else
      integer, parameter :: ntm_cosmo=0
#endif  /* TRACERS_COSMO */
!@var ntm_ocean: Number of TRACERS_OCEAN tracers.
#ifdef TRACERS_OCEAN
      integer, parameter :: ntm_ocean=0
#else
      integer, parameter :: ntm_ocean=0
#endif  /* TRACERS_OCEAN */
!@var ntm_air: Number of TRACERS_AIR tracers.
#if defined TRACERS_AIR || defined HTAP_LIKE_DIAGS
      integer, parameter :: ntm_air=1
#else
      integer, parameter :: ntm_air=0
#endif  /* TRACERS_AIR */
#ifdef TRACERS_AMP
#ifdef TRACERS_AMP_M1
      integer, parameter :: ntmAMP=51
#endif  /* TRACERS_AMP_M1 */
#ifdef TRACERS_AMP_M2
      integer, parameter :: ntmAMP=51
#endif  /* TRACERS_AMP_M2 */
#ifdef TRACERS_AMP_M3
      integer, parameter :: ntmAMP=41
#endif  /* TRACERS_AMP_M3 */
#ifdef TRACERS_AMP_M4
      integer, parameter :: ntmAMP=34
#endif  /* TRACERS_AMP_M4 */
#ifdef TRACERS_AMP_M5
      integer, parameter :: ntmAMP=45
#endif  /* TRACERS_AMP_M5 */
#ifdef TRACERS_AMP_M6
      integer, parameter :: ntmAMP=45
#endif  /* TRACERS_AMP_M6 */
#ifdef TRACERS_AMP_M7
      integer, parameter :: ntmAMP=35
#endif  /* TRACERS_AMP_M7 */
#ifdef TRACERS_AMP_M8
      integer, parameter :: ntmAMP=28
#endif  /* TRACERS_AMP_M8 */
#ifdef TRACERS_SPECIAL_Shindell
      integer, parameter :: ntm_amp=ntmAMP+4
#else
      integer, parameter :: ntm_amp=ntmAMP+5
#endif  /* TRACERS_SPECIAL_Shindell */
#else
      integer, parameter :: ntm_amp=0
#endif  /* TRACERS_AMP */

!@param ntm_chem number of drew-only tracers
      integer, parameter :: ntm_chem=ntm_shindell_trop+
     *                               ntm_terp+
     *                               ntm_shindell_strat+
     *                               ntm_soa
#ifdef TRACERS_AMP
! This is kept seperate, as ntm_dust needs to be set 
c          (in order to calculate dust emissions), but not added to ntm.
      integer, parameter :: oldNTM=ntm_amp+ntm_chem
#else
#ifdef TRACERS_TOMAS
       !constants that have to do with the number of tracers
C    NBS is the number of bulk species (gases and aerosols that don't
C	 have size resolution such as MSA)
C    NAP is the number of size-resolved "prognostic" aerosol species
C	 (ones that undergo transport)
C    NAD is the number of size-resolved "diagnostic" aerosol species
C	 (ones that don't undergo transport or have a complete budget
C	 such as aerosol water and nitrate)
C    NSPECIES is the number of unique chemical species not counting
C	 size-resolved species more than once (the sum of NBS, NAP,
C	 and NAD)
C    NBINS is the number of bins used to resolve the size distribution
C    NTM is the total number of tracer concentrations that the model
C	 tracks.  This counts bulk species, and both prognostic and 
C	 diagnostic aerosols.  Each size-resolved aerosol has a number
C	 of tracers equal to NBINS to resolve its mass distribution.
C	 An additional NBINS are required to resolve the aerosol number
C	 distribution.
C    NTT is the total number of transported tracers.  This is the same
C	 as NTM, but excludes the "diagnostic" aerosol species which
C	 do not undergo transport - ??? will I use this

       !constants that determine the size of diagnostic arrays
C    NXP is the number of transport processes that are tracked
C	 separately
C    NCR is the number of chemical reactions for which data is saved
C    NOPT is the number of aerosol optical properties tracked
C    NFOR is the number of different forcings that are calculated
C    NAERO is the number of aerosol microphysics diagnostics
C    NCONS is the number of conservation quantity diagnostics
C    KCDGN is the number of cloud microphysics and optical depth diagnostics

!yhl - nbs should be deleted, and use shindell gas species..

!@param ntm number of tracers

! yhl - NAP is changed to 6 to exclude mineral dust for now (1/20/2011)
#if (defined TOMAS_12_10NM) 
      integer, parameter :: NBINS=12 
#elif (defined TOMAS_15_10NM) || (defined TOMAS_12_3NM)
      integer, parameter :: NBINS=15 
#endif
      integer, parameter :: NBS=7,NAP=7, NAD=1 !, NXP=7, NCR=2, 
      integer, parameter :: ntm_tomas=NBINS*(NAP+NAD+1)

      integer, parameter :: non_aerosol=ntm_O18+ntm_gasexch+ntm_lerner+
     *                          ntm_water+ntm_koch+ntm_vbs+ntm_het+  !exclude ntm_dust! 
     *                          ntm_nitrate+ntm_cosmo+
     *                          ntm_ocean+ntm_air+ntm_chem+
     *                          ntm_shindell_extra+ntm_ococean+NBS

      integer, parameter :: 
c     &     IDTNUMD = non_aerosol+1,         !NBINS for number distribution
     *     IDTSO4  = non_aerosol+1, !36;NBINS for sulfate mass dist.
     &     IDTNA   = IDTSO4 +NBINS, !66;
     &     IDTECOB = IDTNA+NBINS, !126; NBINS for Hydrophobic EC
     &     IDTECIL = IDTECOB + NBINS, !96; NBINS for Hydrophillic EC
     &     IDTOCOB = IDTECIL+NBINS, !186
     &     IDTOCIL = IDTOCOB+NBINS, !156; OC
     &     IDTDUST = IDTOCIL+NBINS, !216
     &     IDTNUMD = IDTDUST+NBINS,
     &     IDTH2O  = IDTNUMD+NBINS  !246
      double precision, dimension(nbins+1) :: xk
!      integer, parameter :: oldNTM=ntm_tomas
      integer, parameter :: oldNTM=non_aerosol+ntm_tomas !ntm_dust is excluded.      

#else
!@param ntm number of tracers
      integer, parameter :: oldNTM=ntm_O18+ntm_gasexch+ntm_lerner+
     *                          ntm_water+ntm_koch+ntm_vbs+ntm_dust+
     *                          ntm_het+ntm_nitrate+ntm_cosmo+
     *                          ntm_ocean+ntm_air+ntm_chem+
     *                          ntm_shindell_extra+ntm_ococean

#endif  /* TRACERS_TOMAS */
#endif

      integer :: NTM

#ifdef TRACERS_TOMAS
      integer, dimension(nbins) :: n_ANUM =0
      integer, dimension(nbins) :: n_ASO4 =0
      integer, dimension(nbins) :: n_ANACL=0
      integer, dimension(nbins) :: n_AECIL=0
      integer, dimension(nbins) :: n_AECOB=0
      integer, dimension(nbins) :: n_AOCIL=0
      integer, dimension(nbins) :: n_AOCOB=0
      integer, dimension(nbins) :: n_AH2O=0
      integer, dimension(nbins) :: n_ADUST=0  
      integer :: n_SOAgas=0
#endif

!@var N_XXX: variable names of indices for tracers (init = 0)
      integer ::
     *                 n_SF6_c=0,                                        
     *     n_Air=0,    n_SF6=0,   n_Rn222=0, n_CO2=0,      n_N2O=0,
     *     n_CFC11=0,  n_14CO2=0, n_CH4=0,   n_O3=0,       n_water=0,
     *     n_H2O18=0,  n_HDO=0,   n_HTO=0,   n_Ox=0,       n_NOx=0,
     *     n_N2O5=0,   n_HNO3=0,  n_H2O2=0,  n_CH3OOH=0,   n_HCHO=0,
     *     n_HO2NO2=0, n_CO=0,    n_PAN=0,   n_H2O17=0,
     *     n_Isoprene=0, n_AlkylNit=0, n_Alkenes=0, n_Paraffin=0,
     *     n_stratOx=0, n_Terpenes=0,n_codirect=0,
     *     n_isopp1g=0,n_isopp1a=0,n_isopp2g=0,n_isopp2a=0,
     *     n_apinp1g=0,n_apinp1a=0,n_apinp2g=0,n_apinp2a=0,
     *     n_DMS=0,    n_MSA=0,   n_SO2=0,   n_SO4=0,    n_H2O2_s=0,
     *     n_ClOx=0,   n_BrOx=0,  n_HCl=0,   n_HOCl=0,   n_ClONO2=0,
     *     n_HBr=0,    n_HOBr=0,  n_BrONO2=0,n_CFC=0,    n_GLT=0,
     *     n_Pb210 = 0,n_Be7=0,   n_Be10=0,
     .     n_CFCn=0,   n_CO2n=0,  n_Age=0,
     *     n_seasalt1=0,  n_seasalt2=0, n_SO4_d1=0,  n_SO4_d2=0,
     *     n_SO4_d3=0,n_N_d1=0,  n_N_d2=0,  n_N_d3=0,
     &     n_NH3=0,   n_NH4=0,   n_NO3p=0,
     *     n_BCII=0,  n_BCIA=0,  n_BCB=0,
     *     n_OCII=0,  n_OCIA=0,  n_OCB=0,
     *     n_vbsGm2=0, n_vbsGm1=0, n_vbsGz=0,  n_vbsGp1=0, n_vbsGp2=0,
     *     n_vbsGp3=0, n_vbsGp4=0, n_vbsGp5=0, n_vbsGp6=0,
     *     n_vbsAm2=0, n_vbsAm1=0, n_vbsAz=0,  n_vbsAp1=0, n_vbsAp2=0,
     *     n_vbsAp3=0, n_vbsAp4=0, n_vbsAp5=0, n_vbsAp6=0,
     *     n_OCocean=0,
     &     n_clay=0,   n_silt1=0, n_silt2=0, n_silt3=0, n_silt4=0,
     &     n_clayilli=0,n_claykaol=0,n_claysmec=0,n_claycalc=0,
     &     n_clayquar=0,n_sil1quar=0,n_sil1feld=0,n_sil1calc=0,
     &     n_sil1hema=0,n_sil1gyps=0,n_sil2quar=0,n_sil2feld=0,
     &     n_sil2calc=0,n_sil2hema=0,n_sil2gyps=0,n_sil3quar=0,
     &     n_sil3feld=0,n_sil3calc=0,n_sil3hema=0,n_sil3gyps=0,
     &     n_sil1quhe=0,n_sil2quhe=0,n_sil3quhe=0,
     *     n_M_NO3=0,   n_M_NH4=0,   n_M_H2O=0,   n_M_AKK_SU=0,
     *     n_N_AKK_1=0, n_M_ACC_SU=0,n_N_ACC_1=0, n_M_DD1_SU=0,
     *     n_M_DD1_DU=0,n_N_DD1_1=0, n_M_DS1_SU=0,n_M_DS1_DU=0,
     *     n_N_DS1_1 =0,n_M_DD2_SU=0,n_M_DD2_DU=0,n_N_DD2_1 =0,
     *     n_M_DS2_SU=0,n_M_DS2_DU=0,n_N_DS2_1 =0,n_M_SSA_SU=0,
     *     n_M_SSA_SS=0,n_M_SSC_SS=0,
     *     n_M_OCC_SU=0,n_M_OCC_OC=0,n_N_OCC_1 =0,
     *     n_M_BC1_SU=0,n_M_BC1_BC=0,n_N_BC1_1 =0,n_M_BC2_SU=0,
     *     n_M_BC2_BC=0,n_N_BC2_1 =0,n_M_BC3_SU=0,n_M_BC3_BC=0,
     *     n_N_BC3_1 =0,n_M_DBC_SU=0,n_M_DBC_BC=0,n_M_DBC_DU=0,
     *     n_N_DBC_1 =0,n_M_BOC_SU=0,n_M_BOC_BC=0,n_M_BOC_OC=0,
     *     n_N_BOC_1=0, n_M_BCS_SU=0,n_M_BCS_BC=0,n_N_BCS_1 =0,
     *     n_M_MXX_SU=0,n_M_MXX_BC=0,n_M_MXX_OC=0,n_M_MXX_DU=0,
     *     n_M_MXX_SS=0,n_N_MXX_1 =0,n_M_OCS_SU=0,n_M_OCS_OC=0,
     *     n_N_OCS_1=0,n_M_SSS_SS=0,n_M_SSS_SU=0,
     *     n_H2SO4=0, n_N_SSA_1=0, n_N_SSC_1=0
#ifdef TRACERS_AMP
!@var ntmAMPi Index of the first AMP tracer
!@var ntmAMPe Index of the last AMP tracer
      integer :: ntmAMPi=0,ntmAMPe=0
#endif

C**** standard tracer and tracer moment arrays

!@var TRM: Tracer array (kg)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: trm

!@var TRMOM: Second order moments for tracers (kg)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:,:) :: trmom

!@var TRDN1: lowest level downdraft tracer concentration (kg/kg)
       REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: trdn1

#ifdef TRACERS_WATER
!@var TRWM tracer in cloud liquid water amount (kg)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: trwm
#endif

!@var daily_z: altitude of model levels (m), updated once per day
!@+   and used for vertical distribution of 3D emissions
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: daily_z


#ifdef TRACERS_WATER
!@param nWD_TYPES number of tracer types for wetdep purposes
      integer, parameter :: nWD_TYPES=3 !(gas,particle,water)
!@param tr_evap_fact fraction of re-evaporation of raindrops by tracer type
C note, tr_evap_fact is not dimensioned as NTM:
      REAL*8, parameter, dimension(nWD_TYPES) :: tr_evap_fact=
     *     (/1.d0, 0.5d0,  1.d0/)
#endif

C**** Diagnostic indices and meta-data

!@var ntsurfsrcmax maximum number of surface 2D sources/sinks
      integer, parameter :: ntsurfsrcmax=16
!@var nt3Dsrcmax maximum number of 3D tracer sources/sinks
      integer, parameter :: nt3Dsrcmax=7
!@var sfc_src array holds tracer sources that go into trsource( )
!@+ maybe wasteful of memory, but works for now...
      real*8, allocatable, dimension(:,:,:,:) :: sfc_src

!@var MTRACE: timing index for tracers
!@var MCHEM: timing index for chemistry (if needed)
      integer mtrace, mchem
!@var MTRADV: timing index for tracer advection (if requested)
#ifdef TRAC_ADV_CPU
      integer mtradv
#endif
!@dbparam no_emis_over_ice, switch whether (>0) or not (<=0) to
!@+ set surface tracer emissions to zero over >90% ice boxes
      integer :: no_emis_over_ice=0

C**** EVERYTHING BELOW HERE IS TRACER SPECIFIC. PLEASE THINK 
C**** ABOUT MOVING IT ELSEWHERE

!@var 3D on-line radical array for interactive aerosol and gas
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: oh_live
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: no3_live

#ifdef TRACERS_SPECIAL_O18
C**** Water isotope specific parameters

!@dbparam supsatfac factor controlling super saturation for isotopes
      real*8 :: supsatfac = 2d-3
#endif

#if (defined TRACERS_SPECIAL_Shindell) || (defined TRACERS_AEROSOLS_Koch) ||\
    (defined TRACERS_AMP) || (defined TRACERS_TOMAS) 
C**** Chemistry specific 3-D arrays

!@var RSULF1, RSULF2, RSULF3, RSULF4: rate coefficients
c for gas phase sulfur chemistry used by aerosol and chemistry models
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:)::rsulf1,rsulf2,rsulf3,rsulf4
#endif

#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP) ||\
    (defined TRACERS_TOMAS) 
C**** Aerosol specific switches and arrays

!!@dbparam OFFLINE_DMS_SS is 0 for standard case, 1 for offline dms, seasalt emission
      integer :: OFFLINE_DMS_SS = 0
!!@dbparam OFFLINE_SS is 0 for standard case, 1 for offline seasalt emission
      integer :: OFFLINE_SS = 0
!@dbparam aer_int_yr indicates year of emission
      integer :: aer_int_yr = 0
!@var SNFST0,TNFST0 are instantaneous SW, LW aerosol forcings for AEROCOM
c      real*8 SNFST0(2,NTM,IM,JM),TNFST0(2,NTM,IM,JM)
#endif

C**** tracer specific switches

!@dbparam rnsrc is a switch for Radon sources 
!@+       0=standard, 1=Conen&Robertson, 2=modified Schery&Wasiolek
      integer :: rnsrc = 0

C**** arrays that could be general, but are only used by chemistry

!! dbparam trans_emis_overr_yr year for overriding tracer transient emis
!! dbparam trans_emis_overr_day day for overriding tracer transient emis
!@var trans_emis_overr_yr year for overriding tracer transient emis
!@var trans_emis_overr_day day for overriding tracer transient emis
      integer :: trans_emis_overr_yr=0, trans_emis_overr_day=0
! ---- section for altering tracers sources by sector/region ----
!@param n_max_reg  maximum number of regions for emissions altering
      integer, parameter :: n_max_reg=10
!@var num_regions the number of source-altering regions from rundeck
!@var num_sectors the number of source-altering sectors from rundeck
      integer :: num_regions, num_sectors
!@var alter_sources true if any source altering factors are on
      logical :: alter_sources
!@var reg_N the north edge of rectangular regions for emissions altering
!@var reg_S the south edge of rectangular regions for emissions altering
!@var reg_E the east  edge of rectangular regions for emissions altering
!@var reg_W the west  edge of rectangular regions for emissions altering
      real*8, dimension(n_max_reg) :: reg_N,reg_S,reg_E,reg_W
!@var sect_name array hold the sector names (all)
      character*10,dimension(N_MAX_SECT):: sect_name
!@var ef_fact the actual factors that alter sources by region/sector
!@var ef_fact3d factors used to alter 3D sources (these are more
!@+ hard-coded for now...)
      real*8, dimension(N_MAX_SECT,n_max_reg) :: ef_fact,ef_fact3D
! variables for outputting a map of the regions:
      real*8, allocatable, dimension(:,:) :: ef_REG_IJ
! --- end of source-altering section ----------------------------
!@param nChemistry index for tracer chemistry 3D source
!@param nOverwrite index for tracer overwrite 3D source
!@param nOther index for tracer misc. 3D source
!@param nAircraft index for tracer aircraft 3D source
!@param nBiomass index for tracer biomass burning 3D source
!@param nVolcanic index for tracer volcano 3D source
!@param nChemloss index for tracer chemistry 3D loss
! Must be a better way to do this, but for now, it is better than
! hardcoding indicies with a number like "3":
      INTEGER, PARAMETER :: nChemistry = 1, nOverwrite = 2,
     &     nOther = 3, nAircraft = 4, nBiomass = 5,
     &     nVolcanic = 6, nChemloss = 7

#ifdef TRACERS_GASEXCH_ocean
#ifdef TRACERS_GASEXCH_ocean_CFC
!@var ocmip_cfc: CFC-11 emissions estimated from OCMIP surf.conc.
      !60years (1939--1998) OCMIP surfc. concentr. converted to
      !global averaged emission rates
      !each value corresponds to the annual value
      REAL*8, DIMENSION(67,NTM) :: ocmip_cfc
#endif
#endif

#if (defined TRACERS_HETCHEM) || (defined TRACERS_NITRATE)
      integer, parameter :: rhet=3
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: rxts,rxts1,rxts2,rxts3
     *                                         ,rxts4
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:,:) :: krate
#endif

!@var xyz_count,xyz_list count/list of tracers in category xyz.
!@+   A tracer can belong to more than one list.
c note: not applying CPP when declaring counts/lists.
      integer ::  ! counts of tracers meeting the following criteria:
     &      active_count  ! itime_tr0 <= itime
     &     ,gases_count   ! tr_wd_type == nGas
     &     ,aero_count    ! tr_wd_type == nPart
     &     ,water_count   ! tr_wd_type == nWater
     &     ,hlawt_count   ! tr_wd_type == nGas and tr_DHD != 0
     &     ,aqchem_count  ! participates in cloud aqueous chemistry
      integer, dimension(:), allocatable ::
     &     active_list,gases_list,aero_list,water_list,
     &     hlawt_list,aqchem_list

      ! temporary support of legacy interface
      interface ntsurfsrc
        module procedure ntsurfsrc_1
        module procedure ntsurfsrc_all
      end interface ntsurfsrc

      contains

      subroutine initTracerCom()
      use TracerBundle_mod, only: TracerBundle

      tracers = TracerBundle()
      call initTracerMetadata()

      end subroutine initTracerCom

      subroutine remake_tracer_lists()
      use OldTracer_mod, only: trname
!@sum regenerates the counts and lists of tracers in various categories
      use model_com, only : itime
      implicit none
      integer :: n,nactive
      integer, dimension(1000) ::
     &     tmplist_active,tmplist_gases,tmplist_aero,tmplist_water,
     &     tmplist_hlawt,tmplist_aqchem
      active_count = 0
      gases_count = 0
      aero_count = 0
      water_count = 0
      hlawt_count = 0
      aqchem_count = 0
      do n=1,NTM

        if(itime.lt.itime_tr0(n)) cycle

        active_count = active_count + 1
        tmplist_active(active_count) = n

        select case(tr_wd_type(n))
        case(nGAS)
          gases_count = gases_count + 1
          tmplist_gases(gases_count) = active_count
          if(tr_DHD(n).ne.0.) then
            hlawt_count = hlawt_count + 1
            tmplist_hlawt(hlawt_count) = active_count
          endif
        case(nPART)
          aero_count = aero_count + 1
          tmplist_aero(aero_count) = active_count
        case(nWATER)
          water_count = water_count + 1
          tmplist_water(water_count) = active_count
        end select
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)||\
    (defined TRACERS_TOMAS)
        select case (trname(n))
        case('SO2','SO4','H2O2_s','H2O2')
          if ( .not.(
     *         (trname(n).eq."H2O2" .and. coupled_chem.eq.0).or.
     *         (trname(n).eq."H2O2_s" .and. coupled_chem.eq.1)) )
     *         then
            aqchem_count = aqchem_count + 1
            tmplist_aqchem(aqchem_count) = active_count
          end if
        end select
#endif
      enddo
      if(allocated(active_list)) deallocate(active_list)
      allocate(active_list(active_count))
      active_list = tmplist_active(1:active_count)

      if(allocated(gases_list)) deallocate(gases_list)
      allocate(gases_list(gases_count))
      gases_list = tmplist_gases(1:gases_count)

      if(allocated(aero_list )) deallocate(aero_list )
      allocate(aero_list(aero_count))
      aero_list  = tmplist_aero(1:aero_count)

      if(allocated(water_list)) deallocate(water_list)
      allocate(water_list(water_count))
      water_list = tmplist_water(1:water_count)

      if(allocated(hlawt_list)) deallocate(hlawt_list)
      allocate(hlawt_list(hlawt_count))
      hlawt_list = tmplist_hlawt(1:hlawt_count)

#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP) ||\
    (defined TRACERS_TOMAS)
      if(allocated(aqchem_list)) deallocate(aqchem_list)
      allocate(aqchem_list(aqchem_count))
      aqchem_list = tmplist_aqchem(1:aqchem_count)
#endif

      return
      end subroutine remake_tracer_lists

      integer function ntsurfsrc_1(index) result(n)
      use Tracer_mod
      integer, intent(in) :: index

      n = ntsurfsrc(tracers, index)
      
      end function ntsurfsrc_1

      function ntsurfsrc_all() result(n)
      use Tracer_mod
      integer, pointer :: n(:)

      allocate(n(getNumTracers(tracers)))
      n = ntsurfsrc(tracers)
      
      end function ntsurfsrc_all

      subroutine set_ntsurfsrc(index, value)
      use Tracer_mod
      integer, intent(in) :: index
      integer, intent(in) :: value

      call setNtsurfsrc(tracers, index, value)
      end subroutine set_ntsurfsrc

      SUBROUTINE ALLOC_TRACER_COM(grid)
!@sum  To allocate arrays whose sizes now need to be determined at
!@+    run time
!@auth NCCS (Goddard) Development Team
      USE DOMAIN_DECOMP_ATM, ONLY : DIST_GRID, getDomainBounds
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(IN) :: grid

      INTEGER :: J_1H, J_0H, I_1H, I_0H

C****
C**** Extract useful local domain parameters from "grid"
C****
      call getDomainBounds(grid, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)
      I_0H=GRID%I_STRT_HALO
      I_1H=GRID%I_STOP_HALO

      ALLOCATE(   ef_REG_IJ(I_0H:I_1H,J_0H:J_1H) )
      ALLOCATE(     oh_live(I_0H:I_1H,J_0H:J_1H,LM),
     *             no3_live(I_0H:I_1H,J_0H:J_1H,LM),
     *                  trm(I_0H:I_1H,J_0H:J_1H,LM,NTM),
     *                trmom(NMOM,I_0H:I_1H,J_0H:J_1H,LM,NTM),
     *                trdn1(NTM,I_0H:I_1H,J_0H:J_1H),
     *              sfc_src(I_0H:I_1H,J_0H:J_1H,NTM,ntsurfsrcmax))

      ALLOCATE(  daily_z(I_0H:I_1H,J_0H:J_1H,LM) )
      daily_z = 0.

#ifdef TRACERS_WATER
      ALLOCATE(        trwm(I_0H:I_1H,J_0H:J_1H,LM,NTM) )
#endif
#ifdef TRACERS_HETCHEM
      ALLOCATE( rxts(I_0H:I_1H,J_0H:J_1H,LM),
     *          rxts1(I_0H:I_1H,J_0H:J_1H,LM),
     *          rxts2(I_0H:I_1H,J_0H:J_1H,LM),
     *          rxts3(I_0H:I_1H,J_0H:J_1H,LM),
     *          rxts4(I_0H:I_1H,J_0H:J_1H,LM),
     *          krate(I_0H:I_1H,J_0H:J_1H,LM,8,rhet))
#endif
#if (defined TRACERS_SPECIAL_Shindell) || (defined TRACERS_AEROSOLS_Koch) ||\
    (defined TRACERS_AMP) || (defined TRACERS_TOMAS)
      ALLOCATE(  rsulf1(I_0H:I_1H,J_0H:J_1H,LM),
     *           rsulf2(I_0H:I_1H,J_0H:J_1H,LM),
     *           rsulf3(I_0H:I_1H,J_0H:J_1H,LM),
     *           rsulf4(I_0H:I_1H,J_0H:J_1H,LM) )
#endif

      END SUBROUTINE ALLOC_TRACER_COM

      subroutine syncProperty(tracers, property, setValue, values)
      use Dictionary_mod, only: sync_param
      use TracerBundle_mod
      type (TracerBundle_type), intent(inout) :: tracers
      character(len=*) :: property
      interface
        subroutine setValue(n,value)
        integer, intent(in) :: n
        integer, intent(in) :: value
        end subroutine setValue
      end interface
      integer, intent(in) :: values(:)

      integer :: scratch(size(values))
      integer :: n 
      integer :: i

      n = getNumTracers(tracers)
      scratch = values
      call sync_param(property,scratch,n)
      do i = 1, n
         call setValue(i, scratch(i))
      end do

      end subroutine syncProperty

      END MODULE TRACER_COM

