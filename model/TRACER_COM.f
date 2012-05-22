#include "rundeck_opts.h"
      

!@sum  TRACER_COM: Exists alone to minimize the number of dependencies
!@+    This version for simple trace gases/chemistry/isotopes
!@auth Jean Lerner/Gavin Schmidt

      MODULE TRACER_COM
!@sum  TRACER_COM tracer variables
!@auth Jean Lerner
C
      use newTracer_COM, only: getTracerNames, MAXLEN_TRACER_NAME
      USE QUSDEF, only: nmom
      USE RESOLUTION, only: im,jm,lm
      use OldTracer_mod, only: tr_mm
      use OldTracer_mod, only: ntm_power
      use OldTracer_mod, only: t_qlimit
      use OldTracer_mod, only: ntsurfsrc
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

C
      IMPLICIT NONE
      SAVE

!@dbparam COUPLED_CHEM: if 0 => uncoupled, if 1 => coupled
      integer :: COUPLED_CHEM = 0

C**** Each tracer has a variable name and a unique index
!@param NTM number of tracers
!@var TRNAME: Name for each tracer >>> MUST BE LEFT-JUSTIFIED <<<

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
      integer, parameter :: ntm_koch=5
      integer, parameter :: ntm_vbs=0
#elif (defined TRACERS_AEROSOLS_VBS)
      integer, parameter :: ntm_koch=10
      integer, parameter :: ntm_vbs=2*vbs_bins
#else
      integer, parameter :: ntm_koch=13
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
      integer, parameter :: ntm_amp=ntmAMP+5
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

      integer, parameter :: NBS=7,NAP=7, NAD=1,  NBINS=12 !, NXP=7, NCR=2, 
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
      double precision, dimension(nbins) :: xk(nbins+1)
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
      character(len=MAXLEN_TRACER_NAME), allocatable :: trname(:)

#ifdef TRACERS_AMP
#ifdef TRACERS_AMP_M1
      integer, parameter :: AMP_MODES_MAP(ntmAMP)=(/
     *    0 ,0 ,0 ,1 ,1,  !AKK
     *    2 ,2 ,3 ,3 ,3,  !ACC,DD1
     *    4 ,4 ,4 ,5 ,5,  !DS1,DD2
     *    5 ,6 ,6 ,6 ,7,  !DD2,DS2,SSA
     *    7 ,8,           !SSA,SSC
     *    9 ,9 ,9 ,10,10, !OCC,BC1
     *    10,11,11,11,12, !BC1,BC2,BC3
     *    12,12,13,13,13, !BC3,DBC
     *    13,14,14,14,14, !DBC,BOC
     *    15,15,15,16,16, !BCS,MXX
     *    16,16,16,16/)
      integer, parameter :: AMP_NUMB_MAP(ntmAMP)=(/
     *    0 ,0 ,0 ,0 ,1,  !AKK
     *    0 ,2 ,0 ,0 ,3,  !ACC,DD1
     *    0 ,0 ,4 ,0 ,0,  !DS1,DD2
     *    5 ,0 ,0 ,6 ,0,  !DD2,DS2,SSA
     *    0 ,0,           !SSA,SSC
     *    0 ,0 ,9 ,0 ,0 , !OCC,BC1
     *    10,0 ,0 ,11,0 , !BC1,BC2,BC3
     *    0 ,12,0 ,0 ,0 , !BC3,DBC
     *    13,0 ,0 ,0 ,14, !DBC,BOC
     *    0 ,0 ,15,0 ,0 , !BCS,MXX
     *    0 ,0 , 0,16/)
      integer, parameter :: AMP_AERO_MAP(ntmAMP)=(/
     *    1 ,2 ,3 ,4 ,5 ,6 ,7 ,8 ,9 ,10,
     *    11,12,13,14,15,16,17,18,19,20,
     *    21,      24,   26,27,28,29,30,
     *    31,32,33,34,35,36,37,38,39,40,
     *    41,42,43,44,45,46,47,48,49,50,
     *    51,52,53,54   /)

      integer, parameter :: AMP_trm_nm1(ntmAMP)=(/
     *             0 ,         0 ,         0 ,ntm_chem+4 ,ntm_chem+4,  !AKK
     *    ntm_chem+6 ,ntm_chem+6 ,ntm_chem+8 ,ntm_chem+8 ,ntm_chem+8,  !ACC,DD1
     *    ntm_chem+11,ntm_chem+11,ntm_chem+11,ntm_chem+14,ntm_chem+14, !DS1,DD2
     *    ntm_chem+14,ntm_chem+17,ntm_chem+17,ntm_chem+17,ntm_chem+20, !DD2,DS2,SSA
     *    ntm_chem+20,ntm_chem+22,                                     !SSA,SSC
     *    ntm_chem+23,ntm_chem+23,ntm_chem+23,ntm_chem+26,ntm_chem+26, !OCC,BC1
     *    ntm_chem+26,ntm_chem+29,ntm_chem+29,ntm_chem+29,ntm_chem+32, !BC1,BC2,BC3
     *    ntm_chem+32,ntm_chem+32,ntm_chem+35,ntm_chem+35,ntm_chem+35, !BC3,DBC
     *    ntm_chem+35,ntm_chem+39,ntm_chem+39,ntm_chem+39,ntm_chem+39, !DBC,BOC
     *    ntm_chem+43,ntm_chem+43,ntm_chem+43,ntm_chem+46,ntm_chem+46, !BCS,MXX
     *    ntm_chem+46,ntm_chem+46,ntm_chem+46,ntm_chem+46/)
      integer, parameter :: AMP_trm_nm2(ntmAMP)=(/
     *             0 ,         0 ,         0 ,ntm_chem+4 ,ntm_chem+4,  !AKK
     *    ntm_chem+6 ,ntm_chem+6 ,ntm_chem+9 ,ntm_chem+9 ,ntm_chem+9,  !ACC,DD1
     *    ntm_chem+12,ntm_chem+12,ntm_chem+12,ntm_chem+15,ntm_chem+15, !DS1,DD2
     *    ntm_chem+15,ntm_chem+18,ntm_chem+18,ntm_chem+18,ntm_chem+21, !DD2,DS2,SSA
     *    ntm_chem+21,ntm_chem+22,                                     !SSA,SSC
     *    ntm_chem+24,ntm_chem+24,ntm_chem+24,ntm_chem+27,ntm_chem+27, !OCC,BC1
     *    ntm_chem+27,ntm_chem+30,ntm_chem+30,ntm_chem+30,ntm_chem+33, !BC1,BC2,BC3
     *    ntm_chem+33,ntm_chem+33,ntm_chem+37,ntm_chem+37,ntm_chem+37, !BC3,DBC
     *    ntm_chem+37,ntm_chem+41,ntm_chem+41,ntm_chem+41,ntm_chem+41, !DBC,BOC
     *    ntm_chem+44,ntm_chem+44,ntm_chem+44,ntm_chem+50,ntm_chem+50, !BCS,MXX
     *    ntm_chem+50,ntm_chem+50,ntm_chem+50,ntm_chem+50/)
#endif
#ifdef TRACERS_AMP_M2
      integer, parameter :: AMP_MODES_MAP(ntmAMP)=(/
     *    0 ,0 ,0 ,1 ,1,  !AKK
     *    2 ,2 ,3 ,3 ,3,  !ACC,DD1
     *    4 ,4 ,4 ,5 ,5,  !DS1,DD2
     *    5 ,6 ,6 ,6 ,7,  !DD2,DS2,SSA
     *    7 ,8,           !SSA,SSC
     *    9 ,9 ,9 ,10,10, !OCC,BC1
     *    10,11,11,11,12, !BC1,BC2,OSC
     *    12,12,13,13,13, !BC3,DBC
     *    13,14,14,14,14, !DBC,BOC
     *    15,15,15,16,16, !BCS,MXX
     *    16,16,16,16/)
      integer, parameter :: AMP_NUMB_MAP(ntmAMP)=(/
     *    0 ,0 ,0 ,0 ,1,  !AKK
     *    0 ,2 ,0 ,0 ,3,  !ACC,DD1
     *    0 ,0 ,4 ,0 ,0,  !DS1,DD2
     *    5 ,0 ,0 ,6 ,0,  !DD2,DS2,SSA
     *    0 ,0,           !SSA,SSC
     *    0 ,0 ,9 ,0 , 0, !OCC,BC1
     *    10,0 ,0 ,11,0 , !BC1,BC2,OSC
     *    0 ,12,0 ,0 ,0 , !DBC
     *    13,0 ,0 ,0 ,14, !DBC,BOC
     *    0 ,0 ,15,0 ,0 , !BCS,MXX
     *    0 ,0 ,0 ,16/)
      integer, parameter :: AMP_AERO_MAP(ntmAMP)=(/
     *    1 ,2 ,3 ,4 ,5 ,6 ,7 ,8 ,9 ,10,        
     *    11,12,13,14,15,16,17,18,19,20,
     *    21,      24,  26,27,28,29,30,        
     *    31,32,33,34,35,36,37,38,39,40,
     *    41,42,43,44,45,46,47,48,49,50,
     *    51,52,53,54   /)

      integer, parameter :: AMP_trm_nm1(ntmAMP)=(/
     *    0 ,0 ,0 ,4 ,4,  !AKK
     *    6 ,6 ,8 ,8 ,8,  !ACC,DD1
     *    11,11,11,14,14, !DS1,DD2
     *    14,17,17,17,20, !DD2,DS2,SSA
     *    20,22,          !SSA,SSC
     *    23,23,23,26,26, !OCC,BC1
     *    26,29,29,29,32, !BC1,BC2
     *    32,32,35,35,35, !OSC,DBC
     *    35,39,39,39,39, !DBC,BOC
     *    43,43,43,46,46, !BCS,MXX
     *    46,46,46,46/)
      integer, parameter :: AMP_trm_nm2(ntmAMP)=(/
     *    0 ,0 ,0 ,4 ,4,  !AKK
     *    6 ,6 ,9 ,9 ,9,  !ACC,DD1
     *    12,12,12,15,15, !DS1,DD2
     *    15,18,18,18,21, !DD2,DS2,SSA
     *    21,22,          !SSA,SSC
     *    24,24,24,27,27, !OCC,BC1
     *    27,30,30,30,33, !BC1,BC2
     *    33,33,37,37,37, !OCS,DBC
     *    37,41,41,41,41, !DBC,BOC
     *    44,44,44,50,50, !BCS,MXX
     *    50,50,50,50/)
#endif
#ifdef TRACERS_AMP_M3
      integer, parameter :: AMP_MODES_MAP(ntmAMP)=(/
     *    0 ,0 ,0 ,1 ,1,  !AKK
     *    2 ,2 ,3 ,3 ,3,  !ACC,DD1
     *    4 ,4 ,4 ,5 ,5,  !DS1,DD2
     *    5 ,6 ,6 ,6 ,7,  !DD2,DS2,SSA
     *    7 ,8,           !SSA,SSC
     *    9 ,9 ,9 ,10,10, !OCC,BC1
     *    10,11,11,11,12, !BC1,BC2,BOc
     *    12,12,12,13,13, !BOC,MXX
     *    13,13,13,13/)
      integer, parameter :: AMP_NUMB_MAP(ntmAMP)=(/
     *    0 ,0 ,0 ,0 ,1,  !AKK
     *    0 ,2 ,0 ,0 ,3,  !ACC,DD1
     *    0 ,0 ,4 ,0 ,0,  !DS1,DD2
     *    5 ,0 ,0 ,6 ,0,  !DD2,DS2,SSA
     *    0 ,0,           !SSA,SSC
     *    0 ,0 ,9 ,0 ,0 , !OCC,BC1
     *    10,0 ,0 ,11,0 , !BC1,BC2,BOc
     *    0 ,0 ,12,0 ,0 , !BOC,MXX
     *    0 ,0 ,0 ,13/)
      integer, parameter :: AMP_AERO_MAP(ntmAMP)=(/
     *    1 ,2 ,3 ,4 ,5 ,6 ,7 ,8 ,9 ,10,        
     *    11,12,13,14,15      ,18   ,20,
     *    21,22,23,24,25,26,27,28,29,30,        
     *    31,32,33,34,35,36,37,38,39,40,
     *    41,42,43,44  /)

      integer, parameter :: AMP_trm_nm1(ntmAMP)=(/
     *    0 ,0 ,0 ,4 ,4,  !AKK
     *    6 ,6 ,8 ,8 ,8,  !ACC,DD1
     *    11,11,11,14,14, !DS1,DD2
     *    14,17,17,17,20, !DD2,DS2,SSA
     *    20,22,          !SSA,SSC
     *    23,23,23,26,26, !OCC,BC1
     *    26,29,29,29,32, !BC1,BC2,BOC
     *    32,32,32,36,36, !BOC,MXX
     *    36,36,36,36/)

      integer, parameter :: AMP_trm_nm2(ntmAMP)=(/
     *    0 ,0 ,0 ,4 ,4,  !AKK
     *    6 ,6 ,9 ,9 ,9,  !ACC,DD1
     *    12,12,12,15,15, !DS1,DD2
     *    15,18,18,18,21, !DD2,DS2,SSA
     *    21,22,          !SSA,SSC
     *    24,24,24,27,27, !OCC,BC1
     *    27,30,30,30,34, !BC1,BC2,BOC
     *    34,34,34,40,40, !BOC,MXX
     *    40,40,40,40/)
 
#endif
#ifdef TRACERS_AMP_M4
      integer, parameter :: AMP_MODES_MAP(ntmAMP)=(/
     *    0 ,0 ,0 ,1 ,1,  !ACC
     *    2 ,2 ,2 ,3 ,3,  !DD1,DS1
     *    3 ,4 ,4 ,4 ,5,  !DS1,DD2,DS2
     *    5 ,5 ,6 ,6,     !DS2,SSS
     *    7 ,7 ,7 ,8 ,8,  !OCC,BC1
     *    8 ,9 ,9 ,9 ,10, !BC1,BC2,MXX
     *    10,10,10,10,10  !MXX
     *    /)
      integer, parameter :: AMP_NUMB_MAP(ntmAMP)=(/
     *    0 ,0 ,0 ,0 ,1,  !ACC
     *    0 ,0 ,2 ,0 ,0,  !DD1,DS1
     *    3 ,0 ,0 ,4 ,0,  !DS1,DD2,DS2
     *    0 ,5 ,0 ,0,     !DS2,SSS
     *    0 ,0 ,6 ,0 ,0,  !OCC,BC1
     *    7 ,0 ,0 ,8 ,0,  !BC1,BC2,MXX
     *    0, 0, 0 ,0 ,9  !MXX
     *    /)
      integer, parameter :: AMP_AERO_MAP(ntmAMP)=(/
     *    1 ,2 ,3 ,4 ,5 ,6 ,7 ,8 ,9 ,10,        
     *    11,12,13,14,15,16,17,18,19,
     *    21,22,23,24,25,26,27,28,29,30,        
     *    31,32,33,34,35   /)
      integer, parameter :: AMP_trm_nm1(ntmAMP)=(/
     *    0 ,0 ,0 ,4 ,4,  !ACC
     *    6 ,6 ,6 ,9 ,9,  !DD1
     *    9 ,12,12,12,15, !DS1,DD2
     *    15,15,18,18,    !DD2,DS2,SSA
     *    20,20,20,23,23, !SSA,SSC
     *    23,26,26,26,29, !OCC,BC1
     *    29,29,29,29,29/)

      integer, parameter :: AMP_trm_nm2(ntmAMP)=(/
     *    0 ,0 ,0 ,4 ,4,  !ACC
     *    7 ,7 ,7 ,10,10,  !DD1
     *    10,13,13,13,16, !DS1,DD2
     *    16,16,19,19,    !DD2,DS2,SSA
     *    21,21,21,24,24, !SSA,SSC
     *    24,27,27,27,33, !OCC,BC1
     *    33,33,33,33,33/)
 
#endif
#ifdef TRACERS_AMP_M5
      integer, parameter :: AMP_MODES_MAP(ntmAMP)=(/
     *    0 ,0 ,0 ,1 ,1,  !AKK
     *    2 ,2 ,3 ,3 ,3,  !ACC,DD1
     *    4 ,4 ,4 ,       !DS1
     *    5 ,             !SSA
     *    5 ,6 ,          !SSA,SSC
     *    7 ,7 ,7, 8, 8,  !OCC,BC1
     *    8, 9, 9, 9,10,  !BC1,BC2,BC3
     *    10,10,11,11,11, !BC3,DBC
     *    11,12,12,12,12, !DBC,BOC
     *    13,13,13,14,14, !BCS,MXX
     *    14,14,14,14/)
      integer, parameter :: AMP_NUMB_MAP(ntmAMP)=(/
     *    0 ,0 ,0 ,0 ,1,  !AKK
     *    0 ,2 ,0 ,0 ,3,  !ACC,DD1
     *    0 ,0 ,4 ,       !DS1
     *    0 ,             !SSA
     *    0 ,0 ,          !SSA,SSC
     *    0 ,0 ,7, 0, 0,  !OCC,BC1
     *    8, 0, 0, 9, 0,  !BC1,BC2,BC3
     *    0, 10,0 ,0 ,0, !BC3,DBC
     *    11,0 ,0 ,0 ,12, !DBC,BOC
     *    0, 0, 13,0 ,0, !BCS,MXX
     *    0,0,0,14/)
      integer, parameter :: AMP_AERO_MAP(ntmAMP)=(/
     *    1 ,2 ,3 ,4 ,5 ,6 ,7 ,8 ,9 ,10,        
     *    11,12,13,14,15,      18,   20,
     *    21,22,23,24,25,26,27,28,29,30,        
     *    31,32,33,34,35,36,37,38,39,40,
     *    41,42,43,44,45,46,47,48   /)
#endif
#ifdef TRACERS_AMP_M6
      integer, parameter :: AMP_MODES_MAP(ntmAMP)=(/
     *    0 ,0 ,0 ,1 ,1,  !AKK
     *    2 ,2 ,3 ,3 ,3,  !ACC,DD1
     *    4 ,4 ,4 ,       !DS1
     *    5 ,             !SSA
     *    5 ,6 ,          !SSA,SSC
     *    7 ,7 ,7, 8, 8,  !OCC,BC1
     *    8, 9, 9, 9,10,  !BC1,BC2,OCS
     *    10,10,11,11,11, !OCS,DBC
     *    11,12,12,12,12, !DBC,BOC
     *    13,13,13,14,14, !BCS,MXX
     *    14,14,14,14/)
      integer, parameter :: AMP_NUMB_MAP(ntmAMP)=(/
     *    0 ,0 ,0 ,0 ,1,  !AKK
     *    0 ,2 ,0 ,0 ,3,  !ACC,DD1
     *    0 ,0 ,4 ,       !DS1
     *    0 ,             !SSA
     *    0 ,0 ,          !SSA,SSC
     *    0 ,0 ,7, 0, 0,  !OCC,BC1
     *    8, 0, 0, 9, 0,  !BC1,BC2,OCS
     *    0 ,10,0 ,0 ,0 , !OCS,DBC
     *    11,0 ,0 ,0 ,12, !DBC,BOC
     *    0 ,0 ,13,0 ,0 , !BCS,MXX
     *    0 ,0 ,0 ,14/)
      integer, parameter :: AMP_AERO_MAP(ntmAMP)=(/
     *    1 ,2 ,3 ,4 ,5 ,6 ,7 ,8 ,9 ,10,        
     *    11,12,13,14,15,      18,   20,
     *    21,22,23,24,25,26,27,28,29,30,        
     *    31,32,33,34,35,36,37,38,39,40,
     *    41,42,43,44,45,46,47,48   /)
#endif
#ifdef TRACERS_AMP_M7
      integer, parameter :: AMP_MODES_MAP(ntmAMP)=(/
     *    0 ,0 ,0 ,1 ,1,  !AKK
     *    2 ,2 ,3 ,3 ,3,  !ACC,DD1
     *    4 ,4 ,4 ,       !DS1
     *    5 ,             !SSA
     *    5 ,6 ,          !SSA,SSC
     *    7 ,7 ,7, 8, 8,  !OCC,BC1
     *    8, 9, 9, 9,     !BC1,BC2
     *    10,10,10,10,    !OCS,DBC
     *    11,11,          !MXX
     *    11,11,11,11/)   !MXX
      integer, parameter :: AMP_NUMB_MAP(ntmAMP)=(/
     *    0 ,0 ,0 ,0 ,1,  !AKK
     *    0 ,2 ,0 ,0 ,3,  !ACC,DD1
     *    0 ,0 ,4 ,       !DS1
     *    0 ,             !SSA
     *    0 ,0 ,          !SSA,SSC
     *    0 ,0 ,7, 0, 0,  !OCC,BC1
     *    8, 0, 0, 9,     !BC1,BC2
     *    0, 0, 0,10,    !OCS,DBC
     *    0, 0,          !MXX
     *    0, 0, 0,11/)   !MXX
      integer, parameter :: AMP_AERO_MAP(ntmAMP)=(/
     *    1 ,2 ,3 ,4 ,5 ,6 ,7 ,8 ,9 ,10,        
     *    11,12,13,14,15,      18,   20,
     *    21,22,23,24,25,26,27,28,29,30,        
     *    31,32,33,34,35,36,37,38/)
#endif
#ifdef TRACERS_AMP_M8
      integer, parameter :: AMP_MODES_MAP(ntmAMP)=(/
     *    0 ,0 ,0 ,       !
     *    1 ,1 ,2 ,2 ,2,  !ACC,DD1
     *    3 ,3 ,3 ,       !DS1
     *    4 ,             !SSS
     *    4 ,5 ,          !SSS
     *    5 ,5 ,5, 6, 6,  !OCC,BC1
     *    6, 7, 7, 7, 8,  !BC1,BC2
     *    8,8,8,8,8/)    !MXX
      integer, parameter :: AMP_NUMB_MAP(ntmAMP)=(/
     *    0 ,0 ,0 ,       !
     *    0 ,1 ,0 ,0 ,2,  !ACC,DD1
     *    0 ,0 ,3 ,       !DS1
     *    0 ,             !SSS
     *    0 ,0 ,          !SSS
     *    0 ,0 ,5, 0, 0,  !OCC,BC1
     *    6, 0, 0, 7, 0,  !BC1,BC2
     *    0,0,0,0,8/)    !MXX
      integer, parameter :: AMP_AERO_MAP(ntmAMP)=(/
     *    1 ,2 ,3 ,4 ,5 ,        
     *    6 ,7 ,8 ,9 ,10,        
     *    11,12,13,15,16,        
     *    17,18,19,20,21,
     *    22,23,24,25,26,        
     *    27,28,29   /)
#endif
#endif   /* TRACERS_AMP */
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

! The following arrays are currently used in the calculation
! of surface fluxes (via pointers in data structures for
! atm-surface coupling).  More generally, the "layer 1"
! choice could be replaced by the full boundary layer depth.
! TRM1 is currently allocated by the FLUXES allocate routine.
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:), TARGET :: TRM1

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
!@var iso_index indexing taking actual tracer number to isotope
!@+   fractionation number (1=water,2=h2o18,3=hdo,4=hto,5=h2o17)
      integer, allocatable :: iso_index(:)
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
!@dbparam imPI is 0 for industrial simulations, 1 for pre-industrial
      integer :: imPI = 0
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

!@var ssname holds source name, read from file header, e.g. to be
!@+ placed into lname and sname arrays.
!@var freq frequency (annual? monthly?) read from emis file header
!@var res horiz. resolution (M? F?) read from emis file header
!@var nameT tracer name read from emis file header (should match trname)
!@var ty_start starting year/decade for a transient emissions file
!@var ty_end     ending year/decade for a transient emissions file
!! dbparam trans_emis_overr_yr year for overriding tracer transient emis
!! dbparam trans_emis_overr_day day for overriding tracer transient emis
!@var trans_emis_overr_yr year for overriding tracer transient emis
!@var trans_emis_overr_day day for overriding tracer transient emis
      character*30, allocatable, dimension(:,:) :: ssname ! some maybe
      character*10, allocatable, dimension(:,:) :: nameT  ! need not be
      character*1, allocatable, dimension(:,:) :: freq,res ! arrays
      character*9, allocatable, dimension(:,:) :: Tyears   ! here...
      integer, allocatable, dimension(:,:) :: ty_start,ty_end,delTyr
      integer :: trans_emis_overr_yr=0, trans_emis_overr_day=0
! ---- section for altering tracers sources by sector/region ----
!@param n_max_sect maximum number of sectors for emissions altering
!@param n_max_reg  maximum number of regions for emissions altering
      integer, parameter :: n_max_sect=10, n_max_reg=10
!@var num_tr_sectors number of sectors for a particular tracer and source
      integer, allocatable, dimension(:,:) :: num_tr_sectors
!@var num_tr_sectors3D number of sectors for a tracer's 3D source
      integer, allocatable, dimension(:,:) :: num_tr_sectors3D
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
!@var tr_sect_index array hold the sector index for given tracer/source
      integer, allocatable, dimension(:,:,:) :: tr_sect_index
!@var tr_sect_index3D holds 3d source sector index for given tracer/source
      integer, allocatable, dimension(:,:,:) :: tr_sect_index3D
!@var tr_sect_name array hold the sector name for given tracer/source
      character*10,allocatable, dimension(:,:,:):: tr_sect_name
!@var tr_sect_name3D holds 3d source sector name for given tracer/source
      character*10,allocatable, dimension(:,:,:):: tr_sect_name3D
!@var sect_name array hold the sector names (all)
      character*10,dimension(n_max_sect):: sect_name
!@var ef_fact the actual factors that alter sources by region/sector
!@var ef_fact3d factors used to alter 3D sources (these are more
!@+ hard-coded for now...)
      real*8, dimension(n_max_sect,n_max_reg) :: ef_fact,ef_fact3D
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

      contains

      subroutine initTracerCom()
      USE DOMAIN_DECOMP_ATM, only: AM_I_ROOT
      integer :: i

      call getTracerNames(trname)
      NTM = size(trname)

      allocate(ssname(NTM, ntsurfsrcmax))
      allocate(nameT(NTM, ntsurfsrcmax))
      allocate(freq(NTM, ntsurfsrcmax))
      allocate(res(NTM, ntsurfsrcmax))
      allocate(Tyears(NTM, ntsurfsrcmax))
      allocate(ty_start(NTM, ntsurfsrcmax))
      allocate(ty_end(NTM, ntsurfsrcmax))
      allocate(delTyr(NTM, ntsurfsrcmax))
      allocate(num_tr_sectors(NTM, ntsurfsrcmax))
      
      allocate(num_tr_sectors3D(NTM,nt3Dsrcmax))
      allocate(tr_sect_index(NTM,ntsurfsrcmax,n_max_sect))
      allocate(tr_sect_index3D(NTM,nt3Dsrcmax,n_max_sect))
      allocate(tr_sect_name(NTM,ntsurfsrcmax,n_max_sect))
      allocate(tr_sect_name3D(NTM,nt3Dsrcmax,n_max_sect))

#ifdef TRACERS_SPECIAL_O18
      allocate(iso_index(NTM))
#endif

      end subroutine initTracerCom

      subroutine remake_tracer_lists()
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

      END MODULE TRACER_COM

      SUBROUTINE ALLOC_TRACER_COM(grid)
!@sum  To allocate arrays whose sizes now need to be determined at
!@+    run time
!@auth NCCS (Goddard) Development Team
      USE TRACER_COM
      USE DOMAIN_DECOMP_ATM, ONLY : DIST_GRID, GET
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(IN) :: grid

      INTEGER :: J_1H, J_0H, I_1H, I_0H

C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)
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

