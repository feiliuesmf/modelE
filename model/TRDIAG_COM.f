#include "rundeck_opts.h"

      MODULE TRACER_DIAG_COM
!@sum Tracer diagnostic arrays
!@+    Mostly tracer independent, but this may depend on applications
!@auth Jean Lerner
!ver   1.0
      USE MODEL_COM, only: im,jm,lm
      USE DAGCOM, only: npts !npts are conservation quantities
      USE TRACER_COM, only: ntm
#ifdef TRACERS_ON
     *     , ntsurfsrcmax, nt3Dsrcmax
#endif
      IMPLICIT NONE
      SAVE

C**** TAIJS  <<<< KTAIJS and IJTS_xx are Tracer-Dependent >>>>
C**** TAJLS  <<<< KTAJLS and JLS_xx are Tracer-Dependent >>>>

!@dbparam to_volume_MixRat: For printout of tracer concentration
!@+   to_volume_MixRat=1: printout is in Volume Mixing Ratio
!@+   to_volume_MixRat=0: printout is in Mass Mixing Ratio
#ifdef TRACERS_SPECIAL_Shindell
      INTEGER, DIMENSION(NTM) :: to_volume_MixRat=1
#else
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
      real*8, dimension(im,jm,lm,ntm) :: taijln
!@var SNAME_IJT, UNITS_IJT: Names and units of lat-sigma tracer IJ diags
      character(len=30), dimension(lm,ntm) :: sname_ijt,units_ijt
!@var LNAME_IJT: descriptions of tracer IJ diagnostics
      character(len=80), dimension(lm,ntm) :: lname_ijt = 'unused'
!@var SCALE_IJT: printout scaling factor for tracer IJ diagnostics
      REAL*8, dimension(lm,ntm) :: scale_ijt
!@var IR_IJT: range index of IJ diagnosts
      integer, dimension(ntm) :: ir_ijt
!@var IA_IJT: idacc-number for tracer IJ diags (same for all tracers)
      integer ia_ijt
!@var IJTC_POWER: power of 10 used for tracer IJ concentration diags
      integer, dimension(ntm) :: ijtc_power
!@var IJTM_POWER: power of 10 used for tracer IJ mass unit diags
      integer, dimension(ntm) :: ijtm_power
!@var IJT_XX Names for TAIJLN diagnostics (Not currently used)

C**** TAIJN
!@param KTAIJ number of 3D diagnostics for each tracer
#ifdef TRACERS_WATER
#ifdef TRACERS_DRYDEP
      integer, parameter :: ktaij=13
#else
      integer, parameter :: ktaij=11
#endif
#else
#ifdef TRACERS_DRYDEP
      integer, parameter :: ktaij=5
#else
      integer, parameter :: ktaij=3
#endif
#endif
!@var IJT_XX names for taijn diagnostics
      integer tij_conc,tij_surf,tij_mass
#ifdef TRACERS_WATER
!@var IJT_XX names for water-based taijn diagnostics
      integer tij_rvr,tij_seaice,tij_prec,tij_evap,tij_grnd,tij_lk1
     *     ,tij_lk2,tij_soil
#endif
#ifdef TRACERS_DRYDEP
      integer tij_drydep,tij_gsdep
#endif
!@var TAIJN lat/lon tracer diagnostics (all tracers)
      real*8, dimension(im,jm,ktaij,ntm) :: taijn
!@var SCALE_TIJ: printout scaling factor for tracer IJK diagnostics
      REAL*8, dimension(ktaij,ntm) :: scale_tij
!@var SNAME_TIJ,UNITS_TIJ: Names and units of lat-sigma tracer diags
      character(len=30), dimension(ktaij,ntm) :: sname_tij,units_tij
!@var LNAME_TIJ: descriptions of tracer IJK diags
      character(len=80), dimension(ktaij,ntm) :: lname_tij = 'unused'

C**** TAIJS  <<<< KTAIJS and IJTS_xx are Tracer-Dependent >>>>
!@var ijs_XXX index for diags not specific to a certain tracer
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_SPECIAL_Shindell)
      INTEGER ijs_flash,ijs_CtoG    ! ,ijs_OxL1
#ifdef regional_Ox_tracers
      INTEGER ijs_Oxloss, ijs_Oxprod
#endif
#endif

!@param KTAIJS number of special lat/lon tracer diagnostics
      integer, parameter :: ktaijs=90

!@var TAIJS  lat/lon special tracer diagnostics; sources, sinks, etc.
      REAL*8, DIMENSION(IM,JM,ktaijs) :: TAIJS
!@var ijts_source tracer independent array for TAIJS surface src. diags
      INTEGER ijts_source(ntsurfsrcmax,ntm)
!@var ijts_isrc tracer independent array for TAIJS interactive srf. src.
      INTEGER ijts_isrc(ntsurfsrcmax,ntm)
#ifdef TRACERS_AEROSOLS_Koch
!@var ijts_tau tracer independent array for TAIJS hydrated opt. thick.
      INTEGER ijts_tau(ntm)
#endif
!@var ijts_fc tracer independent array for TAIJS SW/LW rad. forcings
      INTEGER ijts_fc(2,ntm)
!@var ijts_3Dsource tracer independent array for TAIJS 3D src. diags
      INTEGER ijts_3Dsource(nt3Dsrcmax,ntm)
!@var SNAME_IJTS, UNITS_IJTS: Names & units of lat-sigma tracer diags
      character(len=30), dimension(ktaijs) :: sname_ijts,units_ijts
!@var LNAME_IJTS: descriptions of tracer IJTS diags
      character(len=80), dimension(ktaijs) :: lname_ijts = 'unused'
!@var SCALE_IJTS: printout scaling factor for tracer IJTS diagnostics
      REAL*8, dimension(ktaijs) :: scale_ijts
!@var IA_IJTS: idacc-number for tracer source/sink IJ diags
      integer ia_ijts(ktaijs)
!@var ijts_power: power of 10 used for tracer IJ source/sink diags
      INTEGER, DIMENSION(ktaijs) :: ijts_power
!@var ijts_index: tracer index associated with a TAIJS diagnostic
      INTEGER, DIMENSION(ktaijs) :: ijts_index

C**** TAJLN
!@param ktajl,ktajlx number of TAJL tracer diagnostics;
!@+          ktajlx includes composites
#ifdef TRACERS_WATER
      INTEGER, PARAMETER :: ktajl=10
#else
      INTEGER, PARAMETER :: ktajl=9
#endif
     &     ,ktajlx=ktajl+2
!@var TAJLN  vertical tracer diagnostics (all tracers)
      REAL*8, DIMENSION(JM,LM,ktajlx,NTM) :: TAJLN
!@var jlnt_xx Names for TAJLN diagnostics
      INTEGER jlnt_conc,jlnt_mass,jlnt_nt_tot,jlnt_nt_mm,jlnt_vt_tot,
     &  jlnt_vt_mm,jlnt_mc,jlnt_turb,jlnt_lscond, jlnt_bebe,
     &  jlnt_bepb
#ifdef TRACERS_WATER
     &     ,jlnt_cldh2o
#endif
!@var SNAME_JLN: Names of lat-sigma tracer JL diagnostics
      character(len=30), dimension(ktajlx,ntm) :: sname_jln
!@var LNAME_JLN,UNITS_JLN: descriptions/units of tracer JL diagnostics
      character(len=50), dimension(ktajlx,ntm) :: lname_jln = 'unused'
     *     ,units_jln
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
!@var jls_XXX index for non-tracer specific or special diags
#ifdef regional_Ox_tracers
      INTEGER jls_Oxloss, jls_Oxprod
#endif
#ifdef TRACERS_AEROSOLS_Koch
      INTEGER jls_OHconk,jls_HO2con,jls_NO3,jls_phot,jls_incloud(2,ntm)
#endif
#ifdef TRACERS_SPECIAL_Shindell
      INTEGER jls_OHcon,jls_H2Omr,jls_N2O5sulf,jls_day
#endif

#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_SPECIAL_Shindell)
#ifdef regional_Ox_tracers
      INTEGER, PARAMETER :: ktajls=108
#else
      INTEGER, PARAMETER :: ktajls=94
#endif
#else
      INTEGER, PARAMETER :: ktajls=36   ! default
#endif
!@var TAJLS  JL special tracer diagnostics for sources, sinks, etc
      REAL*8, DIMENSION(JM,LM,ktajls) :: TAJLS
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
!@var jwt_jls: Weighting index for jls diags 1=by area, 2=simple average
      integer, dimension(ktajls) :: jwt_jls
!@var SNAME_JLS: Names of lat-sigma tracer JL sources/sinks
      character(len=30), dimension(ktajls) :: sname_jls
!@var LNAME_JLS,UNITS_JLS: descriptions/units of tracer JLS diags
      character(len=50), dimension(ktajls) :: lname_jls = 'unused'
     *     ,units_jls
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
      REAL*8, DIMENSION(JM,ktcon,ntmxcon) :: TCONSRV
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
      character(len=30), dimension(ktcon,ntmxcon) ::
     &   name_tconsrv,units_tconsrv
      character(len=80), dimension(ktcon,ntmxcon) :: lname_tconsrv
!@var SCALE_INST,SCALE_CHANGE: Scale factors for tracer conservation
      REAL*8, dimension(ntmxcon) :: SCALE_INST,SCALE_CHANGE
!@var itcon_surf Index array for surface source/sink conservation diags
      INTEGER, DIMENSION(ntsurfsrcmax,ntmxcon) :: itcon_surf
!@var itcon_3D Index array for 3D source/sink conservation diags
      INTEGER, DIMENSION(nt3Dsrcmax,ntmxcon) :: itcon_3Dsrc
!@var itcon_decay Index array for decay conservation diags
      INTEGER, DIMENSION(ntmxcon) :: itcon_decay
!@var itcon_mc Index array for moist convection conserv. diags
      INTEGER, DIMENSION(ntmxcon) :: itcon_mc
!@var itcon_ss Index array for large-scale condensation conserv. diags
      INTEGER, DIMENSION(ntmxcon) :: itcon_ss
#ifdef TRACERS_DRYDEP
!@var itcon_dd Index array for dry deposition conserv. diags
      INTEGER, DIMENSION(ntmxcon,2) :: itcon_dd
#endif
#endif
!@var PDSIGJL temporary storage for mean pressures for jl diags
      REAL*8, DIMENSION(JM,LM) :: PDSIGJL

      END MODULE TRACER_DIAG_COM

#ifdef TRACERS_ON
      SUBROUTINE SET_TCON(QCON,NAME_CON,QSUM,INST_UNIT,SUM_UNIT
     *     ,INST_SC,CHNG_SC, itr,CONPTs)
!@sum  SET_TCON assigns conservation diagnostic array indices
!@auth Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only: sday
      USE MODEL_COM, only: dtsrc,nfiltr
      USE DAGCOM, only: npts,ia_d5d,ia_d5s,ia_filt,ia_12hr,ia_src,conpt0
      USE TRACER_DIAG_COM, only: ktcon,title_tcon,scale_tcon,nsum_tcon
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
      USE MODEL_COM, only: ioread,iowrite,iowrite_mon,iowrite_single
     *     ,irerun,ioread_single,lhead
      USE TRACER_DIAG_COM
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
     *        JM*LM*ktajlx*NTM + JM*LM*ktajls + JM*ntmxcon*ktcon
!@var TA..4(...) dummy arrays for reading diagnostics files
      REAL*4 TAIJLN4(im,jm,lm,ntm),TAIJN4(im,jm,ktaij,ntm)
      REAL*4 TAIJS4(IM,JM,ktaijs),TAJLN4(JM,LM,ktajlx,NTM)
      REAL*4 TAJLS4(JM,LM,ktajls),TCONSRV4(JM,ktcon,ntmxcon)

      write (MODULE_HEADER(lhead+1:80),'(a,i8,a)')
     *   'R8 TACC(',ktacc,'),it'

      SELECT CASE (IACTION)
      CASE (IOWRITE,IOWRITE_MON) ! output to standard restart file
        WRITE (kunit,err=10) MODULE_HEADER,
     *                         TAIJLN,TAIJN,TAIJS,TAJLN,TAJLS,TCONSRV,it
      CASE (IOWRITE_SINGLE)    ! output to acc file
        MODULE_HEADER(LHEAD+1:LHEAD+2) = 'R4'
        WRITE (kunit,err=10) MODULE_HEADER,
     *     REAL(TAIJLN,KIND=4),REAL(TAIJN,KIND=4),REAL(TAIJS,KIND=4),
     *     REAL(TAJLN,KIND=4),REAL(TAJLS,KIND=4),REAL(TCONSRV,KIND=4),it
      CASE (IOREAD:)          ! input from restart file
        SELECT CASE (IACTION)
        CASE (ioread_single)    ! accumulate diagnostic files
          READ (kunit,err=10) HEADER,
     *       TAIJLN4,TAIJN4,TAIJS4,TAJLN4,TAJLS4,TCONSRV4,it_check
          if (it.ne.it_check) then
            PRINT*,"io_trdiag: compare TAIJLN,TAIJLN4, ... dimensions"
            go to 10  ! or should this be just a warning ??
          end if
C****     Accumulate diagnostics (converting back to real*8)
          TAIJLN = TAIJLN+TAIJLN4 ; TAIJN = TAIJN+TAIJN4
          TAIJS  = TAIJS+TAIJS4   ; TAJLN = TAJLN+TAJLN4
          TAJLS  = TAJLS+TAJLS4   ; TCONSRV = TCONSRV+TCONSRV4
          IF (HEADER(1:LHEAD).NE.MODULE_HEADER(1:LHEAD)) THEN
            PRINT*,"Discrepancy in module version ",HEADER
     *           ,MODULE_HEADER
            GO TO 10
          END IF
        CASE (ioread)  ! restarts
          READ (kunit,err=10) HEADER,
     *                         TAIJLN,TAIJN,TAIJS,TAJLN,TAJLS,TCONSRV,it
          IF (HEADER(1:LHEAD).NE.MODULE_HEADER(1:LHEAD)) THEN
            PRINT*,"Discrepancy in module version ",HEADER
     *           ,MODULE_HEADER
            GO TO 10
          end IF
        END SELECT
      END SELECT

      RETURN
 10   IOERR=1
      RETURN
      END SUBROUTINE io_trdiag
#endif

