      MODULE DAGCOM
!@sum  DAGCOM Diagnostic model variables
!@auth Original Development Team
!@var  1.0
      USE MODEL_COM, only : im,jm,lm,imh,ntype
      USE RADNCB, only : LM_REQ

      IMPLICIT NONE
      SAVE
C**** ACCUMULATING DIAGNOSTIC ARRAYS
!@param KAJ number of accumulated zonal budget diagnostics
      INTEGER, PARAMETER :: KAJ=94
!@var AJ zonal budget diagnostics for each surface type
      DOUBLE PRECISION, DIMENSION(JM,KAJ,NTYPE) :: AJ

!@param NREG number of regions for budget diagnostics
      INTEGER, PARAMETER :: NREG=24
!@var AREG regional budget diagnostics
      DOUBLE PRECISION, DIMENSION(NREG,KAJ) :: AREG
!@var TITREG,NAMREG title and names of regions for AREG diagnostics
      CHARACTER*4 TITREG*80,NAMREG(2,23)
!@var JREH lat/lon array defining regions for AREG diagnostics
      INTEGER, DIMENSION(IM,JM) :: JREG

!@param KAPJ number of zonal pressure diagnostics
      INTEGER, PARAMETER :: KAPJ=3  ! why isn't this 2?
!@var APJ zonal pressure diagnostics
      DOUBLE PRECISION, DIMENSION(JM,KAPJ) :: APJ

!@param KAJL number of accumulated AJL diagnostics
      INTEGER, PARAMETER :: KAJL=57
!@var AJL latitude/height diagnostics
      DOUBLE PRECISION, DIMENSION(JM,LM,KAJL) :: AJL

!@param KASJL number of ASJL diagnostics
      INTEGER, PARAMETER :: KASJL=4
!@var ASJL latitude/height supplementary diagnostics (merge with AJL?)
      DOUBLE PRECISION, DIMENSION(JM,LM_REQ,KASJL) :: ASJL

!@param KAIJ number of AIJ diagnostics
      INTEGER, PARAMETER :: KAIJ=150
!@var AIJ latitude/longitude diagnostics
      DOUBLE PRECISION, DIMENSION(IM,JM,KAIJ) :: AIJ

!@param KAIL number of AIL diagnostics
      INTEGER, PARAMETER :: KAIL=16
!@var AIL longitude/height diagnostics
      DOUBLE PRECISION, DIMENSION(IM,LM,KAIL) :: AIL

C NEHIST = (TROPO/STRAT)X(ZKE/EKE/SEKE/ZPE/EPE)X(SH/NH)
!@param NEHIST,HISTDAYS number of energy history diagnostics, and days
      INTEGER, PARAMETER :: NEHIST=20
      INTEGER, PARAMETER :: HIST_DAYS=100
!@var ENERGY energy diagnostics
      DOUBLE PRECISION, DIMENSION(NEHIST,HIST_DAYS) :: ENERGY

!@var NPTS number of points at which standard conserv. diags are called
      INTEGER, PARAMETER :: NPTS = 11  ! 9
!@param NQUANT Number of conserved quantities in conservation diags
      INTEGER, PARAMETER :: NQUANT=18   ! 20?
!@param KCON number of conservation diagnostics
      INTEGER, PARAMETER :: KCON=125    !54
!@var CONSRV conservation diagnostics
      DOUBLE PRECISION, DIMENSION(JM,KCON) :: CONSRV
!@var SCALE_CON scales for conservation diagnostics
      DOUBLE PRECISION, DIMENSION(KCON) :: SCALE_CON
!@var TITLE_CON titles for conservation diagnostics
      CHARACTER*32, DIMENSION(KCON) :: TITLE_CON
!@var NSUM_CON indices for summation of conservation diagnostics
!@var IA_CON IDACC numbers for conservation diagnostics
      INTEGER, DIMENSION(KCON) :: NSUM_CON, IA_CON
!@var NOFM indices for CONSRV array
      INTEGER, DIMENSION(NPTS+1,NQUANT) :: NOFM
!@var icon_xx indexes for conservation quantities
      INTEGER icon_AM,icon_KE,icon_MS,icon_TPE,icon_WM,icon_LKM,icon_LKE
     *     ,icon_EWM,icon_WTG,icon_HTG,icon_OCE,icon_MSI,icon_HSI
     *     ,icon_SSI 
!@var KCMX actual number of conservation diagnostics
      INTEGER :: KCMX = 23 ! take up first 23 indexes for special cases

!@param KSPECA,NSPHER number of spectral diagnostics, and harmonics used
      INTEGER, PARAMETER :: KSPECA=20
      INTEGER, PARAMETER :: NSPHER=8
!@var SPECA spectral diagnostics
      DOUBLE PRECISION, DIMENSION((IMH+1),KSPECA,NSPHER) :: SPECA

!@param KTPE number of spectral diagnostics for pot. enthalpy
      INTEGER, PARAMETER :: KTPE=8
      integer, parameter :: NHEMI=2
!@var ATPE pot. enthalpy spectral diagnostics
      DOUBLE PRECISION, DIMENSION(KTPE,NHEMI) :: ATPE

!@param HR_IN_DAY hours in day
      INTEGER, PARAMETER :: HR_IN_DAY=24
!@param NDLYVAR number of daily diagnostics
      INTEGER, PARAMETER :: NDLYVAR=63
!@param NDLYPT number of points where daily diagnostics are kept
      INTEGER, PARAMETER :: NDLYPT=4
!@var ADAILY daily diagnostics (24 hour cycles at selected points)
      DOUBLE PRECISION, DIMENSION(HR_IN_DAY,NDLYVAR,NDLYPT) :: ADAILY

!@param KAJK number of zonal constant pressure diagnostics
      INTEGER, PARAMETER :: KAJK=51
!@var AJK zonal constant pressure diagnostics
      DOUBLE PRECISION, DIMENSION(JM,LM,KAJK) :: AJK

!@param KAIJK number of lat/lon constant pressure diagnostics
      INTEGER, PARAMETER :: KAIJK=6
!@var KAIJK lat/lon constant pressure diagnostics
      DOUBLE PRECISION, DIMENSION(IM,JM,LM,KAIJK) :: AIJK

!@param KAIJL number of three dimensional diagnostics
      INTEGER, PARAMETER :: KAIJL=5
!@var AIJL three dimensional diagnostics
      DOUBLE PRECISION, DIMENSION(IM,JM,LM,KAIJL) :: AIJL

!@param NWAV_DAG number of components in spectral diagnostics
      INTEGER, PARAMETER :: NWAV_DAG=9   !min(9,imh)
!@param KAJLSP number of spectral diagnostics
      INTEGER, PARAMETER :: KAJLSP=3
!@var AJLSP spectral diagnostics
      DOUBLE PRECISION, DIMENSION(JM,LM,0:NWAV_DAG,KAJLSP) :: AJLSP

!@param Max12HR_sequ,Min12HR_sequ lengths of time series for wave powers
      INTEGER, PARAMETER :: Max12HR_sequ=2*31, Min12HR_sequ=2*28
!@param RE_AND_IM complex components of wave power diagnostics
      INTEGER, PARAMETER :: RE_AND_IM=2
!@param KWP number of wave power diagnostics
      INTEGER, PARAMETER :: KWP=12
!@var WAVE frequency diagnostics (wave power)
      DOUBLE PRECISION,
     &     DIMENSION(RE_AND_IM,Max12HR_sequ,NWAV_DAG,KWP) :: WAVE

!@param KACC total number of diagnostic elements
      INTEGER, PARAMETER :: KACC= JM*KAJ*NTYPE + NREG*KAJ
     *     + JM*KAPJ + JM*LM*KAJL + JM*LM_REQ*KASJL + IM*JM*KAIJ +
     *     IM*LM*KAIL + NEHIST*HIST_DAYS + JM*KCON +
     *     (IMH+1)*KSPECA*NSPHER + KTPE*NHEMI + HR_IN_DAY*NDLYVAR*NDLYPT
     *     + RE_AND_IM*Max12HR_sequ*NWAV_DAG*KWP + JM*LM*KAJK +
     *     IM*JM*LM*KAIJK + IM*JM*LM*KAIJL + JM*LM*(1+NWAV_DAG)*KAJLSP

!@var TSFREZ freezing temperature diagnostics
      integer, parameter :: ktsf=4
      DOUBLE PRECISION, DIMENSION(IM,JM,KTSF) :: TSFREZ

!@param KTD number of diurnal temperature diagnostics
      INTEGER, PARAMETER :: KTD=8
!@var TDIURN diurnal range temperature diagnostics
      DOUBLE PRECISION, DIMENSION(IM,JM,KTD) :: TDIURN

!@var KDIAG array of flags for calling diagnostics print routine
      INTEGER, DIMENSION(12) :: KDIAG

!@param NKEYNR number of key number diagnostics
      INTEGER, PARAMETER :: NKEYNR=42
!@param NKEYMO number of months key diagnostics are saved
      INTEGER, PARAMETER :: NKEYMO=50
!@var KEYNR time-series of key numbers
      INTEGER, DIMENSION(NKEYNR,NKEYMO) :: KEYNR

!@var IWRITE,JWRITE,ITWRITE grid point and surface type for diag. output
      INTEGER :: IWRITE = 0, JWRITE = 0, ITWRITE = 0
!@var QCHECK TRUE for running diagnostic checks
      LOGICAL :: QCHECK = .FALSE.

!@var OA generic diagnostic array for ocean heat transport calculations
C****
C****       DATA SAVED IN ORDER TO CALCULATE OCEAN TRANSPORTS
C****
C****       1  ACE1I+SNOWOI  (INSTANTANEOUS AT NOON GMT)
C****       2  TG1OI  (INSTANTANEOUS AT NOON GMT)
C****       3  TG2OI  (INSTANTANEOUS AT NOON GMT)
C****       4  ENRGP  (INTEGRATED OVER THE DAY)
C****       5  SRHDT  (FOR OCEAN, INTEGRATED OVER THE DAY)
C****       6  TRHDT  (FOR OCEAN, INTEGRATED OVER THE DAY)
C****       7  SHDT   (FOR OCEAN, INTEGRATED OVER THE DAY)
C****       8  EVHDT  (FOR OCEAN, INTEGRATED OVER THE DAY)
C****       9  TRHDT  (FOR OCEAN ICE, INTEGRATED OVER THE DAY)
C****      10  SHDT   (FOR OCEAN ICE, INTEGRATED OVER THE DAY)
C****      11  EVHDT  (FOR OCEAN ICE, INTEGRATED OVER THE DAY)
C****      12  SRHDT  (FOR OCEAN ICE, INTEGRATED OVER THE DAY)
C****
      REAL*8, DIMENSION(IM,JM,12) :: OA

C****
C**** Information about acc-arrays:  
C****      names, indices, units, idacc-numbers, etc.

!@var iparm/dparm int/double global parameters written to acc-file
      integer, parameter :: niparm_max=100
      character(len=20), dimension(niparm_max) :: iparm_name
      integer, dimension(niparm_max) :: iparm
      integer :: niparm=0
      integer, parameter :: ndparm_max=100
      character(len=20), dimension(ndparm_max) :: dparm_name
      double precision, dimension(ndparm_max) :: dparm
      integer :: ndparm=0

!@var J_xxx zonal J diagnostic names
      INTEGER :: J_SRINCP0, J_SRNFP0, J_SRNFP1, J_SRABS, J_SRINCG,
     *     J_SRNFG, J_TRNFP0, J_TRNFP1, J_TRHDT, J_RNFP0, J_RNFP1,
     *     J_RHDT, J_SHDT, J_EVHDT, J_F2DT, J_HZ1, J_TG2, J_TG1, J_EVAP,
     *     J_PRCP, J_TX, J_TX1, J_TSRF, J_DTSGST, J_DTDGTR, J_RICST,
     *     J_RICTR, J_ROSST, J_ROSTR, J_RSI, J_TYPE, J_RSNOW, J_SWCOR,
     *     J_OHT,J_OMLT, J_DTDJS, J_DTDJT, J_LSTR, J_LTRO, J_EPRCP,
     *     J_ERUN1,J_EDIFS, J_F1DT, J_ERUN2, J_HZ0, J_DIFS, J_IMELT,
     *     J_RUN2,J_DWTR2, J_WTR1, J_ACE1, J_WTR2, J_ACE2, J_SNOW,
     *     J_RUN1,J_BRTEMP, J_HZ2, J_PCLDSS, J_PCLDMC, J_PCLD, J_CTOPP
     *     ,J_PRCPSS, J_PRCPMC, J_QP, J_GAM, J_GAMM, J_GAMC, J_TRINCG
     *     ,J_FTHERM, J_HSURF, J_HATM, J_PLAVIS, J_PLANIR, J_ALBVIS
     *     ,J_ALBNIR, J_SRRVIS, J_SRRNIR, J_SRAVIS, J_SRANIR, J_CDLDEP
!@var NAME_J,UNITS_J Names/Units of zonal J diagnostics
      character(len=20), dimension(kaj) :: name_j,units_j
!@var LNAME_J Long names of zonal J diagnostics
      character(len=80), dimension(kaj) :: lname_j
      character(len=20), dimension(kaj) :: name_reg
!@var IA_J IDACC indexes for zonal J diagnostics
      integer, dimension(kaj) :: ia_j

      character(len=20), dimension(kapj) :: name_pj,units_pj
      character(len=80), dimension(kapj) :: lname_pj

!@var IJ_xxx AIJ diagnostic names
      INTEGER :: IJ_RSOI, IJ_RSNW, IJ_SNOW, IJ_SHDT, IJ_PREC, IJ_EVAP,
     *     IJ_SSAT, IJ_BETA,  IJ_SLP1,  IJ_P4UV, IJ_PRES, IJ_PHI1K,
     *     IJ_PHI850,IJ_PHI700, IJ_PHI500, IJ_PHI300, IJ_PHI100,
     *     IJ_PHI30,IJ_T850,IJ_PMCCLD, IJ_CLDTPPR, IJ_CLDCV, IJ_PEV,
     *     IJ_TRNFP0, IJ_SRTR,IJ_NETH, IJ_SRNFP0, IJ_SRINCP0, IJ_SRNFG,
     *     IJ_SRINCG, IJ_TG1,IJ_RSIT, IJ_TDSL, IJ_DTDP, IJ_RUNE, IJ_TS1,
     *     IJ_RUNLI, IJ_WS,IJ_TS, IJ_US, IJ_VS, IJ_SLP,IJ_UJET, IJ_VJET,
     *     IJ_PCLDL, IJ_PCLDM, IJ_PCLDH, IJ_BTMPW, IJ_SRREF, IJ_TOC2,
     *     IJ_TAUS, IJ_TAUUS, IJ_TAUVS, IJ_GWTR, IJ_QS, IJ_STRNGTS,
     *     IJ_ARUNU, IJ_DTGDTS, IJ_PUQ, IJ_PVQ, IJ_TGO, IJ_MSI2, IJ_WLM,
     *     IJ_TGO2,IJ_EVAPO, IJ_EVAPI, IJ_EVAPLI, IJ_EVAPE, IJ_F0OC,
     *     IJ_F0OI,IJ_F0LI, IJ_F0E, IJ_F1LI, IJ_SNWF, IJ_TSLI, IJ_ERUN2,
     *     IJ_SHDTLI, IJ_EVHDT, IJ_TRHDT, IJ_TMAX, IJ_TMIN, IJ_TMNMX,
     *     IJ_PEVAP,IJ_TMAXE, IJ_WMSUM, IJ_PSCLD, IJ_PDCLD, IJ_DCNVFRQ,
     *     IJ_SCNVFRQ, IJ_EMTMOM, IJ_SMTMOM, IJ_FPEU, IJ_FPEV, IJ_FMU,
     *     IJ_FMV,IJ_FQU, IJ_FQV, IJ_FGZU, IJ_FGZV, IJ_ERVR, IJ_MRVR,
     *     IJ_SDRAG,IJ_LKON,IJ_LKOFF
!@var IJ_Gxx names for old AIJG arrays (should be more specific!)
      INTEGER :: IJ_G01,IJ_G02,IJ_G03,IJ_G04,IJ_G05,IJ_G06,IJ_G07,
     *     IJ_G08,IJ_G09,IJ_G10,IJ_G11,IJ_G12,IJ_G13,IJ_G14,IJ_G15,
     *     IJ_G16,IJ_G17,IJ_G18,IJ_G19,IJ_G20,IJ_G21,IJ_G22,IJ_G23,
     *     IJ_G24,IJ_G25,IJ_G26,IJ_G27,IJ_G28,IJ_G29
!@var NAME_IJ,UNITS_IJ Names/Units of lat/lon IJ diagnostics
      character(len=20), dimension(kaij) :: name_ij,units_ij
!@var LNAME_IJ Long names of lat/lon IJ diagnostics
      character(len=80), dimension(kaij) :: lname_ij
!@var IA_IJ IDACC indexes for lat/lon IJ diagnostics
      integer, dimension(kaij) :: ia_ij

      character(len=20), dimension(kajl) :: name_jl,units_jl
      character(len=80), dimension(kajl) :: lname_jl

      character(len=20), dimension(kasjl) :: name_sjl,units_sjl
      character(len=80), dimension(kasjl) :: lname_sjl

      character(len=20), dimension(kajk) :: name_jk,units_jk
      character(len=80), dimension(kajk) :: lname_jk

      character(len=20), dimension(kaijk) :: name_ijk,units_ijk
      character(len=80), dimension(kaijk) :: lname_ijk

      character(len=20), dimension(kaijl) :: name_ijl,units_ijl
      character(len=80), dimension(kaijl) :: lname_ijl

      character(len=20), dimension(kwp) :: name_wave,units_wave
      character(len=80), dimension(kwp) :: lname_wave

      character(len=20), dimension(kajlsp) :: name_jlsp,units_jlsp
      character(len=80), dimension(kajlsp) :: lname_jlsp

      character(len=20), dimension(kcon) :: name_consrv,units_consrv
      character(len=80), dimension(kcon) :: lname_consrv

      character(len=20), dimension(kail) :: name_il,units_il
      character(len=80), dimension(kail) :: lname_il

C**** tf_xxx tsfrez diagnostic names
      INTEGER :: tf_day1,tf_last,tf_lkon,tf_lkoff
      character(len=20), dimension(ktsf) :: name_tsf,units_tsf
      character(len=80), dimension(ktsf) :: lname_tsf

      character(len=8), dimension(ntype) :: stype_names=
     &     (/ 'OCEAN   ','OCEANICE','EARTH   ',
     &        'LANDICE ','LAKE    ','LAKEICE ' /)

c idacc-indices of various processes
      integer, parameter ::
     &     ia_src=1, ia_rad=2, ia_srf=3, ia_dga=4, ia_d4a=5, ia_d5f=6,
     *     ia_d5d=7, ia_d5s=8, ia_12hr=9, ia_filt=10, ia_ocn=11,
     *     ia_nmo=12


      END MODULE DAGCOM

      SUBROUTINE io_diags(kunit,it,iaction,ioerr)
!@sum  io_diag reads and writes diagnostics to file
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : ioread,irerun,irsfic,
     *     iowrite,iowrite_mon,iowrite_single
      USE DAGCOM
      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var HEADER Character string label for individual records
      CHARACTER*8 :: HEADER, MODULE_HEADER = "DIAG01"
!@var it input/ouput value of hour
      INTEGER, INTENT(INOUT) :: it

      SELECT CASE (IACTION)
      CASE (IOWRITE)            ! output to standard restart file
        WRITE (kunit,err=10) MODULE_HEADER,KEYNR,TSFREZ,AJ,AREG,APJ,AJL,
     *       ASJL,AIJ,AIL,ENERGY,CONSRV,SPECA,ATPE,ADAILY,WAVE,
     *       AJK,AIJK,AIJL,AJLSP,TDIURN,OA,it
      CASE (IOWRITE_SINGLE)     ! output in single precision
        WRITE (kunit,err=10) MODULE_HEADER,KEYNR,SNGL(TSFREZ),SNGL(AJ),
     *     SNGL(AREG),SNGL(APJ),SNGL(AJL),
     *     SNGL(ASJL),SNGL(AIJ),SNGL(AIL),SNGL(ENERGY),
     *     SNGL(CONSRV),SNGL(SPECA),SNGL(ATPE),SNGL(ADAILY),SNGL(WAVE),
     *     SNGL(AJK),SNGL(AIJK),SNGL(AIJL),SNGL(AJLSP),
     *     it
      CASE (IOWRITE_MON)        ! output to end-of-month restart file
        WRITE (kunit,err=10) MODULE_HEADER,KEYNR,TSFREZ,it
      CASE (ioread)           ! input from restart file
        READ (kunit,err=10) HEADER,KEYNR,TSFREZ,AJ,AREG,APJ,AJL,ASJL,
     *       AIJ,AIL,ENERGY,CONSRV,SPECA,ATPE,ADAILY,WAVE,AJK,
     *       AIJK,AIJL,AJLSP,TDIURN,OA,it
        IF (HEADER.NE.MODULE_HEADER) THEN
          PRINT*,"Discrepancy in module version",HEADER,MODULE_HEADER
          GO TO 10
        END IF
      CASE (Irsfic)      ! no diag-info at beginning of new run needed
      CASE (irerun)      ! only keynr,tsfrez needed at beg of acc-period
        READ (kunit,err=10) HEADER,KEYNR,TSFREZ  ! 'it' is not read in
        IF (HEADER.NE.MODULE_HEADER) THEN
          PRINT*,"Discrepancy in module version",HEADER,MODULE_HEADER
          GO TO 10
        END IF
      END SELECT

      RETURN
 10   IOERR=1
      RETURN
      END SUBROUTINE io_diags
