      MODULE DAGCOM
!@sum  DAGCOM Diagnostic model variables
!@auth Original Development Team
!@var  1.0
      USE MODEL_COM, only : im,jm,lm,imh,ntype
      USE RADNCB, only : LM_REQ

      IMPLICIT NONE
      SAVE
!@var iu_ij,iu_jl,iu_il,iu_j !  units for selected diag. output
      INTEGER iu_ij,iu_jl,iu_il,iu_j,iu_ijk
C**** Accumulating_period information
!ny   INTEGER, PARAMETER :: NSAMPL  !@var NSAMPL # of sampling schemes
!ny   INTEGER, DIMENSION(nsampl) :: IDACC   !@var IDACC acc-counters
      INTEGER, DIMENSION(12) :: MONACC  !@var MONACC(1)=#Januaries, etc
      CHARACTER*12 ACC_PERIOD           !@var string MONyyr1-yyr2
!@var AMON0,JMON0,JDATE0,JYEAR0,JHOUR0,Itime0  beg.of acc-period
!ny   CHARACTER*4 AMON0
!ny   INTEGER JMON0,JDATE0,JYEAR0,JHOUR0,Itime0

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
!@var J50N,J70N special longitudes for AIL diagnostics
      INTEGER :: J50N,J70N

C NEHIST = (TROPO/STRAT)X(ZKE/EKE/SEKE/ZPE/EPE)X(SH/NH)
!@param NEHIST,HISTDAYS number of energy history diagnostics, and days
      INTEGER, PARAMETER :: NEHIST=20
      INTEGER, PARAMETER :: HIST_DAYS=100
!@var ENERGY energy diagnostics
      DOUBLE PRECISION, DIMENSION(NEHIST,HIST_DAYS) :: ENERGY

!@var NPTS number of points at which standard conserv. diags are called
      INTEGER, PARAMETER :: NPTS = 11  ! 9
!@param NQUANT Number of conserved quantities in conservation diags
      INTEGER, PARAMETER :: NQUANT=18  ! 7 20?
!@param KCON number of conservation diagnostics
      INTEGER, PARAMETER :: KCON=125   ! 54
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
!@var IJD6,NAMD6 (i,j)-coord.,4-char names of boxes w/diurnal cycle diag
!ny   INTEGER, DIMENSION(2,NDLYPT) :: IJD6
!ny   CHARACTER*4, DIMENSION(NDLYPT) :: NAMD6
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
      INTEGER, PARAMETER :: NWAV_DAG=min(9,imh)
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

      COMMON /ACCUM/ AJ,AREG,APJ,AJL,ASJL,AIJ,AIL,
     &  ENERGY,CONSRV,SPECA,ATPE,ADAILY,WAVE,
     &  AJK,AIJK,AIJL,AJLSP
      DOUBLE PRECISION, DIMENSION(KACC) :: ACC
      EQUIVALENCE (ACC,AJ)

!@var TSFREZ freezing temperature diagnostics
      integer, parameter :: ktsf=4
      DOUBLE PRECISION, DIMENSION(IM,JM,KTSF) :: TSFREZ

!@param KTD number of diurnal temperature diagnostics
      INTEGER, PARAMETER :: KTD=8
!@var TDIURN diurnal range temperature diagnostics
      DOUBLE PRECISION, DIMENSION(IM,JM,KTD) :: TDIURN

!@nlparams KDIAG array of flags to control diagnostics printout
      INTEGER, DIMENSION(12) :: KDIAG

!@param NKEYNR number of key number diagnostics
      INTEGER, PARAMETER :: NKEYNR=42
!@param NKEYMO number of months key diagnostics are saved
      INTEGER, PARAMETER :: NKEYMO=50
!@var KEYNR time-series of key numbers
      INTEGER, DIMENSION(NKEYNR,NKEYMO) :: KEYNR = 0
!@var KEYCT next index in KEYNR to be used (1->nkeymo)
      INTEGER :: KEYCT = 1

!@nlparams IWRITE,JWRITE,ITWRITE control rad.debug output (i,j,amount)
      INTEGER :: IWRITE = 0, JWRITE = 0, ITWRITE = 0
!@nlparams QCHECK TRUE for running diagnostic checks
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
!@param KOA number of diagnostics needed for ocean heat transp. calcs
      INTEGER, PARAMETER :: KOA = 12
      REAL*8, DIMENSION(IM,JM,KOA) :: OA

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
!@var SCALE_IJ scaling for weighted AIJ diagnostics
      REAL*8, DIMENSION(KAIJ) :: SCALE_IJ
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

!@var IJK_xxx AIJK diagnostic names
      INTEGER :: IJK_U, IJK_V, IJK_DSE, IJK_DP, IJK_T, IJK_Q
!@var SCALE_IJK scaling for weighted AIJK diagnostics
      REAL*8, DIMENSION(KAIJK) :: SCALE_IJK
!@var OFF_IJK offset for weighted AIJK diagnostics
      REAL*8, DIMENSION(KAIJK) :: OFF_IJK

      character(len=20), dimension(kaijk) :: name_ijk,units_ijk
      character(len=80), dimension(kaijk) :: lname_ijk,scname_ijk

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
      USE MODEL_COM, only : ioread,ioread_single,irerun,irsfic,
     *     iowrite,iowrite_mon,iowrite_single,lhead, idacc,nsampl
      USE DAGCOM
      IMPLICIT NONE
      REAL*4 ACCS(KACC),TSFREZS(IM,JM,KTSF)
c???  add a couple of lines to replace ACCS and avoid 'COMMON BLOCK'
      integer monac1(12),i_ida,i_xtra,l_xtra
!@var Kcomb counts acc-files as they are added up
      INTEGER, SAVE :: Kcomb=0

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER idac1(12)
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var HEADER Character string label for individual records
      CHARACTER*80 :: HEADER, MODULE_HEADER = "DIAG01"
!@var it input/ouput value of hour
      INTEGER, INTENT(INOUT) :: it

      write (MODULE_HEADER(LHEAD+1:LHEAD+14),'(a9,i4,a1)')
     *   'R8: keys(',1+NKEYNR*NKEYMO,')'              ! keyct,keynr(:,:)
      i_ida = Lhead+14+13+1
      write (MODULE_HEADER(LHEAD+15:i_ida-1),'(a10,i2,a1)')
     *   ',TSFR(IJM,',KTSF,')'
      write (MODULE_HEADER(i_ida:i_ida+9),'(a7,i2,a1)')
     *   ',idacc(',nsampl,')'
      write (MODULE_HEADER(i_ida+10:i_ida+9+15),'(a5,i9,a1)')
     *   ',acc(',kacc,')'
      i_xtra = i_ida+9+15+1

      SELECT CASE (IACTION)
      CASE (IOWRITE)            ! output to standard restart file
        write (MODULE_HEADER(i_xtra:80),             '(a7,i2,a)')
     *   ',x(IJM,',KTD+KOA,')'  ! make sure that i_xtra+7+2 < 80
        WRITE (kunit,err=10) MODULE_HEADER,keyct,KEYNR,TSFREZ,
     *     idacc,ACC,
c??? *     idacc,AJ,AREG,APJ,AJL,ASJL,AIJ,AIL,ENERGY,CONSRV,SPECA,ATPE,
c??? *     ADAILY,WAVE,AJK,AIJK,AIJL,AJLSP,
     *     TDIURN,OA,it
      CASE (IOWRITE_SINGLE)     ! output in single precision
        MODULE_HEADER(LHEAD+1:LHEAD+2) = 'R4'
        MODULE_HEADER(i_xtra:80) = ',monacc(12)'
        WRITE (kunit,err=10) MODULE_HEADER,keyct,KEYNR,SNGL(TSFREZ),
     *     idacc,SNGL(ACC),
c??? *     idacc,SNGL(AJ),SNGL(AREG),SNGL(APJ),SNGL(AJL),
c??? *     SNGL(ASJL),SNGL(AIJ),SNGL(AIL),SNGL(ENERGY),
c??? *     SNGL(CONSRV),SNGL(SPECA),SNGL(ATPE),SNGL(ADAILY),SNGL(WAVE),
c??? *     SNGL(AJK),SNGL(AIJK),SNGL(AIJL),SNGL(AJLSP),
     *     monacc,it
      CASE (IOWRITE_MON)        ! output to end-of-month restart file
        MODULE_HEADER(i_ida:80) = ',it '
        WRITE (kunit,err=10) MODULE_HEADER,keyct,KEYNR,TSFREZ,it
      CASE (ioread)           ! input from restart file
        READ (kunit,err=10) HEADER,keyct,KEYNR,TSFREZ,
     *      idacc, ACC,
c??? *      idacc, AJ,AREG,APJ,AJL,ASJL,AIJ,AIL,ENERGY,
c??? *      CONSRV,SPECA,ATPE,ADAILY,WAVE,AJK,AIJK,AIJL,AJLSP,
     *      TDIURN,OA,it
        IF (HEADER(1:LHEAD).NE.MODULE_HEADER(1:LHEAD)) THEN
          PRINT*,"Discrepancy in module version",HEADER,MODULE_HEADER
          GO TO 10
        END IF
      CASE (IOREAD_SINGLE)      !
        READ (kunit,err=10) HEADER,keyct,KEYNR,TSFREZS,
     *      idac1,ACCS,
c??? *  add a couple of lines to avoid 'COMMON BLOCK'
     *      monac1,it
!**** Here we could check the dimensions written into HEADER  ??????
        TSFREZ=TSFREZS
        ACC=ACC+ACCS
c??? *  add many lines to avoid 'COMMON BLOCK' - do I really have to ???
        IDACC = IDACC + IDAC1
!@var idacc(5) is the length of a time series (daily energy history).
!****   If combining acc-files, rather than concatenating these series,
!****   we average their beginnings (up to the length of the shortest)
        Kcomb = Kcomb + 1          ! reverse addition, take min instead
        if (Kcomb.gt.1) IDACC(5) = MIN(IDACC(5)-IDAC1(5),IDAC1(5))
        monacc = monacc + monac1
      CASE (irerun)      ! only keynr,tsfrez needed at beg of acc-period
        READ (kunit,err=10) HEADER,keyct,KEYNR,TSFREZ  ! 'it' not read
        IF (HEADER(1:LHEAD).NE.MODULE_HEADER(1:LHEAD)) THEN
          PRINT*,"Discrepancy in module version",HEADER,MODULE_HEADER
          GO TO 10
        END IF
      CASE (Irsfic)      ! no diag-info at beginning of new run needed
      END SELECT

      RETURN
 10   IOERR=1
      RETURN
      END SUBROUTINE io_diags

      SUBROUTINE aPERIOD (JMON1,JYR1,months,years,  aDATE,LDATE)
!@sum  aPERIOD finds a 7 or 12-character name for an accumulation period
!@auth Reto A. Ruedy
!@ver  1.0
      USE MODEL_COM, only : AMONTH
      implicit none
!@var JMON1,JYR1 month,year of beginning of period 1
      INTEGER JMON1,JYR1
!@var JMONM,JMONL middle,last month of period
      INTEGER JMONM,JMONL
!@var months,years length of 1 period,number of periods
      INTEGER months,years
!@var yr1,yr2 (end)year of 1st and last period
      INTEGER yr1,yr2
!@var aDATE date string: MONyyr1(-yyr2)
      character*12 aDATE
!@var LDATE length of date string (7 or 12)
      INTEGER LDATE

      LDATE = 7                  ! if years=1
      if(years.gt.1) LDATE = 12

      aDATE(1:12)=' '
      aDATE(1:3)=AMONTH(JMON1)        ! letters 1-3 of month IF months=1
      yr1=JYR1
      JMONL=JMON1+months-1
      if(JMONL.GT.12) then
         yr1=yr1+1
         JMONL=JMONL-12
      end if
      yr2=yr1+years-1
      write(aDATE(4:7),'(i4.4)') yr1
      if(years.gt.1) write(aDATE(8:12),'(a1,i4.4)') '-',yr2

      if(months.gt.12) aDATE(1:1)='x'                ! should not happen
      if(months.le.1 .or. months.gt.12) return

!**** 1<months<13: adjust characters 1-3 of aDATE (=beg) if necessary:
!**** beg=F?L where F/L=letter 1 of First/Last month for 2-11 mo.periods
!****    =F+L                                        for 2 month periods
!****    =FML where M=letter 1 of Middle month       for 3 month periods
!****    =FnL where n=length of period if n>3         4-11 month periods
      aDATE(3:3)=AMONTH(JMONL)(1:1)            ! we know: months>1
      IF (months.eq.2) aDATE(2:2)='+'
      JMONM = JMONL-1
      IF(JMONM.eq.0) JMONM=12
      IF (months.eq.3) aDATE(2:2)=AMONTH(JMONM)(1:1)
      IF (months.ge.4.and.months.le.9) write (aDATE(2:2),'(I1)') months
      IF (months.eq.10) aDATE(2:2)='X'         ! roman 10
      IF (months.eq.11) aDATE(2:2)='B'         ! hex   11
      IF (months.eq.6) THEN                    !    exceptions:
         IF (JMON1.eq. 5) aDATE(1:3)='NHW'     ! NH warm season May-Oct
         IF (JMON1.eq.11) aDATE(1:3)='NHC'     ! NH cold season Nov-Apr
      END IF
      IF (months.eq.7) THEN                    !    to avoid ambiguity:
         IF (JMON1.eq. 1) aDATE(1:3)='J7L'     ! Jan-Jul J7J->J7L
         IF (JMON1.eq. 7) aDATE(1:3)='L7J'     ! Jul-Jan J7J->L7J
      END IF
      IF (months.eq.12) THEN
C****    beg=ANn where the period ends with month n if n<10 (except 4)
         aDATE(1:3)='ANN'                      ! regular annual mean
         IF (JMONL.le. 9) WRITE(aDATE(3:3),'(I1)') JMONL
         IF (JMONL.eq. 4) aDATE(1:3)='W+C'     ! NH warm+cold seasons
         IF (JMONL.eq.10) aDATE(1:3)='C+W'     ! NH cold+warm seasons
         IF (JMONL.eq.11) aDATE(1:3)='ANM'     ! meteor. annual mean
      END IF
      return
      end SUBROUTINE aPERIOD
