      MODULE DAGCOM
!@sum  DAGCOM Diagnostic model variables
!@auth Original Development Team
!@var  1.0
      USE E001M12_COM, only : im,jm,lm,imh
      USE RADNCB, only : LM_REQ

      IMPLICIT NONE
      
C**** ACCUMULATING DIAGNOSTIC ARRAYS
!@var KAJ number of accumulated zonal budget diagnostics 
      INTEGER, PARAMETER :: KAJ=94   
!@var AJ,BJ,CJ zonal budget diagnostics for each surface type
      DOUBLE PRECISION, DIMENSION(JM,KAJ) :: AJ,BJ,CJ

!@var NREG number of regions for budget diagnostics 
      INTEGER, PARAMETER :: NREG=24
!@var AREG regional budget diagnostics
      DOUBLE PRECISION, DIMENSION(NREG,KAJ) :: AREG
!@var TITREG,NAMREG title and names of regions for AREG diagnostics
      CHARACTER*4 TITREG*80,NAMREG(2,23)
!@var JREH lat/lon array defining regions for AREG diagnostics
      INTEGER, DIMENSION(IM,JM) :: JREG

!@var KAPJ number of zonal pressure diagnostics 
      INTEGER, PARAMETER :: KAPJ=3  ! why isn't this 2?
!@var APJ zonal pressure diagnostics 
      DOUBLE PRECISION, DIMENSION(JM,KAPJ) :: APJ

!@var KAJL number of accumulated AJL diagnostics 
      INTEGER, PARAMETER :: KAJL=57
!@var AJL latitude/height diagnostics
      DOUBLE PRECISION, DIMENSION(JM,LM,KAJL) :: AJL

!@var KASJL number of ASJL diagnostics 
      INTEGER, PARAMETER :: KASJL=4
!@var ASJL latitude/height supplementary diagnostics 
c**** (merge with AJL?)
      DOUBLE PRECISION, DIMENSION(JM,LM_REQ,KASJL) :: ASJL

!@var KAIJ number of AIJ diagnostics 
      INTEGER, PARAMETER :: KAIJ=100
!@var AIJ latitude/longitude diagnostics 
      DOUBLE PRECISION, DIMENSION(IM,JM,KAIJ) :: AIJ
 
!@var KAIL number of AIL diagnostics 
      INTEGER, PARAMETER :: KAIL=16
!@var AIL longitude/height diagnostics 
      DOUBLE PRECISION, DIMENSION(IM,LM,KAIL) :: AIL

!@var KAIJG number of AIJG diagnostics 
      INTEGER, PARAMETER :: KAIJG=29
!@var AIJG latitude/longitude ground diagnostics (merge with AIJ?)
      DOUBLE PRECISION, DIMENSION(IM,JM,KAIJG) :: AIJG

C NEHIST = (TROPO/STRAT)X(ZKE/EKE/SEKE/ZPE/EPE)X(SH/NH)
!@var NEHIST,HISTDAYS number of energy history diagnostics, and days 
      INTEGER, PARAMETER :: NEHIST=20
      INTEGER, PARAMETER :: HIST_DAYS=100
!@var ENERGY energy diagnostics
      DOUBLE PRECISION, DIMENSION(NEHIST,HIST_DAYS) :: ENERGY

!@var KCON number of conservation diagnostics 
      INTEGER, PARAMETER :: KCON=42
!@var CONSRV conservation diagnostics       
      DOUBLE PRECISION, DIMENSION(JM,KCON) :: CONSRV

!@var KSPECA,NSPHER number of spectral diagnostics, and harmonics used 
      INTEGER, PARAMETER :: KSPECA=20
      INTEGER, PARAMETER :: NSPHER=8
!@var SPECA spectral diagnostics
      DOUBLE PRECISION, DIMENSION((IMH+1),KSPECA,NSPHER) :: SPECA

!@var KTPE number of spectral diagnostics for pot. enthalpy
      INTEGER, PARAMETER :: KTPE=8 
      integer, parameter :: NHEMI=2
!@var ATPE pot. enthalpy spectral diagnostics
      DOUBLE PRECISION, DIMENSION(KTPE,NHEMI) :: ATPE

!@var HR_IN_DAY hours in day 
      integer, parameter :: HR_IN_DAY=24
!@var NDLYVAR number of daily diagnostics 
      integer, parameter :: NDLYVAR=63
!@var NDLYPT number of points where daily diagnostics are kept
      integer, parameter :: NDLYPT=4
!@var ADAILY daily diagnostics (24 hour cycles at selected points)
      double precision, dimension(HR_IN_DAY,NDLYVAR,NDLYPT) :: adaily

!@var KAJK number of zonal constant pressure diagnostics 
      INTEGER, PARAMETER :: KAJK=51
!@var AJK zonal constant pressure diagnostics 
      DOUBLE PRECISION, DIMENSION(JM,LM,KAJK) :: AJK

!@var KAIJK number of lat/lon constant pressure diagnostics 
      INTEGER, PARAMETER :: KAIJK=6
!@var KAIJK lat/lon constant pressure diagnostics 
      DOUBLE PRECISION, DIMENSION(IM,JM,LM,KAIJK) :: AIJK

!@var KAIJL number of three dimensional diagnostics 
      INTEGER, PARAMETER :: KAIJL=5
!@var AIJL three dimensional diagnostics 
      DOUBLE PRECISION, DIMENSION(IM,JM,LM,KAIJL) :: AIJL

!@var NWAV_DAG number of components in spectral diagnostics
      integer, parameter :: NWAV_DAG=10
!@var KAJLSP number of spectral diagnostics
      INTEGER, PARAMETER :: KAJLSP=3
!@var AJLSP spectral diagnostics
      DOUBLE PRECISION, DIMENSION(JM,LM,NWAV_DAG,KAJLSP) :: AJLSP

!@var N12HRS_IN_31DAY number of frequency diagnostics
      integer, parameter :: N12HRS_IN_31DAY=62
!@var RE_AND_IM complex components of frequency diagnostics
      integer, parameter :: RE_AND_IM=2
!@var KWP number of frequency diagnostics
      integer, parameter :: KWP=12
!@var WAVE frequency diagnostics
      DOUBLE PRECISION,
     &     DIMENSION(RE_AND_IM,N12HRS_IN_31DAY,NWAV_DAG,KWP) :: WAVE

!@var KACC total number of diagnostic elements in ACCUM array
      INTEGER, PARAMETER :: KACC= JM*KAJ + JM*KAJ + JM*KAJ + NREG*KAJ
     *     + JM*KAPJ + JM*LM*KAJL + JM*LM_REQ*KASJL + IM*JM*KAIJ +
     *     IM*LM*KAIL + IM*JM*KAIJG + NEHIST*HIST_DAYS + JM*KCON +
     *     (IMH+1)*KSPECA*NSPHER + KTPE*NHEMI + HR_IN_DAY*NDLYVAR*NDLYPT
     *     + RE_AND_IM*N12HRS_IN_31DAY*NWAV_DAG*KWP + JM*LM*KAJK +
     *     IM*JM*LM*KAIJK + IM*JM*LM*KAIJL + JM*LM*NWAV_DAG*KAJLSP

      COMMON /ACCUM/ AJ,BJ,CJ,AREG,APJ,AJL,ASJL,AIJ,AIL,
     &  AIJG,ENERGY,CONSRV,SPECA,ATPE,ADAILY,WAVE,
     &  AJK,AIJK,AIJL,AJLSP

!@var TSFREZ freezing temperature diagnostics
      DOUBLE PRECISION, DIMENSION(IM,JM,2) :: TSFREZ

!@var KTD number of diurnal temperature diagnostics 
      INTEGER, PARAMETER :: KTD=8
!@var TDIURN diurnal range temperature diagnostics 
      DOUBLE PRECISION, DIMENSION(IM,JM,KTD) :: TDIURN

!@var KDIAG array of flags for calling diagnostics print routine
      INTEGER, DIMENSION(12) :: KDIAG

!@var NKEYNR number of key number diagnostics
      INTEGER, PARAMETER :: NKEYNR=42
!@var NKEYMO number of months key diagnostics are saved
      INTEGER, PARAMETER :: NKEYMO=50
!@var KEYNR time-series of key numbers  
      INTEGER, DIMENSION(NKEYNR,NKEYMO) :: KEYNR

!@var IWRITE,JWRITE,ITWRITE grid point and surface type for diag. output
      INTEGER :: IWRITE = 0, JWRITE = 0, ITWRITE = 0
!@var QCHECK TRUE for running diagnostic checks
      LOGICAL :: QCHECK = .FALSE.

!@var IJ_xxx AIJ diagnostic names
      INTEGER :: IJ_RSOI, IJ_RSNW, IJ_SNOW, IJ_SHDT, IJ_PREC, IJ_EVAP,
     *     IJ_SSAT, IJ_BETA,  IJ_SLP1,  IJ_P4UV, IJ_PRES, IJ_PHI1K,
     *     IJ_PHI850,IJ_PHI700, IJ_PHI500, IJ_PHI300, IJ_PHI100,
     *     IJ_PHI30,IJ_T850,IJ_PMCCLD, IJ_CLDTPPR, IJ_CLDCV, IJ_PEV,
     *     IJ_TRNFP0, IJ_SRTR,IJ_NETH, IJ_SRNFP0, IJ_SRINCP0, IJ_SRNFG,
     *     IJ_SRINCG, IJ_TG1,IJ_RSIT, IJ_TDSL, IJ_DTDP, IJ_RUNE, IJ_TS1,
     *     IJ_RUNLI, IJ_WS,IJ_TS, IJ_US, IJ_VS, IJ_SLP,IJ_UJET, IJ_VJET,
     *     IJ_PCLDL, IJ_PCLDM, IJ_PCLDH, IJ_BTMPW, IJ_SRREF, IJ_ODATA4,
     *     IJ_TAUS, IJ_TAUUS, IJ_TAUVS, IJ_GWTR, IJ_QS, IJ_STRNGTS,
     *     IJ_ARUNU, IJ_DTGDTS, IJ_PUQ, IJ_PVQ, IJ_TGO, IJ_MSI2, IJ_WLM,
     *     IJ_TGO2,IJ_EVAPO, IJ_EVAPI, IJ_EVAPLI, IJ_EVAPE, IJ_F0OC,
     *     IJ_F0OI,IJ_F0LI, IJ_F0E, IJ_F1LI, IJ_SNWF, IJ_TSLI, IJ_ERUN2,
     *     IJ_SHDTLI, IJ_EVHDT, IJ_TRHDT, IJ_TMAX, IJ_TMIN, IJ_TMNMX,
     *     IJ_PEVAP,IJ_TMAXE, IJ_WMSUM, IJ_PSCLD, IJ_PDCLD, IJ_DCNVFRQ,
     *     IJ_SCNVFRQ, IJ_EMTMOM, IJ_SMTMOM, IJ_FPEU, IJ_FPEV, IJ_FMU,
     *     IJ_FMV,IJ_FQU, IJ_FQV, IJ_FGZU, IJ_FGZV, IJ_ERVR, IJ_MRVR
!@var NAME_IJ Names of IJ diagnostics
      CHARACTER*7 NAME_IJ(KAIJ)
!@var IA_IJ IDACC indexes for IJ diagnostics
      INTEGER IA_IJ(KAIJ)

      END MODULE DAGCOM
