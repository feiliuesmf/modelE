#include "rundeck_opts.h"
      MODULE DIAG_COM
!@sum  DIAG_COM Diagnostic model variables
!@auth Original Development Team
!@ver  1.0
      USE CONSTANT, only : twopi
      USE MODEL_COM, only : im,jm,lm,imh,fim,ntype,kep,istrat,lm_req
      USE GEOM, only : dlon

      IMPLICIT NONE
      SAVE
      private

      public LM_REQ,im,jm,lm,imh,ntype,istrat,kep

C**** Accumulating_period information
      INTEGER, DIMENSION(12), public :: MONACC  !@var MONACC(1)=#Januaries, etc
      CHARACTER*12, public :: ACC_PERIOD='PARTIAL'    !@var string MONyyr1-yyr2
!@var AMON0,JMON0,JDATE0,JYEAR0,JHOUR0,Itime0  beg.of acc-period

!!  WARNING: if new diagnostics are added, change io_diags/reset_DIAG !!
C**** ACCUMULATING DIAGNOSTIC ARRAYS

C**** Define a budget grid (for AJ,AREG,CONSRV diags)
C**** This is fixed for all grid types and resolutions
!@param JM_BUDG grid size for budget page diags
! will be 46 based on basic 4x5 grid
      INTEGER, PARAMETER, public :: JM_BUDG=JM  ! for temporary continuity
c!@var DXYP_BUDG area array on budget grid
c      REAL*8, DIMENSION(JM_BUDG) :: DXYP_BUDG


!@param KAJ number of accumulated zonal budget diagnostics
      INTEGER, PARAMETER, public :: KAJ=80
!@var AJ zonal budget diagnostics for each surface type
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:), public :: AJ,AJ_loc

!@var SQRTM moved from DIAG5A where it was a saved local array to this
!@var place so its size could be allocated dynamically and still have
!@var it preserved from call to call of DIAG5A
      REAL*8, ALLOCATABLE, DIMENSION(:,:), public :: SQRTM
!@param NREG number of regions for budget diagnostics
      INTEGER, PARAMETER, public :: NREG=24
!@var AREGJ regional budget diagnostics
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:), public :: AREGJ,AREGJ_loc

!@var TITREG,NAMREG title and names of regions for AREG diagnostics
      CHARACTER*4, public :: TITREG*80,NAMREG(2,23)
!@var JREH lat/lon array defining regions for AREG diagnostics
cgsfc      INTEGER, ALLOCATABLE, DIMENSION(:,:), public :: JREG
      INTEGER, DIMENSION(IM,JM), public :: JREG

!@param KAPJ number of zonal pressure diagnostics
      INTEGER, PARAMETER, public :: KAPJ=2
!@var APJ zonal pressure diagnostics
      REAL*8, ALLOCATABLE, DIMENSION(:,:), public :: APJ,APJ_loc

!@param KAJL,KAJLX number of AJL diagnostics,KAJLX includes composites
      INTEGER, PARAMETER, public :: KAJL=70+KEP, KAJLX=KAJL+50
!@var AJL latitude/height diagnostics
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:), public :: AJL,AJL_loc

!@param KASJL number of ASJL diagnostics
      INTEGER, PARAMETER, public :: KASJL=4
!@var ASJL latitude/height supplementary diagnostics (merge with AJL?)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:), public :: ASJL,ASJL_loc

!@param KAIJ,KAIJX number of AIJ diagnostics, KAIJX includes composites
      INTEGER, PARAMETER, public :: KAIJ=328 , KAIJX=KAIJ+400
!@var AIJ latitude/longitude diagnostics
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:), public :: AIJ,AIJ_loc

!@param KAIL number of AIL diagnostics
      INTEGER, PARAMETER, public :: KAIL=15
!@var AIL longitude/height diagnostics
      REAL*8, DIMENSION(IM,LM,KAIL), public :: AIL
!@var J50N,J70N,J5NUV,J5SUV,J5S,J5N special latitudes for AIL diags
      INTEGER, PARAMETER, public :: J50N  = (50.+90.)*(JM-1)/180.+1.5
      INTEGER, PARAMETER, public :: J70N  = (70.+90.)*(JM-1)/180.+1.5
      INTEGER, PARAMETER, public :: J5NUV = (90.+5.)*(JM-1.)/180.+2.
      INTEGER, PARAMETER, public :: J5SUV = (90.-5.)*(JM-1.)/180.+2.
      INTEGER, PARAMETER, public :: J5N   = (90.+5.)*(JM-1.)/180.+1.5
      INTEGER, PARAMETER, public :: J5S   = (90.-5.)*(JM-1.)/180.+1.5

C NEHIST=(TROPO/L STRAT/M STRAT/U STRAT)X(ZKE/EKE/SEKE/ZPE/EPE)X(SH/NH)
!@param NED number of different energy history diagnostics
!@param NEHIST,HIST_DAYS number of energy history columns,rows (max)
      INTEGER, PARAMETER, public :: NED=10
      INTEGER, PARAMETER, public :: NEHIST=NED*(2+ISTRAT)
      INTEGER, PARAMETER, public :: HIST_DAYS=100
!@var ENERGY energy diagnostics
      REAL*8, DIMENSION(NEHIST,HIST_DAYS), public :: ENERGY

!@var NPTS number of points at which standard conserv. diags are called
      INTEGER, PARAMETER, public :: NPTS = 11
!@param NQUANT Number of conserved quantities in conservation diags
      INTEGER, PARAMETER, public :: NQUANT=22
!@param KCON number of conservation diagnostics
      INTEGER, PARAMETER, public :: KCON=170
!@var CONSRV conservation diagnostics
      REAL*8, ALLOCATABLE, DIMENSION(:,:), public :: CONSRV,CONSRV_loc
!@var SCALE_CON scales for conservation diagnostics
      REAL*8, DIMENSION(KCON), public :: SCALE_CON
!@var TITLE_CON titles for conservation diagnostics
      CHARACTER*32, DIMENSION(KCON), public :: TITLE_CON
!@var NSUM_CON indices for summation of conservation diagnostics
!@var IA_CON IDACC numbers for conservation diagnostics
      INTEGER, DIMENSION(KCON), public :: NSUM_CON, IA_CON
!@var NOFM indices for CONSRV array
      INTEGER, DIMENSION(NPTS+1,NQUANT), public :: NOFM
!@var icon_xx indexes for conservation quantities
      INTEGER, public ::
     &     icon_AM,icon_KE,icon_MS,icon_TPE,icon_WM,icon_LKM
     *     ,icon_LKE,icon_EWM,icon_WTG,icon_HTG,icon_OCE,icon_OMSI
     *     ,icon_OHSI,icon_OSSI,icon_LMSI,icon_LHSI,icon_MLI,icon_HLI
!@var KCMX actual number of conservation diagnostics
      INTEGER, public :: KCMX = 25 ! take up first 25 indexes for special cases
!@var CONPT0 default titles for each point where conserv diags. are done
      CHARACTER*10, DIMENSION(NPTS), public :: CONPT0 = (/
     *     "DYNAMICS  ","CONDENSATN","RADIATION ","PRECIPITAT",
     *     "LAND SURFC","SURFACE   ","FILTER    ","OCEAN     ",
     *     "DAILY     ","SRF OCN FL","OCN DYNAM "/)

!@param KSPECA,NSPHER number of spectral diagnostics, and harmonics used
      INTEGER, PARAMETER, public :: KSPECA=20
      INTEGER, PARAMETER, public :: NSPHER=4*(2+ISTRAT)
!@var SPECA spectral diagnostics
      REAL*8, DIMENSION((IMH+1),KSPECA,NSPHER), public :: SPECA
!@var KLAYER index for dividing up atmosphere into layers for spec.anal.
      INTEGER, DIMENSION(LM), public :: KLAYER
!@param PSPEC pressure levels at which layers are seperated and defined
C**** 1000 - 150: troposphere           150 - 10 : low strat.
C****   10 - 1: mid strat               1 and up : upp strat.
      REAL*8, DIMENSION(4), PARAMETER, public ::
     &     PSPEC = (/ 150., 10., 1., 0. /)
!@var LSTR level of interface between low and mid strat. (approx 10 mb)
      INTEGER, public :: LSTR = LM   ! defaults to model top.

!@param KTPE number of spectral diagnostics for pot. enthalpy
      INTEGER, PARAMETER, public :: KTPE=8
      integer, parameter, public :: NHEMI=2
!@var ATPE pot. enthalpy spectral diagnostics
      REAL*8, DIMENSION(KTPE,NHEMI), public :: ATPE

!@param HR_IN_DAY hours in day
      INTEGER, PARAMETER, public :: HR_IN_DAY=24
!@param NDIUVAR number of diurnal diagnostics
#ifdef TRACERS_DUST
      INTEGER, PARAMETER, public :: NDIUVAR=302
#else
#if (defined TRACERS_MINERALS) || (defined TRACERS_QUARZHEM)
      INTEGER, PARAMETER, public :: NDIUVAR=61
#else
      INTEGER, PARAMETER, public :: NDIUVAR=56
#endif
#endif
!@param NDIUPT number of points where diurnal diagnostics are kept
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
      INTEGER, PARAMETER, public :: NDIUPT=34
#else
      INTEGER, PARAMETER, public :: NDIUPT=4
#endif
!@dpparam adiurn_dust  flag to switch on/off intra daily diagnostics for dust
!@+                    default=0 (off)
      INTEGER, public :: adiurn_dust=0
!@dbparam IJDD,NAMDD (i,j)-coord.,names of boxes w/diurnal cycle diag
!@+       defaults set in DIAG_RES (depends on resolution)
      INTEGER, DIMENSION(2,NDIUPT), public :: IJDD
      CHARACTER*4, DIMENSION(NDIUPT), public :: NAMDD

!@var ADIURN diurnal diagnostics (24 hour cycles at selected points)
      REAL*8, DIMENSION(HR_IN_DAY,NDIUVAR,NDIUPT), public :: ADIURN
!@param HR_IN_MONTH hours in month
      INTEGER, PARAMETER, public :: HR_IN_MONTH=HR_IN_DAY*31
#ifndef NO_HDIURN
!@var HDIURN hourly diagnostics (hourly value at selected points)
!@+     Same quantities as ADIURN but not averaged over the month
      REAL*8, DIMENSION(HR_IN_MONTH,NDIUVAR,NDIUPT), public :: HDIURN
#endif
!@param KAJK number of zonal constant pressure diagnostics
!@param KAJKX number of zonal constant pressure composit diagnostics
      INTEGER, PARAMETER, public :: KAJK=51, KAJKX=KAJK+100
!@var AJK zonal constant pressure diagnostics
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:), public :: AJK,AJK_loc

!@param KAIJK,KAIJX number of lat/lon constant pressure diagnostics
      INTEGER, PARAMETER, public :: KAIJK=28, kaijkx=kaijk+400
!@var KAIJK lat/lon constant pressure diagnostics
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:), public :: AIJK,AIJK_loc

!@param NWAV_DAG number of components in spectral diagnostics
      INTEGER, PARAMETER, public :: NWAV_DAG=min(9,imh)
!@param Max12HR_sequ,Min12HR_sequ lengths of time series for wave powers
      INTEGER, PARAMETER, public :: Max12HR_sequ=2*31, Min12HR_sequ=2*28
!@param RE_AND_IM complex components of wave power diagnostics
      INTEGER, PARAMETER, public :: RE_AND_IM=2
!@param KWP number of wave power diagnostics
      INTEGER, PARAMETER, public :: KWP=12
!@var WAVE frequency diagnostics (wave power)
      REAL*8, DIMENSION(RE_AND_IM,Max12HR_sequ,NWAV_DAG,KWP), public ::
     &     WAVE

C**** parameters and variables for ISCCP diags
!@param ntau,npress number of ISCCP optical depth,pressure categories
      integer, parameter, public :: ntau=7,npres=7
!@param nisccp number of ISCCP histogram regions
      integer, parameter, public :: nisccp = 5
!@var isccp_press pressure mid points for isccp histogram
      INTEGER, PARAMETER, public ::
     &     isccp_press(npres) = (/ 90, 245, 375, 500,
     *     630, 740, 900 /)
!@var isccp_tau lower bound of optical depth for each isccp tau category
      REAL*8, PARAMETER, public ::
     &     isccp_tau(ntau) = (/ 0d0,.1d0,1.3d0,3.6d0,
     *     9.4d0,23d0,60d0 /)
!@var isccp_taum mid point of optical depth for each isccp tau category
      REAL*8, PARAMETER, public :: isccp_taum(ntau-1) = (/ .05d0,0.7d0,2
     *     .95d0,6.5d0,16.2d0,41.5d0 /)
!@var isccp_late edge latitudes for each isccp lat category (region)
!@var isccp_lat midpoint latitudes for each isccp lat category (region)
! calculation of midpoints is hard-coded until all fortran compilers
! allow array arithmetic in PARAMETERS
      REAL*8, PARAMETER, public ::
     &  isccp_late(nisccp+1)=(/-60,-30,-15,15,30,60/)
     & ,isccp_lat(nisccp)=(/-45.,-22.5,0.,22.5,45./)
!@var isccp_reg latitudinal index for ISCCP histogram regions
      integer, public :: isccp_reg(JM)
!@var AISCCP accumlated array of ISCCP histogram
      real*8, public :: AISCCP(ntau,npres,nisccp)

!@param KGZ number of pressure levels for some diags
      INTEGER, PARAMETER, public :: KGZ = 13
!@param kgz_max is the actual number of geopotential heights saved
      INTEGER, public :: kgz_max
!@param PMB pressure levels for geopotential heights (extends to strat)
!@param GHT ~mean geopotential heights at PMB level (extends to strat)
!@param PMNAME strings describing PMB pressure levels
      REAL*8, DIMENSION(KGZ), PARAMETER, public ::
     &     PMB=(/1000d0,850d0,700d0,500d0,300d0,100d0,30d0,10d0,
     *     3.4d0,.7d0,.16d0,.07d0,.03d0/),
     *     GHT=(/0.,1500.,3000.,5600.,9500.,16400.,24000.,30000.,
     *     40000.,50000.,61000.,67000.,72000./)
      CHARACTER*4, DIMENSION(KGZ), PARAMETER, public :: PMNAME=(/
     *     "1000","850 ","700 ","500 ","300 ","100 ","30  ","10  ",
     *     "3.4 ","0.7 ",".16 ",".07 ",".03 " /)

C**** Instantaneous constant pressure level fields
!@var Z_inst saved instantaneous height field (at PMB levels)
!@var RH_inst saved instantaneous relative hum (at PMB levels)
!@var T_inst saved instantaneous temperature(at PMB levels)
!@var P_acc accumulated precip (special for SUBDD)
!@var PM_acc accumulated moist convective precip (special for SUBDD)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:), public ::
     &     Z_inst,RH_inst,T_inst
      REAL*8, ALLOCATABLE, DIMENSION(:,:), public :: P_acc,PM_acc

      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:), public :: AFLX_ST

!@param KTSF number of freezing temperature diagnostics
      integer, parameter, public :: ktsf=4
!@var TSFREZ freezing temperature diagnostics
C****   1  FIRST DAY OF GROWING SEASON (JULIAN DAY)
C****   2  LAST DAY OF GROWING SEASON (JULIAN DAY)
C****   3  LAST DAY OF ICE-FREE LAKE (JULIAN DAY)
C****   4  LAST DAY OF ICED-UP LAKE  (JULIAN DAY)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:), public :: TSFREZ,TSFREZ_loc

!@param KTD number of diurnal temperature diagnostics
      INTEGER, PARAMETER, public :: KTD=9
!@var TDIURN diurnal range temperature diagnostics
C****   1  MIN TG1 OVER EARTH FOR CURRENT DAY (C)
C****   2  MAX TG1 OVER EARTH FOR CURRENT DAY (C)
C****   3  MIN TS OVER EARTH FOR CURRENT DAY (K)
C****   4  MAX TS OVER EARTH FOR CURRENT DAY (K)
C****   5  SUM OF COMPOSITE TS OVER TIME FOR CURRENT DAY (C)
C****   6  MAX COMPOSITE TS FOR CURRENT DAY (K)
C****   7  MAX TG1 OVER OCEAN ICE FOR CURRENT DAY (C)
C****   8  MAX TG1 OVER LAND ICE FOR CURRENT DAY (C)
C****   9  MIN COMPOSITE TS FOR CURRENT DAY (K)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:), public :: TDIURN

!@nlparam KDIAG array of flags to control diagnostics printout
      INTEGER, DIMENSION(13), public :: KDIAG

!@param NKEYNR number of key number diagnostics
      INTEGER, PARAMETER, public :: NKEYNR=42
!@param NKEYMO number of months key diagnostics are saved
      INTEGER, PARAMETER, public :: NKEYMO=50
!@var KEYNR time-series of key numbers
      INTEGER, DIMENSION(NKEYNR,NKEYMO), public :: KEYNR = 0
!@var KEYCT next index in KEYNR to be used (1->nkeymo)
      INTEGER, public :: KEYCT = 1

!@nlparam IWRITE,JWRITE,ITWRITE control rad.debug output (i,j,amount)
      INTEGER, public :: IWRITE = 0, JWRITE = 0, ITWRITE = 23
!@nlparam QDIAG TRUE for outputting binary diagnostics
      LOGICAL, public :: QDIAG = .FALSE.
!@nlparam QDIAG_ratios TRUE for forming ratios if title="q1 x q2"
      LOGICAL, public :: QDIAG_ratios = .TRUE.

!@var OA generic diagnostic array for ocean heat transport calculations
C****
C****       DATA SAVED IN ORDER TO CALCULATE OCEAN TRANSPORTS
C****
C****       1  ACE1I+SNOWOI  (INSTANTANEOUS AT NOON GMT)
C****       2  MSI2   (INSTANTANEOUS AT NOON GMT)
C****       3  HSIT   (INSTANTANEOUS AT NOON GMT)
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
C**** Extra array needed for dealing with advected ice
C****      13  HCHSI  (HORIZ CONV SEA ICE ENRG, INTEGRATED OVER THE DAY)
C****
!@param KOA number of diagnostics needed for ocean heat transp. calcs
      INTEGER, PARAMETER, public :: KOA = 13
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:), public :: OA

C****
C**** Information about acc-arrays:
C****      names, indices, units, idacc-numbers, etc.

!@var iparm/dparm int/double global parameters written to acc-file
      integer, parameter, public :: niparm_max=100
      character(len=20), dimension(niparm_max), public :: iparm_name
      integer, dimension(niparm_max), public :: iparm
      integer, public :: niparm=0
      integer, parameter, public :: ndparm_max=100
      character(len=20), dimension(ndparm_max), public :: dparm_name
      REAL*8, dimension(ndparm_max), public :: dparm
      integer, public :: ndparm=0

!@var J_xxx zonal J diagnostic names
      INTEGER, public ::
     &     J_SRINCP0, J_SRNFP0, J_SRNFP1, J_SRABS, J_SRINCG,
     *     J_SRNFG, J_TRNFP0, J_TRNFP1, J_TRHDT, J_RNFP0, J_RNFP1,
     *     J_RHDT, J_SHDT, J_EVHDT, J_HZ1, J_TG2, J_TG1, J_EVAP,
     *     J_PRCP, J_TX, J_TX1, J_TSRF, J_DTSGST, J_DTDGTR, J_RICST,
     *     J_RICTR, J_ROSST, J_ROSTR, J_RSI, J_TYPE, J_RSNOW,
     *     J_OHT, J_DTDJS, J_DTDJT, J_LSTR, J_LTRO, J_EPRCP,
     *     J_RUN, J_ERUN, J_HZ0, J_H2OCH4,
     *     J_RVRD,J_ERVR,J_IMELT, J_HMELT, J_SMELT,J_IMPLM, J_IMPLH,
     *     J_WTR1,J_ACE1, J_WTR2,J_ACE2, J_SNOW, J_BRTEMP, J_HZ2,
     *     J_PCLDSS,J_PCLDMC, J_PCLD,J_CTOPP, J_PRCPSS, J_PRCPMC, J_QP,
     *     J_GAM,J_GAMM, J_GAMC,J_TRINCG, J_FTHERM, J_HSURF, J_HATM,
     *     J_PLAVIS,J_PLANIR,J_ALBVIS, J_ALBNIR, J_SRRVIS, J_SRRNIR,
     *     J_SRAVIS,J_SRANIR,J_CLDDEP, J_CLRTOA, J_CLRTRP, J_TOTTRP
!@var NAME_J,UNITS_J Names/Units of zonal J diagnostics
      character(len=20), dimension(kaj), public :: name_j,units_j
!@var LNAME_J Long names of zonal J diagnostics
      character(len=80), dimension(kaj), public :: lname_j
!@var STITLE_J short titles for print out for zonal J diagnostics
      character(len=16), dimension(kaj), public :: stitle_j
!@var SCALE_J scale for zonal J diagnostics
      real*8, dimension(kaj), public :: scale_j
!@var IA_J IDACC indexes for zonal J diagnostics
      integer, dimension(kaj), public :: ia_j
!@var k_j_out number of directly printed out budget diags
      integer, public :: k_j_out

      character(len=20), dimension(kaj), public :: name_reg
      character(len=20), dimension(kapj), public :: name_pj,units_pj
      character(len=80), dimension(kapj), public :: lname_pj

!@var IJ_xxx AIJ diagnostic names
      INTEGER, public ::
     &     IJ_RSOI, IJ_RSNW, IJ_SNOW, IJ_SHDT, IJ_PREC, IJ_EVAP,
     *     IJ_SSAT, IJ_BETA,  IJ_SLP1,  IJ_P4UV, IJ_PRES, IJ_PHI1K,
     *     IJ_PHI850, IJ_PHI700, IJ_PHI500, IJ_PHI300, IJ_PHI100,
     *     IJ_PHI30, IJ_PHI10, IJ_PHI3p4, IJ_PHI0p7, IJ_PHI0p16,
     *     IJ_PHI0p07, IJ_PHI0p03, IJ_T850, IJ_T500, IJ_T300, IJ_Q850,
     *     IJ_Q500, IJ_Q300, IJ_PMCCLD, IJ_CLDTPPR, IJ_CLDCV, IJ_DSEV,
     *     IJ_CLDTPT, IJ_CLDCV1, IJ_CLDT1T,IJ_CLDT1P, IJ_T700, IJ_Q700,
     *     ij_wtrcld,ij_icecld,ij_optdw,ij_optdi,IJ_PBLHT, IJ_RH700,
     *     IJ_RH1, IJ_RH850, IJ_RH500, IJ_RH300, IJ_SWCRF, IJ_LWCRF,
     *     IJ_TRNFP0, IJ_SRTR, IJ_NETH, IJ_SRNFP0, IJ_SRINCP0, IJ_SRNFG,
     *     IJ_SRINCG, IJ_TG1, IJ_RSIT, IJ_TDSL, IJ_TDCOMP, IJ_DTDP,
     *     IJ_RUNE, IJ_TS1, IJ_RUNLI, IJ_WS, IJ_TS, IJ_US, IJ_VS,
     *     IJ_SLP, IJ_UJET, IJ_VJET, IJ_PCLDL, IJ_PCLDM, IJ_PCLDH,
     *     IJ_BTMPW, IJ_SRREF, IJ_SRVIS, IJ_TOC2, IJ_TAUS, IJ_TAUUS,
     *     IJ_TAUVS, IJ_GWTR, IJ_QS, IJ_STRNGTS, IJ_ARUNU, IJ_DTGDTS,
     *     IJ_PUQ, IJ_PVQ, IJ_TGO, IJ_MSI, IJ_TGO2, IJ_EVAPO,
     *     IJ_EVAPI, IJ_EVAPLI,IJ_EVAPE, IJ_F0OC,IJ_F0OI,IJ_F0LI,IJ_F0E,
     *     IJ_F1LI, IJ_SNWF, IJ_TSLI, IJ_ERUN2, IJ_SHDTLI, IJ_EVHDT,
     *     IJ_TRHDT, IJ_TMAX, IJ_TMIN, IJ_TMNMX, IJ_PEVAP, IJ_TMAXE,
     *     IJ_WMSUM, IJ_PSCLD, IJ_PDCLD, IJ_DCNVFRQ, IJ_SCNVFRQ,
     *     IJ_EMTMOM, IJ_SMTMOM, IJ_FMU, IJ_FMV, IJ_SSTABX,
     *     IJ_FGZU, IJ_FGZV, IJ_ERVR, IJ_MRVR, IJ_SSS, IJ_PRECMC,
     *     IJ_LKON, IJ_LKOFF, IJ_LKICE, IJ_PTROP, IJ_TTROP, IJ_TSI,
     *     IJ_SSI1,IJ_SSI2,IJ_SMFX,     ! IJ_MSU2,IJ_MSU2R,
     *     IJ_MLTP,IJ_FRMP, IJ_P850, IJ_CLR_SRINCG,
     *     IJ_GPP, IJ_RAUTO, IJ_CLAB, IJ_DLEAF, !VEG DIAGNOSTICS
     *     IJ_SOILRESP, IJ_SOILCPOOLSUM, !additional veg diags (soil bgc)
     *     IJ_GICE, IJ_GWTR1, IJ_ZSNOW, IJ_AFLMLT, IJ_AERUNS, IJ_AERUNU,
     *     IJ_HTSOIL, IJ_HTSNOW, IJ_AINTRCP,
     *     IJ_SRNTP,IJ_TRNTP,IJ_CLR_SRNTP,IJ_CLR_TRNTP, IJ_TRSDN,
     *     IJ_TRSUP, IJ_CLR_SRNFG,IJ_CLR_TRDNG,IJ_CLR_SRUPTOA,
     *     IJ_CLR_TRUPTOA, IJ_CLDW, IJ_CLDI, IJ_QM, IJ_SSH, IJ_FWOC,
     *     IJ_FWIO, IJ_HTIO, IJ_STIO, IJ_DSKIN,
     *     IJ_SWDCLS,IJ_SWNCLS,IJ_LWDCLS,IJ_SWNCLT,IJ_LWNCLT,
     *     IJ_P1000,IJ_P925,IJ_P700,IJ_P600,IJ_P500, IJ_LI, IJ_LK,
     &     IJ_FVEG,IJ_GUSTI, IJ_MCCON, IJ_SRVDIR, IJ_SRVISSURF
     *     ,IJ_WISUM
     *     ,ij_3dnwm,ij_3dnim,ij_3dnws,ij_3dnis
     *     ,ij_3drwm,ij_3drim,ij_3drws,ij_3dris
     *     ,ij_3dlwm,ij_3dlim,ij_3dlws,ij_3dlis
     *     ,ij_ssprec,ij_mcprec,IJ_WMCLWP,IJ_WMCTWP
     &     ,ij_wdry,ij_wtke,ij_wmoist,ij_wsgcm,ij_wspdf
     &     ,ij_flam
      INTEGER, public, dimension(LM) :: IJ_MCamFX,IJ_TEMPL,IJ_GRIDH
     &     ,IJ_HUSL
!@var IJ_Gxx names for old AIJG arrays (should be more specific!)
      INTEGER, public ::
     &   IJ_G01,IJ_G02,IJ_G03,IJ_G04,IJ_G05,IJ_G06,IJ_G07,
     *     IJ_G08,IJ_G09,IJ_G10,IJ_G11,IJ_G12,IJ_G13,IJ_G14,IJ_G15,
     *     IJ_G16,IJ_G17,IJ_G18,IJ_G19,IJ_G20,IJ_G21,IJ_G22,IJ_G23,
     *     IJ_G24,IJ_G25,IJ_G26,IJ_G27,IJ_G28,IJ_G29
!@var IJ_GWx names for gravity wave diagnostics
      INTEGER, public ::
     &     IJ_GW1,IJ_GW2,IJ_GW3,IJ_GW4,IJ_GW5,IJ_GW6,IJ_GW7,IJ_GW8
     *     ,IJ_GW9
!@var IJ_xxxI names for ISCCP diagnostics
      INTEGER, public ::
     &     IJ_CTPI,IJ_TAUI,IJ_LCLDI,IJ_MCLDI,IJ_HCLDI,IJ_TCLDI,IJ_SCLDI

!@param LEGEND "contour levels" for ij-maps
      CHARACTER(LEN=40), DIMENSION(25), PARAMETER, public :: LEGEND=(/ !
     1  '0=0,1=5...9=45,A=50...K=100             ', ! ir_pct    fac=.2
     2  '0=0...9=90,A=100...I=180...R=270        ', ! ir_angl       .1
     3  '1=.5...9=4.5,A=5...Z=17.5,+=MORE        ', ! ir_0_18        2
     4  '1=.1...9=.9,A=1...Z=3.5,+=MORE          ', ! ir_0_4        10
     5  '1=2...9=18,A=20...Z=70,+=MORE           ', ! ir_0_71       .5
     6  '1=50...9=450,A=500...Z=1750,+=MORE      ', ! ir_0_1775     .02
     7  '1=100...9=900,A=1000...Z=3500,+=MORE    ', ! ir_0_3550     .01
     8  '1=20...9=180,A=200...Z=700,+=MORE       ', ! ir_0_710      .05
     9  'A=1...Z=26,3=30...9=90,+=100-150,*=MORE ', ! ir_0_26_150    1
     O  '0=0,A=.1...Z=2.6,3=3...9=9,+=10-15      ', ! ir_0_3_15     10
     1  '-=LESS,Z=-78...0=0...9=27,+=MORE        ', ! ir_m80_28     .33
     2  '-=LESS,Z=-260...0=0...9=90,+=MORE       ', ! ir_m265_95    .1
     3  '-=LESS,Z=-520...0=0...9=180,+=MORE      ', ! ir_m530_190   .05
     4  '-=LESS,Z=-1300...0=0...9=450,+=MORE     ', ! ir_m1325_475  .02
     5  '-=LESS,Z=-2600...0=0...9=900,+=MORE     ', ! ir_m2650_950  .01
     6  '-=LESS,Z=-3900...0=0...9=1350,+=MORE    ', ! ir_m3975_1425 .007
     7  '-=LESS,Z=-5200...0=0...9=1800,+=MORE    ', ! ir_m5300_1900 .005
     8  '-=LESS,9=-.9...0=0,A=.1...Z=2.6,+=MORE  ', ! ir_m1_3       10
     9  '-=LESS,9=-45...0=0,A=5...I=45...+=MORE  ', ! ir_m45_130    .2
     O  '-=LESS,9=-90...0=0,A=10...Z=260,+=MORE  ', ! ir_m95_265    .1
     1  '-=LESS,9=-180...A=20...Z=520,+=MORE     ', ! ir_m190_530   .05
     2  '-=LESS,9=-9...0=0,A=1...Z=26,+=MORE     ', ! ir_m9_26       1
     3  '-=LESS,9=-36...0=0,A=4...Z=104,+=MORE   ', ! ir_m38_106    .25
     4  '1=5...9=45,A=50...Z=175,+=MORE          ', ! ir_0_180      .2
     5  '9=-512...1=-2,0=0,A=2,B=4,C=8...+=MORE  '/)! ir_log2       1.
!@var ir_xxxx names for indices to LEGEND indicating the (rounded) range
      integer, parameter, public ::
     &     ir_pct=1, ir_angl=2, ir_0_18=3, ir_0_4=4,
     * ir_0_71=5, ir_0_1775=6, ir_0_3550=7, ir_0_710=8, ir_0_26_150=9,
     * ir_0_3_15=10, ir_m80_28=11, ir_m265_95=12, ir_m530_190=13,
     * ir_m1325_475=14, ir_m2650_950=15, ir_m3975_1425=16,
     * ir_m5300_1900=17, ir_m1_3=18, ir_m45_130=19, ir_m95_265=20,
     * ir_m190_530=21, ir_m9_26=22, ir_m38_106=23, ir_0_180=24,
     * ir_log2=25
!@var fac_legnd = 1/(range_of_1_colorbox)
      real*8, dimension(25), public :: fac_legnd=(/
     1      1d0/5,  1d0/10,    2.d0,   10.d0,   1d0/2,
     6     1d0/50, 1d0/100,  1d0/20,    1.d0,   10.d0,
     1      1d0/3,  1d0/10,  1d0/20,  1d0/50, 1d0/100,
     6    1d0/150, 1d0/200,   10.d0,   1d0/5,  1d0/10,
     1     1d0/20,    1.d0,   1d0/4,   1d0/5,     1d0  /)

!@param CBAR "color bars" for ij-maps
      CHARACTER(LEN=38), PARAMETER, DIMENSION(5), public :: CBAR=(/
     &     ' 0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ+',  ! ib_pos
     &     ' 0123456789ABCDEFGHIJKX               ',  ! ib_pct
     &     '-9876543210ABCDEFGHIJKLMNOPQRSTUVWXYZ+',  ! ib_npp,ib_ntr
     &     ' 0ABCDEFGHIJKLMNOPQRSTUVWXYZ3456789+* ',  ! ib_hyb
     &     '-ZYXWVUTSRQPONMLKJIHGFEDCBA0123456789+'/) ! ib_nnp
!@var ib_xxx indices for color bars
      integer, parameter, public ::
     &     ib_pos=1,ib_pct=2,ib_npp=3,ib_hyb=4,ib_nnp=5
     &     ,ib_ntr=6

!@dbparam isccp_diags: if 1 accumulate ISCCP cloud data (default 0)
      INTEGER, public :: isccp_diags = 0

!@dbparam lh_diags: if 1 accumulate 3D latent heating profiles (default 0)
      INTEGER, public :: lh_diags = 0

!@var SCALE_IJ scaling for weighted AIJ diagnostics
      REAL*8, DIMENSION(KAIJ), public :: SCALE_IJ
!@var NAME_IJ,UNITS_IJ Names/Units of lat/lon IJ diagnostics
      character(len=30), dimension(kaijx), public :: name_ij,units_ij
!@var LNAME_IJ Long names of lat/lon IJ diagnostics
      character(len=80), dimension(kaijx), public :: lname_ij
!@var IW_IJ weighting indices for IJ diagnostics
      integer, dimension(kaij), public :: iw_ij
!@var nwts_ij = number of weight-ij-arrays used in IJ-diagnostics
      integer, parameter, public :: nwts_ij = 8
!@var wt_ij various weight-arrays use in ij-diagnostics
      real*8, dimension(IM,JM,NWTS_IJ), public :: wt_ij
!@var IW_xxx index for weight-array
      integer, parameter, public :: iw_all=1 , iw_ocn=2 , iw_lake=3,
     *   iw_lice=4 , iw_soil=5 , iw_bare=6 , iw_veg=7, iw_land=8
!@var IR_IJ range indices for IJ diagnostics
      integer, dimension(kaij), public :: ir_ij
!@var IA_IJ IDACC indexes for lat/lon IJ diagnostics
      integer, dimension(kaij), public :: ia_ij
!@var [ij]grid_ij 1=primary grid  2=secondary grid
      integer, dimension(kaij), public :: igrid_ij,jgrid_ij

!@var JL_xxx names for JL diagnostic indices
      INTEGER, public ::
     &     jl_mcmflx,jl_srhr,jl_trcr,jl_sshr,jl_trbhr,jl_mchr,jl_totntlh
     *     ,jl_zmfntlh,jl_totvtlh,jl_zmfvtlh,jl_ape,jl_dtdyn,jl_dudfmdrg
     *     ,jl_totcld,jl_dumtndrg,jl_dushrdrg,jl_dumcdrgm10,jl_mcdflx
     *     ,jl_dumcdrgp10,jl_dumcdrgm40,jl_dumcdrgp40,jl_dumcdrgm20
     *     ,jl_dumcdrgp20,jl_sscld,jl_mccld,jl_sdifcoef,jl_dudtsdif
     *     ,jl_gwfirst,jl_dtdtsdrg,jl_epflxv,jl_rhe,jl_epflxn,jl_damdc
     *     ,jl_dammc,jl_40,jl_uepac,jl_vepac,jl_wepac,jl_uwpac,jl_vwpac
     *     ,jl_wwpac,jl_47,jl_zmfntmom,jl_totntmom,jl_mchphas,jl_mcdtotw
     *     ,jl_dudtsdrg,jl_mcldht,jl_trbke,jl_trbdlht,jl_mcheat,jl_mcdry
     *     ,jl_mcdeep,jl_mcshlw
     *     ,jl_cldmc,jl_cldss,jl_csizmc,jl_csizss,jl_dudtvdif
     *     ,jl_wcld,jl_icld,jl_wcod,jl_icod,jl_wcsiz,jl_icsiz
     *     ,jl_cnumwm,jl_cnumim,jl_cnumws,jl_cnumis

!@var JGRID_U, JGRID_KE latitudes at which U-wind and KE diags are defined
!@+   (1 for primary latitudes, 2 for secondary latitudes)
!@+   In future, can be other than the values in GEOM if requested.
      integer, public :: jgrid_u, jgrid_ke
      integer, parameter, public ::
     &     igridc=0,igride=1,jgridc=0,jgride=2,kgridc=0,kgride=4,
     &     ijkgridc=igridc+jgridc+kgridc

!@var SNAME_JL Names of lat-sigma JL diagnostics
      character(len=30), dimension(kajlx), public :: sname_jl
!@var LNAME_JL,UNITS_JL Descriptions/Units of JL diagnostics
      character(len=50), dimension(kajlx), public :: lname_jl,units_jl
!@var SCALE_JL printout scaling factors for JL diagnostics
      REAL*8, dimension(kajlx), public :: scale_jl
!@var IA_JL,JGRID_JL idacc-numbers,gridtypes for JL diagnostics
      integer, dimension(kajlx), public :: ia_jl,jgrid_jl
!@var POW_JL printed output scaled by 10**(-pow_jl)
      integer, dimension(kajlx), public :: pow_jl

!@var NAME_SJL Names of radiative-layer-only SJL diagnostics
      character(len=30), dimension(kasjl), public :: name_sjl
!@var LNAME_SJL,UNITS_SJL Descriptions/Units of SJL diagnostics
      character(len=50), dimension(kasjl), public :: lname_sjl,units_sjl
!@var SCALE_SJL printout scaling factors for SJL diagnostics
      REAL*8, dimension(kasjl), public :: scale_sjl
!@var IA_SJL idacc-numbers for SJL diagnostics
      integer, dimension(kasjl), public :: ia_sjl

!@var JK_xxx names for JK diagnostic indices
      INTEGER, public ::
     &     JK_dpa ,JK_dpb ,JK_temp ,JK_hght
     &    ,JK_q ,JK_theta ,JK_rh ,JK_u
     &    ,JK_v ,JK_zmfke ,JK_totke ,JK_zmfntsh
     &    ,JK_totntsh ,JK_zmfntgeo ,JK_totntgeo ,JK_zmfntlh
     &    ,JK_totntlh ,JK_zmfntke ,JK_totntke ,JK_zmfntmom
     &    ,JK_totntmom ,JK_p2kedpgf ,JK_dpsqr ,JK_nptsavg
     &    ,JK_vvel ,JK_zmfvtdse ,JK_totvtdse ,JK_zmfvtlh
     &    ,JK_totvtlh ,JK_vtgeoeddy ,JK_barekegen ,JK_potvort
     &    ,JK_vtpv ,JK_vtpveddy ,JK_nptsavg1 ,JK_totvtke
     &    ,JK_vtameddy ,JK_totvtam ,JK_sheth ,JK_dudtmadv
     &    ,JK_dtdtmadv ,JK_dudttem ,JK_dtdttem ,JK_epflxncp
     &    ,JK_epflxvcp ,JK_uinst ,JK_totdudt ,JK_tinst
     &    ,JK_totdtdt ,JK_eddvtpt ,JK_cldh2o

!@var SNAME_JK Names of lat-pressure JK diagnostics
      character(len=30), dimension(kajkx), public :: sname_jk
!@var LNAME_JK,UNITS_JK Descriptions/Units of JK diagnostics
      character(len=50), dimension(kajkx), public :: lname_jk,units_jk
!@var SCALE_JK printout scaling factors for JK diagnostics
      REAL*8, dimension(kajkx), public :: scale_jk
!@var IA_JK,JGRID_JK idacc-numbers,gridtypes for JK diagnostics
      integer, dimension(kajkx), public :: ia_jk,jgrid_jk
!@var POW_JK printed output scaled by 10**(-pow_jk)
      integer, dimension(kajkx), public :: pow_jk

!@var IJK_xxx AIJK diagnostic names
      INTEGER, public ::
     &     IJK_U, IJK_V, IJK_DSE, IJK_DP, IJK_T, IJK_Q, IJK_R,
     *     IJK_W, IJK_PF, IJL_CF ,IJK_UV, IJK_VQ, IJK_VT, IJK_UU,
     *     IJK_VV, IJK_TT
     *    ,IJL_LLH,IJL_MCTLH,IJL_MCDLH,IJL_MCSLH
     *    ,IJL_REWM,IJL_REWS,IJL_CDWM,IJL_CDWS,IJL_CWWM,IJL_CWWS
     *    ,IJL_REIM,IJL_REIS,IJL_CDIM,IJL_CDIS,IJL_CWIM,IJL_CWIS

!@var SCALE_IJK scaling for weighted AIJK diagnostics
      REAL*8, DIMENSION(KAIJKx), public :: SCALE_IJK
!@var OFF_IJK offset for weighted AIJK diagnostics
      REAL*8, DIMENSION(KAIJKx), public :: OFF_IJK

!@var NAME_IJK Names of lon-lat-pressure IJK diagnostics
      character(len=30), dimension(kaijkx), public :: name_ijk
!@var LNAME_IJK,UNITS_IJK Descriptions/Units of IJK diagnostics
      character(len=50), dimension(kaijkx), public ::
     &     lname_ijk,units_ijk
!@var jgrid_ijk 1=primary grid  2=secondary grid
      integer, dimension(KAIJKx), public :: jgrid_ijk

      character(len=20), dimension(kwp), public :: name_wave,units_wave
      character(len=80), dimension(kwp), public :: lname_wave

      character(len=20), dimension(kcon), public ::
     &     name_consrv,units_consrv
      character(len=80), dimension(kcon), public :: lname_consrv

      character(len=20), dimension(kail), public :: name_il,units_il
      character(len=80), dimension(kail), public :: lname_il
      real*8, dimension(kail), public :: scale_il
      integer, dimension(kail), public :: ia_il
!@var IL_xxx names for longitude height diagnostics
      INTEGER, public ::
     &     IL_UEQ,IL_VEQ,IL_WEQ,IL_TEQ,IL_QEQ,IL_MCEQ,IL_REQ
     *     ,IL_W50N,IL_T50N,IL_R50N,IL_U50N,IL_W70N,IL_T70N,IL_R70N
     *     ,IL_U70N

      character(len=20), dimension(ndiuvar), public :: name_dd,units_dd
      character(len=80), dimension(ndiuvar), public :: lname_dd
      real*8, dimension(ndiuvar), public :: scale_dd

!@var IDD_xxx names for diurnal diagnostics
      INTEGER, public ::
     &     IDD_ISW, IDD_PALB, IDD_GALB, IDD_ABSA, IDD_ECND,
     *     IDD_SPR, IDD_PT5, IDD_TS, IDD_TG1, IDD_Q5, IDD_QS,
     *     IDD_QG, IDD_SWG, IDD_LWG, IDD_SH, IDD_LH, IDD_HZ0, IDD_UG,
     *     IDD_VG, IDD_WG, IDD_US, IDD_VS, IDD_WS, IDD_CIA, IDD_RIS,
     *     IDD_RIG, IDD_CM, IDD_CH, IDD_CQ, IDD_EDS, IDD_DBL, IDD_DCF,
     *     IDD_LDC, IDD_PR, IDD_EV, IDD_DMC, IDD_SMC, IDD_CL7, IDD_W,
     *     IDD_CCV, IDD_SSP, IDD_MCP,
     &     idd_wtke,idd_wd,idd_wm,idd_wsgcm,idd_wspdf,idd_wtrsh,
     &     idd_u1,idd_v1,idd_uv1,idd_t1,idd_qq1,idd_p1,idd_w1,idd_phi1,
     *     idd_load1,idd_conc1,idd_emis,idd_emis2,idd_tau1,idd_tau_cs1,
     *     idd_sr1,idd_tr1,idd_ws2,idd_ustar,idd_us3,idd_stress,
     *     idd_lmon,idd_rifl,idd_zpbl1,idd_uabl1,idd_vabl1,idd_uvabl1,
     *     idd_tabl1,idd_qabl1,idd_zhat1,idd_e1,idd_km1,idd_ri1,
     *     idd_wet,idd_grav,idd_turb

!@var tf_xxx tsfrez diagnostic names
      INTEGER, public :: tf_day1,tf_last,tf_lkon,tf_lkoff
      character(len=20), dimension(ktsf), public :: name_tsf,units_tsf
      character(len=80), dimension(ktsf), public :: lname_tsf

      character(len=8), dimension(ntype), public :: stype_names=
     &     (/ 'OCEAN   ','OCEANICE','EARTH   ',
     &        'LANDICE ','LAKE    ','LAKEICE ' /)

c idacc-indices of various processes
      integer, parameter, public ::
     &     ia_src=1, ia_rad=2, ia_srf=3, ia_dga=4, ia_d4a=5, ia_d5f=6,
     *     ia_d5d=7, ia_d5s=8, ia_12hr=9, ia_filt=10, ia_ocn=11,
     *     ia_inst=12

!@var PLE,PLM, PLE_DN ref pressures at upper, middle and lower edge
      REAL*8, DIMENSION(LM), public :: PLE
      REAL*8, DIMENSION(LM), public :: PLE_DN
      REAL*8, DIMENSION(LM+LM_REQ), public :: PLM
!@var P1000K scaling to change reference pressure from 1mb to 1000mb
      REAL*8, public :: P1000K
CXXXX inci,incj NOT GRID-INDPENDENT
!@var inci,incj print increments for i and j, so maps/tables fit on page
      integer, parameter, public ::
     &     inci=(im+35)/36,incj=(JM+23)/24, jmby2=jm/2
!@var linect = current line on page of print out
      integer, public :: linect

!@var XWON scale factor for diag. printout needed for Wonderland model
      REAL*8, public :: XWON = TWOPI/(DLON*FIM)

!@var LMOMAX max no. of layers in any ocean
      INTEGER, PARAMETER, public :: LMOMAX=50
!@var ZOC, ZOC1 ocean depths for diagnostics (m) (ONLY FOR DEEP OCEAN)
      REAL*8, public :: ZOC(LMOMAX) = 0. , ZOC1(LMOMAX+1) = 0.


      END MODULE DIAG_COM

      SUBROUTINE ALLOC_DIAG_COM(grid)
!@sum  To allocate arrays whose sizes now need to be determined at
!@+    run time
!@auth NCCS (Goddard) Development Team
!@ver  1.0
      USE DOMAIN_DECOMP, ONLY : DIST_GRID
      USE DOMAIN_DECOMP, ONLY : GET
      USE RESOLUTION, ONLY : IM,LM
      USE MODEL_COM, ONLY : NTYPE,lm_req
      USE DIAG_COM, ONLY : KAJ,KAPJ,KCON,KAJL,KASJL,KAIJ,KAJK,KAIJK,
     &                   KGZ,KOA,KTSF,nwts_ij,KTD,NREG
      USE DIAG_COM, ONLY : SQRTM,AJ_loc,AREGJ_loc,JREG,APJ_loc,AJL_loc
     *     ,ASJL_loc,AIJ_loc,CONSRV_loc,AJK_loc, AIJK_loc, AFLX_ST
     *     ,Z_inst,RH_inst,T_inst,TDIURN,TSFREZ_loc,OA,P_acc,PM_acc


      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(IN) :: grid
      INTEGER :: I_1H, I_0H, J_1H, J_0H
      INTEGER :: IER
      LOGICAL, SAVE :: init = .false.


      If (init) Then
         Return ! Only invoke once
      End If
      init = .true.

      CALL GET( grid, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H  )
      I_0H = grid%I_STRT_HALO
      I_1H = grid%I_STOP_HALO

      ALLOCATE( ! JREG(I_0H:I_1H, J_0H:J_1H),
     &         SQRTM(I_0H:I_1H, J_0H:J_1H),
     &         STAT = IER)

      ALLOCATE( APJ_loc(J_0H:J_1H, KAPJ),
     &         CONSRV_loc(J_0H:J_1H, KCON),
     &         STAT = IER)

      ALLOCATE(AJ_loc(J_0H:J_1H, KAJ, NTYPE),
     &         AREGJ_loc(NREG,J_0H:J_1H,KAJ),
     &         AJL_loc(J_0H:J_1H, LM, KAJL),
     &         ASJL_loc(J_0H:J_1H,LM_REQ,KASJL),
     &         AIJ_loc(I_0H:I_1H,J_0H:J_1H,KAIJ),
     &         AJK_loc(J_0H:J_1H,LM,KAJK),
     &         Z_inst(KGZ,I_0H:I_1H,J_0H:J_1H),
     &         RH_inst(KGZ,I_0H:I_1H,J_0H:J_1H),
     &         T_inst(KGZ,I_0H:I_1H,J_0H:J_1H),
     &         TSFREZ_loc(I_0H:I_1H,J_0H:J_1H,KTSF),
     &         P_acc(I_0H:I_1H,J_0H:J_1H),
     &         PM_acc(I_0H:I_1H,J_0H:J_1H),
     &         TDIURN(I_0H:I_1H,J_0H:J_1H,KTD),
     &         OA(I_0H:I_1H,J_0H:J_1H,KOA),
     &         STAT = IER)

      ALLOCATE( AIJK_loc(I_0H:I_1H,J_0H:J_1H,LM,KAIJK),
     &         AFLX_ST(LM+LM_REQ+1,I_0H:I_1H,J_0H:J_1H,5),
     &         STAT = IER)

      RETURN
      END SUBROUTINE ALLOC_DIAG_COM

      SUBROUTINE ALLOC_DIAG_COM_glob
!@sum  To allocate arrays global arrays only when needed
!@auth NCCS (Goddard) Development Team
!@ver  1.0
      USE RESOLUTION, ONLY : IM,JM,LM
      USE MODEL_COM, ONLY : NTYPE,lm_req
      USE DOMAIN_DECOMP, Only : AM_I_ROOT
      USE DIAG_COM, ONLY : KAJ,KAPJ,KCON,KAJL,KASJL,KAIJ,KAJK,KAIJK,
     &                   KGZ,KOA,KTSF,nwts_ij,KTD,NREG
      USE DIAG_COM, ONLY : AJ,AREGJ,JREG,APJ,AJL
     *     ,ASJL,AIJ,CONSRV,AJK, AIJK
     *     ,TSFREZ

      IMPLICIT NONE
      INTEGER :: IER

      if(.not.AM_I_ROOT()) return

c     ALLOCATE( ! JREG(IM, JM))
      ALLOCATE( APJ(JM, KAPJ),
     &         CONSRV(JM, KCON),
     &         STAT = IER)

      ALLOCATE(AJ(JM, KAJ, NTYPE),
     &         AREGJ(NREG,JM,KAJ),
     &         AJL(JM, LM, KAJL),
     &         ASJL(JM,LM_REQ,KASJL),
     &         AIJ(IM,JM,KAIJ),
     &         AJK(JM,LM,KAJK),
     &         TSFREZ(IM,JM,KTSF),
     &         STAT = IER)


      ALLOCATE( AIJK(IM,JM,LM,KAIJK),
     &         STAT = IER)

      RETURN
      END SUBROUTINE ALLOC_DIAG_COM_glob

      SUBROUTINE DEALLOC_DIAG_COM_glob
!@sum  To deallocate global arrays not currently needed
!@auth NCCS (Goddard) Development Team
!@ver  1.0
      USE DOMAIN_DECOMP, Only : AM_I_ROOT
      USE DIAG_COM, ONLY : AJ,AREGJ,JREG,APJ,AJL
     *     ,ASJL,AIJ,CONSRV,AJK, AIJK
     *     ,TSFREZ

      IMPLICIT NONE

      if(.not.AM_I_ROOT()) return

c     DEALLOCATE( ! JREG)
      DEALLOCATE( AJ,AREGJ,APJ,AJL,ASJL,AIJ,CONSRV,AJK, AIJK,TSFREZ )

      RETURN
      END SUBROUTINE DEALLOC_DIAG_COM_glob

      SUBROUTINE io_diags(kunit,it,iaction,ioerr)
!@sum  io_diag reads and writes diagnostics to file
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : ioread,ioread_single,irerun
     *    ,iowrite,iowrite_mon,iowrite_single,lhead, idacc,nsampl
     *    ,Kradia
      USE DIAG_COM
      USE DOMAIN_DECOMP, Only : grid, GET, PACK_DATA, UNPACK_DATA
      USE DOMAIN_DECOMP, Only : PACK_DATAj, UNPACK_DATAj, LOG_PARALLEL
      USE DOMAIN_DECOMP, Only : PACK_COLUMN, UNPACK_COLUMN
      USE DOMAIN_DECOMP, Only : AM_I_ROOT, CHECKSUM, CHECKSUMj
      USE DOMAIN_DECOMP, Only : ESMF_BCAST
      IMPLICIT NONE

!@param KACC total number of diagnostic elements
      INTEGER, PARAMETER :: KACC= JM*KAJ*NTYPE + JM*NREG*KAJ + JM*KAPJ
     *     + JM*LM*KAJL + JM*LM_REQ*KASJL + IM*JM*KAIJ +
     *     IM*LM*KAIL + NEHIST*HIST_DAYS + JM*KCON +
     *     (IMH+1)*KSPECA*NSPHER + KTPE*NHEMI + HR_IN_DAY*NDIUVAR*NDIUPT
     *     + RE_AND_IM*Max12HR_sequ*NWAV_DAG*KWP + JM*LM*KAJK +
     *     IM*JM*LM*KAIJK+ntau*npres*nisccp
#ifndef NO_HDIURN
     *     + HR_IN_MONTH*NDIUVAR*NDIUPT
#endif
!@var AJ4,...,AFLX4 real*4 dummy arrays needed for postprocessing only
      ! REAL*4 AJ4(JM,KAJ,NTYPE),AREGJ4(NREG,JM,KAJ),APJ4(JM,KAPJ)
      ! REAL*4 AJL4(JM,LM,KAJL),ASJL4(JM,LM_REQ,KASJL),AIJ4(IM,JM,KAIJ)
      ! REAL*4 AIL4(IM,LM,KAIL),ENERGY4(NEHIST,HIST_DAYS)
      ! REAL*4 CONSRV4(JM,KCON),SPECA4(IMH+1,KSPECA,NSPHER)
      ! REAL*4 ATPE4(KTPE,NHEMI),ADIURN4(HR_IN_DAY,NDIUVAR,NDIUPT)
      ! REAL*4 WAVE4(RE_AND_IM,Max12HR_sequ,NWAV_DAG,KWP)
      ! REAL*4 AJK4(JM,LM,KAJK),AIJK4(IM,JM,LM,KAIJK)
      ! REAL*4 AISCCP4(ntau,npres,nisccp)
      ! REAL*4 TSFREZ4(IM,JM,KTSF),AFLX4(LM+LM_REQ+1,IM,JM,5)
      REAL*4,allocatable, dimension(:,:,:) :: AJ4,AREGJ4, AJL4,ASJL4
     *            ,AIJ4,AIL4,AJK4, SPECA4, ADIURN4, AISCCP4, TSFREZ4
      REAL*4,allocatable, dimension(:,:) :: APJ4,ENERGY4,CONSRV4,ATPE4
      REAL*4,allocatable, dimension(:,:,:,:) :: WAVE4, AIJK4, AFLX4
#ifndef NO_HDIURN
      ! REAL*4 HDIURN4(HR_IN_MONTH,NDIUVAR,NDIUPT)
      REAL*4,allocatable ::  HDIURN4(:,:,:)
#endif
      integer monac1(12),i_ida,i_xtra,it_check
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
      ! REAL*8 :: AFLX_ST_GLOB(LM+LM_REQ+1,IM,JM,5)
      REAL*8,allocatable :: AFLX_ST_GLOB(:,:,:,:)

      ! REAL*8 :: TDIURN_glob(IM, JM, KTD)
      ! REAL*8 :: OA_glob    (IM, JM, KOA)
      REAL*8,allocatable :: TDIURN_glob(:,:,:), OA_glob(:,:,:)
      INTEGER :: J_0, J_1

      CALL GET( grid, J_STRT=J_0, J_STOP=J_1  )

      if(kradia.gt.0) then
        write (MODULE_HEADER(LHEAD+1:80),'(a6,i8,a20,i3,a7)')
     *   '#acc(=',idacc(2),') R8:SU.SD.TU.TD.dT(',lm+lm_req+1,',ijM,5)'

        IF (AM_I_ROOT()) allocate (AFLX_ST_glob(LM+LM_REQ+1,IM,JM,5))

        SELECT CASE (IACTION)
        CASE (IOWRITE)            ! output to standard restart file
          CALL PACK_COLUMN(grid, AFLX_ST, AFLX_ST_glob)
          IF (AM_I_ROOT()) THEN
            WRITE (kunit,err=10) MODULE_HEADER,idacc(2),AFLX_ST_glob,it
          END IF
        CASE (IOWRITE_SINGLE)     ! output in single precision
          MODULE_HEADER(LHEAD+18:LHEAD+18) = '4'
          MODULE_HEADER(LHEAD+44:80) = ',monacc(12)'
          CALL PACK_COLUMN(grid, AFLX_ST, AFLX_ST_glob)
          IF (AM_I_ROOT()) THEN
            WRITE (kunit,err=10) MODULE_HEADER,idacc(2),
     *           REAL(AFLX_ST_glob,KIND=4), monacc,it
          END IF
        CASE (IOWRITE_MON)        ! output to end-of-month restart file
          MODULE_HEADER(LHEAD+1:80) = 'itime '
          IF (AM_I_ROOT()) THEN
            WRITE (kunit,err=10) MODULE_HEADER,it
          END IF
        CASE (ioread)           ! input from restart file
          if (AM_I_ROOT()) THEN
           READ (kunit,err=10) HEADER,idacc(2),AFLX_ST_glob,it
           IF (HEADER(1:LHEAD).NE.MODULE_HEADER(1:LHEAD)) THEN
            PRINT*,"Discrepancy in module version ",HEADER,MODULE_HEADER
            GO TO 10
           END IF
          endif
          CALL UNPACK_COLUMN(grid, AFLX_ST_glob, AFLX_ST)
          CALL ESMF_BCAST(grid, idacc)
          CALL ESMF_BCAST(grid, it   )

        CASE (IOREAD_SINGLE)      !
!ESMF-- Allow all processes to read to avoid scattering monac1 and idac1.
          CALL PACK_COLUMN(grid, AFLX_ST, AFLX_ST_glob)
          if (AM_I_ROOT()) then
             allocate (AFLX4(LM+LM_REQ+1,IM,JM,5))
             READ (kunit,err=10) HEADER,idac1(2),AFLX4,monac1
             AFLX_ST_glob = AFLX_ST_glob + AFLX4
             deallocate (AFLX4)
          end if
          CALL UNPACK_COLUMN(grid, AFLX_ST_glob, AFLX_ST)
          CALL ESMF_BCAST(grid,idac1)
          CALL ESMF_BCAST(grid,monac1)
          IDACC(2) = IDACC(2) + IDAC1(2)
          monacc = monacc + monac1
        END SELECT

        if (AM_I_ROOT()) deallocate (AFLX_ST_glob)

        return
      end if

C**** The regular model (Kradia le 0)
      write (MODULE_HEADER(LHEAD+1:LHEAD+15),'(a10,i4,a1)')
     *   'I/R8 keys(',1+NKEYNR*NKEYMO,')'             ! keyct,keynr(:,:)
      i_ida = Lhead + 10+4+1 + 10+2+1 + 1
      write (MODULE_HEADER(LHEAD+10+4+1+1:i_ida-1),'(a10,i2,a1)')
     *   ',TSFR(IJM,',KTSF,')'
      write (MODULE_HEADER(i_ida:i_ida+9),'(a7,i2,a1)')
     *   ',idacc(',nsampl,')'
      write (MODULE_HEADER(i_ida+9+1:i_ida+9 + 5+8+1),'(a5,i8,a1)')
     *   ',acc(',kacc,')'
      i_xtra = i_ida+9 + 5+8+1 + 1

      call alloc_diag_com_glob
      if (AM_I_ROOT()) allocate ( TDIURN_glob(IM, JM, KTD)
     *                           ,OA_glob(IM, JM, KOA))

      SELECT CASE (IACTION)
      CASE (IOWRITE)            ! output to standard restart file
        write (MODULE_HEADER(i_xtra:80),             '(a7,i2,a)')
     *   ',x(IJM,',KTD+KOA,')'  ! make sure that i_xtra+7+2 < 80

        Call Gather_Diagnostics()

        If (AM_I_ROOT()) THEN
          WRITE (kunit,err=10) MODULE_HEADER,keyct,KEYNR,TSFREZ,
     *     idacc, AJ,AREGJ,APJ,AJL,ASJL,AIJ,
     *     AIL, ENERGY,CONSRV,
     *     SPECA,ATPE,ADIURN,WAVE,AJK,AIJK,AISCCP,
#ifndef NO_HDIURN
     *     HDIURN,
#endif
     *     TDIURN_glob,OA_glob,it
        END IF
      CASE (IOWRITE_SINGLE)     ! output in single precision
        MODULE_HEADER(LHEAD+1:LHEAD+4) = 'I/R4'
        MODULE_HEADER(i_xtra:80) = ',monacc(12)'

        Call Gather_Diagnostics()

        If (AM_I_ROOT()) THEN
          WRITE (kunit,err=10) MODULE_HEADER,
     *     keyct,KEYNR,REAL(TSFREZ,KIND=4),   idacc,
     *     REAL(AJ,KIND=4),REAL(AREGJ,KIND=4),REAL(APJ,KIND=4),
     *     REAL(AJL,KIND=4),REAL(ASJL,KIND=4),
     *     REAL(AIJ,KIND=4),REAL(AIL,KIND=4),
     *     REAL(ENERGY,KIND=4), REAL(CONSRV,KIND=4),
     *     REAL(SPECA,KIND=4),REAL(ATPE,KIND=4),REAL(ADIURN,KIND=4),
     *     REAL(WAVE,KIND=4),REAL(AJK,KIND=4),
     *     REAL(AIJK,KIND=4),REAL(AISCCP,KIND=4),
#ifndef NO_HDIURN
     *     REAL(HDIURN,KIND=4),
#endif
     *     monacc,it
        END IF
      CASE (IOWRITE_MON)        ! output to end-of-month restart file
        MODULE_HEADER(i_ida:80) = ',it '

        CALL PACK_DATA(grid, TSFREZ_loc, TSFREZ)
        If (AM_I_ROOT()) THEN
          WRITE (kunit,err=10) MODULE_HEADER,keyct,KEYNR,TSFREZ,it
        END IF
      CASE (ioread)           ! input from restart file
        if (AM_I_ROOT()) Then
          READ (kunit,err=10) HEADER,keyct,KEYNR,TSFREZ,
     *       idacc, AJ,AREGJ,APJ,AJL,ASJL,AIJ,AIL,
     *       ENERGY,CONSRV,SPECA,ATPE,ADIURN,WAVE,AJK,AIJK,AISCCP,
#ifndef NO_HDIURN
     *       HDIURN,
#endif
     *       TDIURN_glob,OA_glob,it
        IF (HEADER(1:LHEAD).NE.MODULE_HEADER(1:LHEAD)) THEN
          PRINT*,"Discrepancy in module version ",HEADER,MODULE_HEADER
          GO TO 10
        END IF
        END IF

        Call BCAST_Scalars()
        Call Scatter_Diagnostics()

      CASE (IOREAD_SINGLE)      !
        call Gather_Diagnostics()  ! to keep global arrays in up-to-date
        If (AM_I_ROOT()) Then
          call alloc_diag_r4
          READ (kunit,err=10) HEADER,keyct,KEYNR,TSFREZ4,
     *         idac1, AJ4,AREGJ4,APJ4,AJL4,ASJL4,AIJ4,AIL4,ENERGY4
     *         ,CONSRV4,SPECA4,ATPE4,ADIURN4,WAVE4,AJK4,AIJK4,AISCCP4,
#ifndef NO_HDIURN
     *       HDIURN4,
#endif
     *       monac1,it_check
          if(it.ne.it_check) then
            PRINT*,"io_diags: compare aj,aj4, ... dimensions"
            GO TO 10            ! or should that be just a warning ??
          end if

        ! copy or add in to full precision global variables
        ! First "non-distributed" arrays
          AIL=AIL+AIL4
          ENERGY=ENERGY+ENERGY4
          SPECA=SPECA+SPECA4 ; ATPE=ATPE+ATPE4 ; ADIURN=ADIURN+ADIURN4
          WAVE=WAVE+WAVE4
          AISCCP=AISCCP+AISCCP4
#ifndef NO_HDIURN
          HDIURN=HDIURN+HDIURN4
#endif
          ! Now for the global versions of the distributed arrays
          TSFREZ= TSFREZ4       ! not accumulated
          AJ    = AJ     + AJ4
          AREGJ = AREGJ  + AREGJ4
          APJ   = APJ    + APJ4
          AJL   = AJL    + AJL4
          ASJL  = ASJL   + ASJL4
          AIJ   = AIJ    + AIJ4
          CONSRV= CONSRV + CONSRV4
          AJK   = AJK    + AJK4
          AIJK  = AIJK   + AIJK4

          IDACC = IDACC + IDAC1
          call dealloc_diag_r4
        End If

        ! Send to other processors - update distributed arrays
        Call BCAST_Scalars()

        Call Scatter_Diagnostics()

!@var idacc(5) is the length of a time series (daily energy history).
!****   If combining acc-files, rather than concatenating these series,
!****   we average their beginnings (up to the length of the shortest)
        Kcomb = Kcomb + 1          ! reverse addition, take min instead
        if (Kcomb.gt.1) IDACC(5) = MIN(IDACC(5)-IDAC1(5),IDAC1(5))
        monacc = monacc + monac1
      CASE (irerun)      ! only keynr,tsfrez needed at beg of acc-period
        If (AM_I_ROOT()) Then
          READ (kunit,err=10) HEADER,keyct,KEYNR,TSFREZ  ! 'it' not read
          IF (HEADER(1:LHEAD).NE.MODULE_HEADER(1:LHEAD)) THEN
            PRINT*,"Discrepancy in module version ",HEADER,MODULE_HEADER
            GO TO 10
          END IF
        End If
        CALL ESMF_BCAST(grid, keyct )
        CALL ESMF_BCAST(grid, KEYNR )
        CALL UNPACK_DATA(grid,  TSFREZ, TSFREZ_loc)
      END SELECT

      call dealloc_diag_com_glob
      if ( am_i_root() ) deallocate(TDIURN_glob, OA_glob)
      RETURN

 10   IOERR=1
      call dealloc_diag_com_glob
      if ( am_i_root() ) deallocate(TDIURN_glob, OA_glob)
      RETURN

      Contains

      Subroutine BCAST_Scalars()
        CALL ESMF_BCAST(grid, keyct )
        CALL ESMF_BCAST(grid, KEYNR )
        CALL ESMF_BCAST(grid, idacc )
        CALL ESMF_BCAST(grid, AIL   )
        CALL ESMF_BCAST(grid, ENERGY)
        CALL ESMF_BCAST(grid, SPECA )
        CALL ESMF_BCAST(grid, ATPE  )
        CALL ESMF_BCAST(grid, ADIURN)
        CALL ESMF_BCAST(grid, WAVE  )
        CALL ESMF_BCAST(grid, AISCCP)
#ifndef NO_HDIURN
        CALL ESMF_BCAST(grid, HDIURN)
#endif
        CALL ESMF_BCAST(grid, it    )
      End Subroutine BCAST_Scalars

      Subroutine Scatter_Diagnostics()
        CALL UNPACK_DATA(grid,  TSFREZ, TSFREZ_loc)
        CALL UNPACK_DATAj(grid, AJ,     AJ_loc)
        CALL UNPACK_DATA(grid, AREGJ,  AREGJ_loc)
        CALL UNPACK_DATAj(grid, APJ,    APJ_loc)
        CALL UNPACK_DATAj(grid, AJL,    AJL_loc)
        CALL UNPACK_DATAj(grid, ASJL,   ASJL_loc)
        CALL UNPACK_DATA(grid,  AIJ,    AIJ_loc)
        CALL UNPACK_DATAj(grid, CONSRV, CONSRV_loc)
        CALL UNPACK_DATAj(grid, AJK,    AJK_loc)
        CALL UNPACK_DATA(grid,  AIJK,   AIJK_loc)
        CALL UNPACK_DATA(grid,  TDIURN_glob, TDIURN)
        CALL UNPACK_DATA(grid,  OA_glob,     OA)
      End Subroutine Scatter_Diagnostics

      Subroutine Gather_Diagnostics()
        CALL PACK_DATA(grid,  TSFREZ_loc, TSFREZ)
        CALL PACK_DATAj(grid, AJ_loc,     AJ)
        CALL PACK_DATA(grid, AREGJ_loc,  AREGJ)
        CALL PACK_DATAj(grid, APJ_loc,    APJ)
        CALL PACK_DATAj(grid, AJL_loc,    AJL)
        CALL PACK_DATAj(grid, ASJL_loc,   ASJL)
        CALL PACK_DATA(grid,  AIJ_loc,    AIJ)
        CALL PACK_DATAj(grid, CONSRV_loc, CONSRV)
        CALL PACK_DATAj(grid, AJK_loc,    AJK)
        CALL PACK_DATA(grid,  AIJK_loc,   AIJK)
        CALL PACK_DATA(grid,  TDIURN, TDIURN_glob)
        CALL PACK_DATA(grid,  OA, OA_glob)
      End Subroutine Gather_Diagnostics

      Subroutine alloc_diag_r4
        allocate (AJ4(JM,KAJ,NTYPE),AREGJ4(NREG,JM,KAJ),APJ4(JM,KAPJ))
        allocate (AJL4(JM,LM,KAJL),ASJL4(JM,LM_REQ,KASJL))
        allocate (AIJ4(IM,JM,KAIJ))
        allocate (AIL4(IM,LM,KAIL),ENERGY4(NEHIST,HIST_DAYS))
        allocate (CONSRV4(JM,KCON),SPECA4(IMH+1,KSPECA,NSPHER))
        allocate (ATPE4(KTPE,NHEMI),ADIURN4(HR_IN_DAY,NDIUVAR,NDIUPT))
        allocate (WAVE4(RE_AND_IM,Max12HR_sequ,NWAV_DAG,KWP))
        allocate (AJK4(JM,LM,KAJK),AIJK4(IM,JM,LM,KAIJK))
        allocate (AISCCP4(ntau,npres,nisccp))
        allocate (TSFREZ4(IM,JM,KTSF))
#ifndef NO_HDIURN
        allocate (HDIURN4(HR_IN_MONTH,NDIUVAR,NDIUPT))
#endif
      End Subroutine alloc_diag_r4

      Subroutine dealloc_diag_r4
        deallocate (AJ4,AREGJ4,APJ4,  AJL4,ASJL4,AIJ4,AIL4)
        deallocate (ENERGY4,CONSRV4,SPECA4,ATPE4,ADIURN4,WAVE4)
        deallocate (AJK4,AIJK4, AISCCP4,TSFREZ4)
#ifndef NO_HDIURN
        deallocate (HDIURN4)
#endif
      End Subroutine dealloc_diag_r4


      END SUBROUTINE io_diags

      SUBROUTINE aPERIOD (JMON1,JYR1,months,years,moff,  aDATE,LDATE)
!@sum  aPERIOD finds a 7 or 12-character name for an accumulation period
!@+   if the earliest month is NOT the beginning of the 2-6 month period
!@+   the name will reflect that fact ONLY for 2 or 3-month periods
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
!@var moff = # of months from beginning of period to JMON1 if months<12
      integer moff
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
      if (moff.gt.0.and.months.le.3) then  ! earliest month is NOT month
        JMONL = 1 + mod(10+jmon1,12)       ! 1 of the 2-3 month period
        yr1=JYR1
        if (jmon1.gt.1) yr1=yr1+1
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
      IF (months.eq.2) then
        aDATE(2:2)='+'
        return
      end if
      if (months.eq.3) then
        JMONM = JMONL-1
        if (moff.eq.1) jmonm = jmon1+1
        if (jmonm.gt.12) jmonm = jmonm-12
        if (jmonm.le.0 ) jmonm = jmonm+12
        aDATE(2:2)=AMONTH(JMONM)(1:1)
        return
      end if
      if (moff.gt.0) then  ! can't tell non-consec. from consec. periods
        jmon1 = jmon1-moff
        if (jmon1.le.0) jmon1 = jmon1+12
        JMONL=JMON1+months-1
        if (jmonl.gt.12) jmonl = jmonl-12
        aDATE(1:1)=AMONTH(JMON1)(1:1)
        aDATE(3:3)=AMONTH(JMONL)(1:1)
      end if
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


C**** Should this be here?
      SUBROUTINE SET_J_BUDG
!@sum set_j_budg definition for grid points map to budget-grid zonal means
!@auth Gavin Schmidt
      USE GEOM, only : j_budg, lat_dg
      USE DIAG_COM, only : jm_budg
      USE DOMAIN_DECOMP, only :GRID,GET
      IMPLICIT NONE
!@var I,J are atm grid point values for the accumulation
      INTEGER :: I,J,J_0,J_1,I_0,I_1,J_0H,J_1H,I_0H,I_1H
      INTEGER :: IER
 
C**** define atmospheric grid
      CALL GET(grid, J_STRT=J_0,J_STOP=J_1, J_STRT_HALO=J_0H,
     *     J_STOP_HALO=J_1H  ) 
      I_0H = grid%I_STRT_HALO ; I_1H = grid%I_STOP_HALO
      I_0 = grid%I_STRT ; I_1 = grid%I_STOP

      ALLOCATE( J_BUDG(I_0H:I_1H, J_0H:J_1H), STAT = IER) 

      DO J=J_0,J_1
        J_BUDG(I_0:I_1,J)=J   ! temporary
C**** this should be valid for all lat/lon grids
c        J_BUDG(I_0:I_1,J)=NINT(1+(lat_dg(J,1)+90)*(JM_BUDG-1)/180.)
      END DO

      RETURN
      END SUBROUTINE SET_J_BUDG

      SUBROUTINE INC_AJ(I,J,ITYPE,J_DIAG,ACC)
!@sum inc_aj grid dependent incrementer for zonal mean budget diags
!@auth Gavin Schmidt
      USE DIAG_COM, only : aj=>aj_loc
      USE GEOM, only : j_budg
      IMPLICIT NONE
!@var I,J are atm grid point values for the accumulation
!@var ITYPE is the surface type
      INTEGER, INTENT(IN) :: I,J,ITYPE
!@var J_DIAG is placing of diagnostic element
      INTEGER, INTENT(IN) :: J_DIAG
!@var ACC value of the diagnostic to be accumulated
      REAL*8, INTENT(IN) :: ACC

C**** accumulate I,J value on the budget grid using j_budg to assign
C**** each point to a zonal mean (not bitwise reprodcible for MPI).
      AJ(J_BUDG(I,J),J_DIAG,ITYPE) = AJ(J_BUDG(I,J),J_DIAG,ITYPE)+ACC
      
      RETURN
      END SUBROUTINE INC_AJ

      SUBROUTINE INC_AREG(I,J,JR,J_DIAG,ACC)
!@sum inc_areg grid dependent incrementer for regional budget diags
!@auth Gavin Schmidt
      USE DIAG_COM, only : aregj=>aregj_loc
      USE GEOM, only : j_budg
      IMPLICIT NONE
!@var I,J are atm grid point values for the accumulation
!@var JR is the region
      INTEGER, INTENT(IN) :: I,J,JR
!@var J_DIAG is placing of diagnostic element
      INTEGER, INTENT(IN) :: J_DIAG
!@var ACC value of the diagnostic to be accumulated
      REAL*8, INTENT(IN) :: ACC

C**** accumulate I,J value on the budget grid using j_budg to assign
C**** each point to a zonal mean (not bitwise reproducible for MPI).
      AREGJ(JR,J_BUDG(I,J),J_DIAG) = AREGJ(JR,J_BUDG(I,J),J_DIAG)+ACC
      
      RETURN
      END SUBROUTINE INC_AREG
