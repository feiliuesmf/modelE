#include "rundeck_opts.h"

!@sum  TRACER_COM: Exists alone to minimize the number of dependencies
!@+    This version for simple trace gases/chemistry/isotopes
!@auth Jean Lerner/Gavin Schmidt
!@ver  1.0

      MODULE TRACER_COM
!@sum  TRACER_COM tracer variables
!@auth Jean Lerner
!@ver  1.0
C                 
      USE QUSDEF, only: nmom
      USE MODEL_COM, only: im,jm,lm
C      
      IMPLICIT NONE
      SAVE

C**** Each tracer has a variable name and a unique index
!@param NTM number of tracers
!@var TRNAME: Name for each tracer >>> MUST BE LEFT-JUSTIFIED <<<
#ifdef TRACERS_SPECIAL_O18
      integer, parameter :: ntm=3
      character*8, parameter :: trname(ntm)=(/
     *     'Water   ','H2O18   ','HDO     '/)  !,'HTO     '/)
#else     
#if (defined TRACERS_SPECIAL_Lerner) && (defined TRACERS_WATER) 
      integer, parameter :: ntm=10
      character*8, parameter :: trname(ntm)= (/
     *     'Air     ','SF6     ','Rn222   ','CO2     ','N2O     ',
     *     'CFC11   ','14CO2   ','CH4     ','O3      ','Water   '/)
#else 
#if defined TRACERS_SPECIAL_Lerner
      integer, parameter :: ntm=9
      character*8, parameter :: trname(ntm)= (/
     *     'Air     ','SF6     ','Rn222   ','CO2     ','N2O     ',
     *     'CFC11   ','14CO2   ','CH4     ','O3      '/)
#else
#if (defined TRACERS_DUST) && (defined TRACERS_SPECIAL_Shindell) && (defined TRACERS_AEROSOLS_Koch)
!@var ntm_chem number of drew-only tracers
#ifdef regional_Ox_tracers
      integer, parameter :: ntm=32,ntm_chem=21,ntm_dust=4
C Note: please always put the regional Ox tracers at the end,
C starting with OxREG1 to facilitate loops. Also, Ox must be tracer.
      character*8, parameter :: trname(ntm)=(/
     *    'Ox      ','NOx     ','N2O5    ','HNO3    ','H2O2    ',
     *    'CH3OOH  ','HCHO    ','HO2NO2  ','CO      ','CH4     ',
     *    'PAN     ','Isoprene','AlkylNit','Alkenes ','Paraffin',
     *    'OxREG1  ','OxREG2  ','OxREG3  ','OxREG4  ','OxREG5  ',
     *    'OxREG6  ',
     *    'DMS     ','MSA     ','SO2     ','SO4     ','H2O2_s  ',
     *    'seasalt1','seasalt2',
     *    'Clay    ','Silt1   ','Silt2   ','Silt3   '/)
#else
      integer, parameter :: ntm=26,ntm_chem=15,ntm_dust=4
      character*8, parameter :: trname(ntm)=(/
     *    'Ox      ','NOx     ','N2O5    ','HNO3    ','H2O2    ',
     *    'CH3OOH  ','HCHO    ','HO2NO2  ','CO      ','CH4     ',
     *    'PAN     ','Isoprene','AlkylNit','Alkenes ','Paraffin',
     *    'DMS     ','MSA     ','SO2     ','SO4     ','H2O2_s  ',
     *    'seasalt1','seasalt2',
     *    'Clay    ','Silt1   ','Silt2   ','Silt3   '/)
#endif
#else
#if (defined TRACERS_DUST) && (defined TRACERS_SPECIAL_Shindell)
!@var ntm_chem number of drew-only tracers
#ifdef regional_Ox_tracers
      integer, parameter :: ntm=25,ntm_chem=21,ntm_dust=4
C Note: please always put the regional Ox tracers at the end,
C starting with OxREG1 to facilitate loops. Also, Ox must be tracer.
      character*8, parameter :: trname(ntm)=(/
     *    'Ox      ','NOx     ','N2O5    ','HNO3    ','H2O2    ',
     *    'CH3OOH  ','HCHO    ','HO2NO2  ','CO      ','CH4     ',
     *    'PAN     ','Isoprene','AlkylNit','Alkenes ','Paraffin',
     *    'OxREG1  ','OxREG2  ','OxREG3  ','OxREG4  ','OxREG5  ',
     *    'OxREG6  ',
     *    'Clay    ','Silt1   ','Silt2   ','Silt3   '/)
#else
      integer, parameter :: ntm=19,ntm_chem=15,ntm_dust=4
      character*8, parameter :: trname(ntm)=(/
     *    'Ox      ','NOx     ','N2O5    ','HNO3    ','H2O2    ',
     *    'CH3OOH  ','HCHO    ','HO2NO2  ','CO      ','CH4     ',
     *    'PAN     ','Isoprene','AlkylNit','Alkenes ','Paraffin',
     *    'Clay    ','Silt1   ','Silt2   ','Silt3   '/)
#endif
#else
#if (defined TRACERS_AEROSOLS_Koch) && (defined TRACERS_SPECIAL_Shindell)
!@var ntm_chem number of drew-only tracers
#ifdef regional_Ox_tracers
      integer, parameter :: ntm=26,ntm_chem=21
C Note: please always put the regional Ox tracers at the end,
C starting with OxREG1 to facilitate loops. Also, Ox must be tracer.
      character*8, parameter :: trname(ntm)=(/
     *    'Ox      ','NOx     ','N2O5    ','HNO3    ','H2O2    ',
     *    'CH3OOH  ','HCHO    ','HO2NO2  ','CO      ','CH4     ',
     *    'PAN     ','Isoprene','AlkylNit','Alkenes ','Paraffin',
     *    'OxREG1  ','OxREG2  ','OxREG3  ','OxREG4  ','OxREG5  ',
     *    'OxREG6  ',
     *    'DMS     ','MSA     ','SO2     ','SO4     ','H2O2_s  '/)
#else
      integer, parameter :: ntm=20,ntm_chem=15
      character*8, parameter :: trname(ntm)=(/
     *    'Ox      ','NOx     ','N2O5    ','HNO3    ','H2O2    ',
     *    'CH3OOH  ','HCHO    ','HO2NO2  ','CO      ','CH4     ',
     *    'PAN     ','Isoprene','AlkylNit','Alkenes ','Paraffin',
     *    'DMS     ','MSA     ','SO2     ','SO4     ','H2O2_s  '/)
#endif
#else
#ifdef TRACERS_SPECIAL_Shindell
#ifdef SHINDELL_STRAT_CHEM
#ifdef regional_Ox_tracers
      integer, parameter :: ntm=31,ntm_chem=31
C Note: please always put the regional Ox tracers at the end,
C starting with OxREG1 to facilitate loops. Also, Ox must be tracer.
      character*8, parameter :: trname(ntm)=(/
     *    'Ox      ','NOx     ','ClOx    ','BrOx    ','N2O5    ',
     *    'HNO3    ','H2O2    ','CH3OOH  ','HCHO    ','HO2NO2  ',
     *    'CO      ','CH4     ','PAN     ','Isoprene','AlkylNit',
     *    'Alkenes ','Paraffin','HCl     ','HOCl    ','ClONO2  ',
     *    'HBr     ','HOBr    ','BrONO2  ','N2O     ','CFC     ',
     *    'OxREG1  ','OxREG2  ','OxREG3  ','OxREG4  ','OxREG5  ',
     *    'OxREG6  '/)
#else
      integer, parameter :: ntm=25,ntm_chem=25
      character*8, parameter :: trname(ntm)=(/
     *    'Ox      ','NOx     ','ClOx    ','BrOx    ','N2O5    ',
     *    'HNO3    ','H2O2    ','CH3OOH  ','HCHO    ','HO2NO2  ',
     *    'CO      ','CH4     ','PAN     ','Isoprene','AlkylNit',
     *    'Alkenes ','Paraffin','HCl     ','HOCl    ','ClONO2  ',
     *    'HBr     ','HOBr    ','BrONO2  ','N2O     ','CFC     '/)
#endif
#else
#ifdef regional_Ox_tracers
      integer, parameter :: ntm=21,ntm_chem=21
C Note: please always put the regional Ox tracers at the end,
C starting with OxREG1 to facilitate loops. Also, Ox must be tracer.
      character*8, parameter :: trname(ntm)=(/
     *    'Ox      ','NOx     ','N2O5    ','HNO3    ','H2O2    ',
     *    'CH3OOH  ','HCHO    ','HO2NO2  ','CO      ','CH4     ',
     *    'PAN     ','Isoprene','AlkylNit','Alkenes ','Paraffin',
     *    'OxREG1  ','OxREG2  ','OxREG3  ','OxREG4  ','OxREG5  ',
     *    'OxREG6  '/)
#else
      integer, parameter :: ntm=15,ntm_chem=15
      character*8, parameter :: trname(ntm)=(/
     *    'Ox      ','NOx     ','N2O5    ','HNO3    ','H2O2    ',
     *    'CH3OOH  ','HCHO    ','HO2NO2  ','CO      ','CH4     ',
     *    'PAN     ','Isoprene','AlkylNit','Alkenes ','Paraffin'/)
#endif
#endif
#else
#if (defined TRACERS_AEROSOLS_Koch) && (defined TRACERS_HETCHEM)
      integer, parameter :: ntm=23,Ntm_dust=4
      character*8, parameter :: trname(ntm)=(/
     *    'DMS     ','MSA     ','SO2     ','SO4     ','H2O2_s  ',
     *    'seasalt1','seasalt2','BCII    ','BCIA    ','BCB     ',
     *    'OCII    ','OCIA    ','OCB     ','Clay    ','Silt1   ',
     *    'Silt2   ','Silt3   ',
     *    'SO4_s1  ','SO4_s2  ','SO4_d1  ','SO4_d2  ','SO4_d3  ',
     *    'SO4_d4  '/)
#else
#if (defined TRACERS_AEROSOLS_Koch) && (defined TRACERS_DUST)
      integer, parameter :: ntm=17,ntm_dust=4
      character*8, parameter :: trname(ntm)=(/
     *    'DMS     ','MSA     ','SO2     ','SO4     ','H2O2_s  ',
     *    'seasalt1','seasalt2','BCII    ','BCIA    ','BCB     ',
     *    'OCII    ','OCIA    ','OCB     ','Clay    ','Silt1   ',
     *    'Silt2   ','Silt3   '/)
#else
#if (defined TRACERS_AEROSOLS_Koch) && (defined TRACERS_COSMO)
      integer, parameter :: ntm=11
      character*8, parameter :: trname(ntm)=(/
     *    'DMS     ','MSA     ','SO2     ','SO4     ','H2O2_s  ',
     *    'seasalt1','seasalt2',
     *     'Pb210   ','Be7     ','Be10    ','Rn222   '/)
#else
#ifdef TRACERS_COSMO
      integer, parameter :: ntm=4
      character*8, parameter :: trname(ntm)=(/
     *     'Pb210   ','Be7     ','Be10    ','Rn222   '/)
#else
#ifdef TRACERS_AEROSOLS_Koch
      integer, parameter :: ntm=13
      character*8, parameter :: trname(ntm)=(/
     *    'DMS     ','MSA     ','SO2     ','SO4     ','H2O2_s  ',
     *    'seasalt1','seasalt2','BCII    ','BCIA    ','BCB     ',
     *    'OCB     ','OCII    ','OCIA    '/)
c     integer, parameter :: ntm=3
c     character*8, parameter :: trname(ntm)=(/
c    *    'OCII    ','OCIA    ','OCB     '/)
c    *    'DMS     ','MSA     ','SO2     ','SO4     ','H2O2_s  ',
c    *    'OCII    ','OCB     '/)
    
c use these for Rutgers runs
c     integer, parameter :: ntm=4
c     character*8, parameter :: trname(ntm)=(/
c    *    'DMS     ','SO2     ','SO4     ','H2O2_s  '/)
#else
#ifdef TRACERS_OM_SP
      integer, parameter :: ntm=7
      character*8, parameter :: trname(ntm)=(/
     *    'OCI1    ','OCI2    ','OCI3    ','OCA1    ','OCA2    ',
     *    'OCA3    ','OCA4    '/)
#else
#ifdef TRACERS_DUST
!@var Ntm_dust number of dust tracers
      integer,parameter :: ntm=5,ntm_dust=5
      character*8, parameter :: trname(ntm)=(/'Clay    ','Silt1   ',
     &     'Silt2   ','Silt3   ','Silt4   '/)
#else
#if (defined TRACERS_MINERALS) && (defined TRACERS_QUARZHEM)
      INTEGER,PARAMETER :: Ntm=23,Ntm_dust=23,Ntm_min=20,Ntm_quhe=3
      CHARACTER*8,PARAMETER :: TrName(Ntm)=(/
     &     'ClayIlli','ClayKaol','ClaySmec','ClayCalc','ClayQuar',
     &     'Sil1Quar','Sil1Feld','Sil1Calc','Sil1Hema','Sil1Gyps',
     &     'Sil2Quar','Sil2Feld','Sil2Calc','Sil2Hema','Sil2Gyps',
     &     'Sil3Quar','Sil3Feld','Sil3Calc','Sil3Hema','Sil3Gyps',
     &     'Sil1QuHe','Sil2QuHe','Sil3QuHe'/)
#else
#ifdef TRACERS_MINERALS
      INTEGER,PARAMETER :: Ntm=20,Ntm_dust=20,Ntm_min=20
      CHARACTER*8,PARAMETER :: TrName(Ntm)=(/
     &     'ClayIlli','ClayKaol','ClaySmec','ClayCalc','ClayQuar',
     &     'Sil1Quar','Sil1Feld','Sil1Calc','Sil1Hema','Sil1Gyps',
     &     'Sil2Quar','Sil2Feld','Sil2Calc','Sil2Hema','Sil2Gyps',
     &     'Sil3Quar','Sil3Feld','Sil3Calc','Sil3Hema','Sil3Gyps'/)
#else
#ifdef TRACERS_QUARZHEM
      INTEGER,PARAMETER :: Ntm=3,Ntm_dust=3,Ntm_quhe=3
      CHARACTER*8,PARAMETER :: TrName(Ntm)=(/'Sil1QuHe','Sil2QuHe',
     &     'Sil3QuHe'/)
#else
#ifdef TRACERS_WATER
      integer, parameter :: ntm=2
      character*8, parameter :: trname(ntm)=(/'Air     ','Water   '/)
#else
#ifdef TRACERS_OCEAN
      integer, parameter :: ntm=1
      character*8, parameter :: trname(ntm)=(/'Water   '/)
#else ! default for TRACERS_ON
      integer, parameter :: ntm=1
      character*8, parameter :: trname(ntm)=(/'Air     '/)
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif 
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#ifdef regional_Ox_tracers
!@var NregOx number of regional Ox tracers
      integer, parameter :: NregOx=6
!@var regOx_s southern limit of regional Ox tracers (deg)
!@var regOx_n northern limit of regional Ox tracers (deg)
!@var regOx_e eastern limit of regional Ox tracers (deg)
!@var regOx_w western limit of regional Ox tracers (deg)
!@var regOx_t upper (top) limit of regional Ox tracers (hPa)
!@var regOx_b lower (bottom) limit of regional Ox tracers (hPa)
      real*8, dimension(NregOx)::regOx_s,regOx_n,regOx_t,regOx_b,
     &                           regOx_e,regOx_w
#endif
!@var N_XXX: variable names of indices for tracers (init = 0)
      integer :: 
     *     n_Air=0,    n_SF6=0,   n_Rn222=0, n_CO2=0,      n_N2O=0,
     *     n_CFC11=0,  n_14CO2=0, n_CH4=0,   n_O3=0,       n_water=0,
     *     n_H2O18=0,  n_HDO=0,   n_HTO=0,   n_Ox=0,       n_NOx=0, 
     *     n_N2O5=0,   n_HNO3=0,  n_H2O2=0,  n_CH3OOH=0,   n_HCHO=0,
     *     n_HO2NO2=0, n_CO=0,    n_PAN=0,
     *     n_Isoprene=0, n_AlkylNit=0, n_Alkenes=0, n_Paraffin=0,
     *     n_DMS=0,    n_MSA=0,   n_SO2=0,   n_SO4=0,    n_H2O2_s=0,
     *     n_ClOx=0,   n_BrOx=0,  n_HCl=0,   n_HOCl=0,   n_ClONO2=0,
     *     n_HBr=0,    n_HOBr=0,  n_BrONO2=0,n_CFC=0,
     *     n_Pb210 = 0,n_Be7=0,   n_Be10=0,
     *     n_seasalt1=0,  n_seasalt2=0, n_SO4_d1=0,  n_SO4_d2=0,
     *     n_SO4_d3=0, n_SO4_d4=0, n_SO4_s1=0,  n_SO4_s2=0,
     *     n_BCII=0,  n_BCIA=0,  n_BCB=0,
     *     n_OCII=0,  n_OCIA=0,  n_OCB=0,
     *     n_OCI1=0,  n_OCI2=0,  n_OCI3=0,
     *     n_OCA1=0,  n_OCA2=0,  n_OCA3=0, n_OCA4,
     *     n_OxREG1=0,n_OxREG2=0,n_OxREG3=0,
     *     n_OxREG4=0,n_OxREG5=0,n_OxREG6=0,
     &     n_clay=0,   n_silt1=0, n_silt2=0, n_silt3=0, n_silt4=0,
     &     n_clayilli=0,n_claykaol=0,n_claysmec=0,n_claycalc=0,
     &     n_clayquar=0,n_sil1quar=0,n_sil1feld=0,n_sil1calc=0,
     &     n_sil1hema=0,n_sil1gyps=0,n_sil2quar=0,n_sil2feld=0,
     &     n_sil2calc=0,n_sil2hema=0,n_sil2gyps=0,n_sil3quar=0,
     &     n_sil3feld=0,n_sil3calc=0,n_sil3hema=0,n_sil3gyps=0,
     &     n_sil1quhe=0,n_sil2quhe=0,n_sil3quhe=0

!@var 3D on-line radical array for interactive aerosol and gas
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: oh_live
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: no3_live

C****    The following are set in tracer_IC
!@var T_QLIMIT: if t_qlimit=.true. tracer is maintained as positive
      logical, dimension(ntm) :: t_qlimit
!@var trdecay radioactive decay constant (1/s) (=0 for stable tracers)
      real*8, dimension(ntm) :: trdecay
!@dbparam ITIME_TR0: start time for each tracer (hours)
      integer, dimension(ntm) :: itime_tr0
!@var MTRACE: timing index for tracers
      integer mtrace
#if (defined TRACERS_SPECIAL_Shindell) || (defined TRACERS_AEROSOLS_Koch)
!@var MCHEM: timing index for chemistry
      integer mchem
!@var RSULF1, RSULF2, RSULF3, RSULF4: rate coefficients
c for gas phase sulfur chemistry used by aerosol and chemistry models
      real*8, DIMENSION(IM,JM,LM):: rsulf1, rsulf2,rsulf3,rsulf4 
#endif
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_OM_SP)
!@dbparam imAER 0 determines emission choice: 0 present-day, 1 AEROCOM,
c  2 Bond/Streets past, present, future using sector inputs; 3 historic
      integer :: imAER = 1
!@dbparam imPI is 0 for industrial simulations, 1 for pre-industrial
      integer :: imPI = 0
!@dbparam aer_int_yr indicates year of emission
      integer :: aer_int_yr = 0
!@var SNFST0,TNFST0 are instantaneous SW, LW aerosol forcings for AEROCOM
      real*8 SNFST0(2,NTM,IM,JM),TNFST0(2,NTM,IM,JM)
#endif
!@dbparam rnsrc is 0=standard, 1=Conen&Robertson, 2=modified Schery&Wasiolek
      integer :: rnsrc = 0
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
!@dbparam imDUST is 1 for AEROCOM-prescribed simulations, 0 interactive
      INTEGER :: imDUST=0
!@param lim dimension 1 of lookup table for mean surface wind speed integration
!@param ljm dimension 2 of lookup table for mean surface wind speed integration
      INTEGER,PARAMETER :: kim=234,kjm=234
!@var table1 array for lookup table for calculation of mean surface wind speed
!@var x11 index of table1 for GCM surface wind speed from 0 to 50 m/s
!@var x21 index of table1 for sub grid scale velocity scale (sigma)
      REAL*8 :: table1(kim,kjm),x11(kim),x21(kjm)
#endif
#ifdef TRACERS_ON
!@var NTM_POWER: Power of 10 associated with each tracer (for printing)
      integer, dimension(ntm) :: ntm_power
!@var TR_MM: molecular mass of each tracer (g/mole)
      real*8, dimension(ntm) :: tr_mm
!@var needtrs: true if surface tracer value from PBL is required
      logical, dimension(ntm) :: needtrs
#ifdef TRACERS_DRYDEP
!@var dodrydep: true if tracer should undergo dry deposition
      logical, dimension(ntm) :: dodrydep
!@var F0 reactivity factor for oxidation of biological substances
      real*8, dimension(ntm)  :: F0
!@var HSTAR Henry's Law const for tracer dry deposition. mole/(L atm)
!@+   Same as the tr_RKD wet dep variable, except for units.
!@+   If F0 & HSTAR both 0, & tracer not particulate, then no drydep.
      real*8, dimension(ntm)  :: HSTAR
#endif
#endif
C**** Note units for these parameters!
C**** Example: clay dust; trpdens=2.5d3, trradius=0.73d-6 
C****          silt dust; trpdens=2.65d3, trradius=6.1d-6 
C****
!@var trpdens tracer particle density (kg/m^3) 
!@+               (=0 for non-particle tracers)
      real*8, dimension(ntm) :: trpdens
#if (defined TRACERS_MINERALS) || (defined TRACERS_QUARZHEM)
!@param DenQuarz particle density of quartz
!@param DenHema particle density of hematite
      REAL*8,PARAMETER :: DenQuarz=2.62D3,DenHema=5.3D3
#endif
#ifdef TRACERS_QUARZHEM
!@dbparam FreeFe free iron to total iron (free+structural) ratio in minerals
      REAL*8 :: FreeFe=0.5D0
!@dbparam FrHeQu fraction of hematite in quartz/hematite aggregate
      REAL*8 :: FrHeQu=0.1D0
#endif
!@var trradius tracer effective radius (m) (=0 for non particle tracers)
      real*8, dimension(ntm) :: trradius

!@var TRM: Tracer array (kg)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: trm

!@var TRMOM: Second order moments for tracers (kg)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:,:) :: trmom

!@var ntsurfsrcmax maximum number of surface 2D sources/sinks
      integer, parameter :: ntsurfsrcmax=14
!@var ntsurfsrc no. of non-interactive surface sources for each tracer
      integer, dimension(ntm) :: ntsurfsrc 
!@var ntisurfsrc no. of interactive surface sources for each tracer
      integer, dimension(ntm) :: ntisurfsrc 
!@var nt3Dsrcmax maximum number of 3D tracer sources/sinks
      integer, parameter :: nt3Dsrcmax=6

!@param nGAS   index for wetdep tracer type = gas
!@param nPART  index for wetdep tracer type = particle/aerosol
!@param nWATER index for wetdep tracer type = water 
!@+       (original condense method)
      integer, parameter :: nGAS=1, nPART=2, nWATER=3
!@var tr_wd_TYPE: tracer wet dep type (gas, particle, water)
      integer, dimension(ntm) :: tr_wd_TYPE
!@var tr_RKD: Henry's Law coefficient (in mole/Joule please !)
      real*8, dimension(ntm) :: tr_RKD
!@var tr_DHD: coefficient of temperature-dependence term of Henry's
!@+   Law coefficient (in Joule/mole please !)
      real*8, dimension(ntm) :: tr_DHD
!@var fq_aer fraction of aerosol that condenses
      real*8 fq_aer(ntm)
!@param nChemistry index for tracer chemistry 3D source
!@param nStratwrite index for tracer stratosphic overwrite 3D source
!@param nLightning index for tracer lightning 3D source
!@param nAircraft index for tracer aircraft 3D source
      INTEGER, PARAMETER :: nChemistry  = 1, nStratwrite = 2,
     &                      nLightning  = 3, nAircraft   = 4

#ifdef TRACERS_WATER   
!@param nWD_TYPES number of tracer types for wetdep purposes  
      integer, parameter :: nWD_TYPES=3 !(gas,particle,water)
!@param tr_evap_fact fraction of re-evaporation by tracer type
C note, tr_evap_fact is not dimensioned as NTM:
      REAL*8, parameter, dimension(nWD_TYPES) :: tr_evap_fact=
     *     (/1.d0, 0.5d0,  1.d0/)
!@var tr_H2ObyCH4 conc. of tracer in water from methane oxidation 
      real*8, dimension(ntm) :: tr_H2ObyCH4
!@var dowetdep true if tracer has some form of wet deposition
      logical, dimension(ntm) :: dowetdep
!@var TRWM tracer in cloud liquid water amount (kg)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: trwm
#ifdef TRACERS_SPECIAL_O18
!@dbparam supsatfac factor controlling super saturation for isotopes
      real*8 :: supsatfac = 2d-3
#endif
#endif

#if (defined TRACERS_WATER) || (defined TRACERS_OCEAN)
!@var TRW0 default tracer concentration in water (kg/kg)
      real*8, dimension(ntm) :: trw0
!@var NTROCN scaling power factor for ocean/ice tracer concentrations
      integer, dimension(ntm) :: ntrocn
#endif

#ifdef TRACERS_OCEAN
!@var TRGLAC tracer ratio in glacial runoff to ocean (kg/kg)
      real*8, dimension(ntm) :: trglac
#endif
#ifdef TRACERS_HETCHEM
      integer, parameter :: rhet=3
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: rxts,rxts1,rxts2,rxts3
     *                                         ,rxts4,rxtss1,rxtss2
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: krate
#endif

      END MODULE TRACER_COM

      SUBROUTINE ALLOC_TRACER_COM(grid)
!@sum  To allocate arrays whose sizes now need to be determined at
!@+    run time
!@auth NCCS (Goddard) Development Team
!@ver  1.0
      USE TRACER_COM
      USE DOMAIN_DECOMP, ONLY : DIST_GRID, GET
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(IN) :: grid

      INTEGER :: J_1H, J_0H
      INTEGER :: IER

C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)

      ALLOCATE(     oh_live(IM,J_0H:J_1H,LM),
     *             no3_live(IM,J_0H:J_1H,LM),
     *                  trm(IM,J_0H:J_1H,LM,NTM),
     *                trmom(NMOM,IM,J_0H:J_1H,LM,NTM) )
  
#ifdef TRACERS_WATER
      ALLOCATE(        trwm(IM,J_0H:J_1H,LM,NTM) )
#endif
#ifdef TRACERS_HETCHEM
      ALLOCATE( rxts(IM,J_0H:J_1H,LM),rxts1(IM,J_0H:J_1H,LM),
     *         rxtss1(IM,J_0H:J_1H,LM),rxtss2(IM,J_0H:J_1H,LM),
     *         rxts2(IM,J_0H:J_1H,LM),rxts3(IM,J_0H:J_1H,LM),
     *         rxts4(IM,J_0H:J_1H,LM),krate(IM,J_0H:J_1H,LM,rhet) )
#endif

      END SUBROUTINE ALLOC_TRACER_COM

