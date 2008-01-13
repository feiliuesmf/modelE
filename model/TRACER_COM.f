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
      integer, parameter :: ntm=4
      character*8, parameter :: trname(ntm)=(/
     *     'Water   ','H2O18   ','HDO     ', 'H2O17   '/)  !,'HTO     '/)
#else
#ifdef TRACERS_GASEXCH_Natassa
#ifdef TRACERS_GASEXCH_CFC_Natassa
      integer, parameter :: ntm=1
      character*8, parameter :: trname(ntm)=(/'CFCn    '/)
#endif  /* TRACERS_GASEXCH_CFC_Natassa */
#ifdef TRACERS_GASEXCH_CO2_Natassa
      integer, parameter :: ntm=1
      character*8, parameter :: trname(ntm)=(/'CO2n    '/)
#endif  /* TRACERS_GASEXCH_CO2_Natassa */
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
#if (defined TRACERS_DUST) && (defined SHINDELL_STRAT_CHEM) &&\
    (defined TRACERS_AEROSOLS_Koch) && (defined TRACERS_HETCHEM) &&\
    (defined TRACERS_NITRATE)
!@var ntm_chem number of drew-only tracers

      integer, parameter :: ntm=51,ntm_chem=25,ntm_dust=4
      character*8, parameter :: trname(ntm)=(/
     *    'Ox      ','NOx     ','ClOx    ','BrOx    ','N2O5    ',
     *    'HNO3    ','H2O2    ','CH3OOH  ','HCHO    ','HO2NO2  ',
     *    'CO      ','CH4     ','PAN     ','Isoprene','AlkylNit',
     *    'Alkenes ','Paraffin','HCl     ','HOCl    ','ClONO2  ',
     *    'HBr     ','HOBr    ','BrONO2  ','N2O     ','CFC     ',
     *    'DMS     ','MSA     ','SO2     ','SO4     ','H2O2_s  ',
     *    'seasalt1','seasalt2','BCII    ','BCIA    ','BCB     ',
     *    'OCII    ','OCIA    ','OCB     ',
     *    'Clay    ','Silt1   ','Silt2   ','Silt3   ','NH3     ',
     *    'NH4     ','NO3p    ','SO4_d1  ','SO4_d2  ','SO4_d3  ',
     *    'N_d1    ','N_d2    ','N_d3    '/)
#else
#if (defined TRACERS_DUST) && (defined SHINDELL_STRAT_CHEM) &&\
    (defined TRACERS_AEROSOLS_Koch) && (defined TRACERS_NITRATE)
!@var ntm_chem number of drew-only tracers

      integer, parameter :: ntm=45,ntm_chem=25,ntm_dust=4
      character*8, parameter :: trname(ntm)=(/
     *    'Ox      ','NOx     ','ClOx    ','BrOx    ','N2O5    ',
     *    'HNO3    ','H2O2    ','CH3OOH  ','HCHO    ','HO2NO2  ',
     *    'CO      ','CH4     ','PAN     ','Isoprene','AlkylNit',
     *    'Alkenes ','Paraffin','HCl     ','HOCl    ','ClONO2  ',
     *    'HBr     ','HOBr    ','BrONO2  ','N2O     ','CFC     ',
     *    'DMS     ','MSA     ','SO2     ','SO4     ','H2O2_s  ',
     *    'seasalt1','seasalt2','BCII    ','BCIA    ','BCB     ',
     *    'OCII    ','OCIA    ','OCB     ',
     *    'Clay    ','Silt1   ','Silt2   ','Silt3   ','NH3     ',
     *    'NH4     ','NO3p    '/)
#else
#if (defined SHINDELL_STRAT_CHEM) && (defined TRACERS_AEROSOLS_Koch)
!@var ntm_chem number of drew-only tracers

#ifdef SULF_ONLY_AEROSOLS
      integer, parameter :: ntm=30,ntm_chem=25
#else
      integer, parameter :: ntm=38,ntm_chem=25
#endif /* SULF_ONLY_AEROSOLS */
      character*8, parameter :: trname(ntm)=(/
     *    'Ox      ','NOx     ','ClOx    ','BrOx    ','N2O5    ',
     *    'HNO3    ','H2O2    ','CH3OOH  ','HCHO    ','HO2NO2  ',
     *    'CO      ','CH4     ','PAN     ','Isoprene','AlkylNit',
     *    'Alkenes ','Paraffin','HCl     ','HOCl    ','ClONO2  ',
     *    'HBr     ','HOBr    ','BrONO2  ','N2O     ','CFC     ',
     *    'DMS     ','MSA     ','SO2     ','SO4     ','H2O2_s  '
#ifdef SULF_ONLY_AEROSOLS
     *    /)
#else
     *   ,'seasalt1','seasalt2','BCII    ','BCIA    ','BCB     ',
     *    'OCII    ','OCIA    ','OCB     '/)
#endif /* SULF_ONLY_AEROSOLS */
#else
#if (defined TRACERS_DUST) && (defined TRACERS_SPECIAL_Shindell) &&\
    (defined TRACERS_AEROSOLS_Koch) && (defined TRACERS_HETCHEM) &&\
    (defined TRACERS_NITRATE)
!@var ntm_chem number of drew-only tracers
      integer, parameter :: ntm=41,ntm_chem=15,ntm_dust=4
      character*8, parameter :: trname(ntm)=(/
     *    'Ox      ','NOx     ','N2O5    ','HNO3    ','H2O2    ',
     *    'CH3OOH  ','HCHO    ','HO2NO2  ','CO      ','CH4     ',
     *    'PAN     ','Isoprene','AlkylNit','Alkenes ','Paraffin',
     *    'DMS     ','MSA     ','SO2     ','SO4     ','H2O2_s  ',
     *    'seasalt1','seasalt2','BCII    ','BCIA    ','BCB     ',
     *    'OCII    ','OCIA    ','OCB     ',
     *    'Clay    ','Silt1   ','Silt2   ','Silt3   ','NH3     ',
     *    'NH4     ','NO3p    ','SO4_d1  ','SO4_d2  ','SO4_d3  ',
     *    'N_d1    ','N_d2    ','N_d3    '/)
#else
#if (defined TRACERS_DUST) && (defined TRACERS_SPECIAL_Shindell) &&\
    (defined TRACERS_AEROSOLS_Koch) && (defined TRACERS_HETCHEM)
!@var ntm_chem number of drew-only tracers
      integer, parameter :: ntm=38,ntm_chem=15,ntm_dust=4
      character*8, parameter :: trname(ntm)=(/
     *    'Ox      ','NOx     ','N2O5    ','HNO3    ','H2O2    ',
     *    'CH3OOH  ','HCHO    ','HO2NO2  ','CO      ','CH4     ',
     *    'PAN     ','Isoprene','AlkylNit','Alkenes ','Paraffin',
     *    'DMS     ','MSA     ','SO2     ','SO4     ','H2O2_s  ',
     *    'seasalt1','seasalt2','BCII    ','BCIA    ','BCB     ',
     *    'OCII    ','OCIA    ','OCB     ','Clay    ','Silt1   ',
     *    'Silt2   ','Silt3   ','SO4_d1  ','SO4_d2  ','SO4_d3  ',
     *    'N_d1    ','N_d2    ','N_d3    '/)
#else
#if (defined TRACERS_SPECIAL_Shindell) && (defined TRACERS_AEROSOLS_Koch) &&\
    (defined TRACERS_HETCHEM)
!@var ntm_chem number of drew-only tracers
      integer, parameter :: ntm=34,ntm_chem=15
      character*8, parameter :: trname(ntm)=(/
     *    'Ox      ','NOx     ','N2O5    ','HNO3    ','H2O2    ',
     *    'CH3OOH  ','HCHO    ','HO2NO2  ','CO      ','CH4     ',
     *    'PAN     ','Isoprene','AlkylNit','Alkenes ','Paraffin',
     *    'DMS     ','MSA     ','SO2     ','SO4     ','H2O2_s  ',
     *    'seasalt1','seasalt2','BCII    ','BCIA    ','BCB     ',
     *    'OCII    ','OCIA    ','OCB     ',
     *    'SO4_d1  ','SO4_d2  ','SO4_d3  ',
     *    'N_d1    ','N_d2    ','N_d3    '/)
#else
#if (defined TRACERS_DUST) && (defined TRACERS_SPECIAL_Shindell) &&\
    (defined TRACERS_AEROSOLS_Koch) && (defined TRACERS_NITRATE)
!@var ntm_chem number of drew-only tracers
      integer, parameter :: ntm=35,ntm_chem=15,ntm_dust=4
      character*8, parameter :: trname(ntm)=(/
     *    'Ox      ','NOx     ','N2O5    ','HNO3    ','H2O2    ',
     *    'CH3OOH  ','HCHO    ','HO2NO2  ','CO      ','CH4     ',
     *    'PAN     ','Isoprene','AlkylNit','Alkenes ','Paraffin',
     *    'DMS     ','MSA     ','SO2     ','SO4     ','H2O2_s  ',
     *    'seasalt1','seasalt2','BCII    ','BCIA    ','BCB     ',
     *    'OCII    ','OCIA    ','OCB     ','NH3     ','NH4     ',
     *    'NO3p    ','Clay    ','Silt1   ','Silt2   ','Silt3   '/)
#else
#if (defined TRACERS_DUST) && (defined TRACERS_AEROSOLS_Koch) &&\
    (defined TRACERS_NITRATE)
      integer, parameter :: ntm=20,ntm_dust=4
      character*8, parameter :: trname(ntm)=(/
     *    'DMS     ','MSA     ','SO2     ','SO4     ','H2O2_s  ',
     *    'seasalt1','seasalt2','BCII    ','BCIA    ','BCB     ',
     *    'OCII    ','OCIA    ','OCB     ','NH3     ','NH4     ',
     *    'NO3p    ','Clay    ','Silt1   ','Silt2   ','Silt3   '/)
#else
#if (defined TRACERS_AEROSOLS_Koch) && (defined TRACERS_NITRATE)
      integer, parameter :: ntm=13,ntm_dust=4
      character*8, parameter :: trname(ntm)=(/
     *    'DMS     ','SO2     ','SO4     ','H2O2_s  ',
     *    'BCII    ','BCIA    ','BCB     ',
     *    'OCII    ','OCIA    ','OCB     ','NH3     ','NH4     ',
     *    'NO3p    '/)
#else
#if (defined TRACERS_DUST) && (defined TRACERS_SPECIAL_Shindell) &&\
    (defined TRACERS_AEROSOLS_Koch)
!@var ntm_chem number of drew-only tracers
      integer, parameter :: ntm=26,ntm_chem=15,ntm_dust=4
      character*8, parameter :: trname(ntm)=(/
     *    'Ox      ','NOx     ','N2O5    ','HNO3    ','H2O2    ',
     *    'CH3OOH  ','HCHO    ','HO2NO2  ','CO      ','CH4     ',
     *    'PAN     ','Isoprene','AlkylNit','Alkenes ','Paraffin',
     *    'DMS     ','MSA     ','SO2     ','SO4     ','H2O2_s  ',
     *    'seasalt1','seasalt2',
     *    'Clay    ','Silt1   ','Silt2   ','Silt3   '/)
#else
#if (defined TRACERS_DUST) && (defined TRACERS_SPECIAL_Shindell)
!@var ntm_chem number of drew-only tracers
      integer, parameter :: ntm=19,ntm_chem=15,ntm_dust=4
      character*8, parameter :: trname(ntm)=(/
     *    'Ox      ','NOx     ','N2O5    ','HNO3    ','H2O2    ',
     *    'CH3OOH  ','HCHO    ','HO2NO2  ','CO      ','CH4     ',
     *    'PAN     ','Isoprene','AlkylNit','Alkenes ','Paraffin',
     *    'Clay    ','Silt1   ','Silt2   ','Silt3   '/)
#else
#if (defined TRACERS_AEROSOLS_Koch) && (defined TRACERS_SPECIAL_Shindell)
!@var ntm_chem number of drew-only tracers
      integer, parameter :: ntm=20,ntm_chem=15
      character*8, parameter :: trname(ntm)=(/
     *    'Ox      ','NOx     ','N2O5    ','HNO3    ','H2O2    ',
     *    'CH3OOH  ','HCHO    ','HO2NO2  ','CO      ','CH4     ',
     *    'PAN     ','Isoprene','AlkylNit','Alkenes ','Paraffin',
     *    'DMS     ','MSA     ','SO2     ','SO4     ','H2O2_s  '/)
#else
#ifdef TRACERS_SPECIAL_Shindell
#ifdef SHINDELL_STRAT_CHEM
#ifdef SHINDELL_STRAT_EXTRA
      integer, parameter :: ntm=26,ntm_chem=25
CCC   integer, parameter :: ntm=28,ntm_chem=25
      character*8, parameter :: trname(ntm)=(/
     *    'Ox      ','NOx     ','ClOx    ','BrOx    ','N2O5    ',
     *    'HNO3    ','H2O2    ','CH3OOH  ','HCHO    ','HO2NO2  ',
     *    'CO      ','CH4     ','PAN     ','Isoprene','AlkylNit',
     *    'Alkenes ','Paraffin','HCl     ','HOCl    ','ClONO2  ',
     *    'HBr     ','HOBr    ','BrONO2  ','N2O     ','CFC     ',
     &                          'GLT     '/)
CCC  &    'Be7     ','Be10    ','GLT     '/)
#else  /* STRAT_CHEM but not STRAT_EXTRA */
      integer, parameter :: ntm=26,ntm_chem=25
      character*8, parameter :: trname(ntm)=(/
     *    'Ox      ','NOx     ','ClOx    ','BrOx    ','N2O5    ',
     *    'HNO3    ','H2O2    ','CH3OOH  ','HCHO    ','HO2NO2  ',
     *    'CO      ','CH4     ','PAN     ','Isoprene','AlkylNit',
     *    'Alkenes ','Paraffin','HCl     ','HOCl    ','ClONO2  ',
     *    'HBr     ','HOBr    ','BrONO2  ','N2O     ','CFC     ',
     *    'Water   '/)
#endif   /* not STRAT_EXTRA */
#else    /* SPECIAL_Shindell but not SHINDELL_STRAT_CHEM */
      integer, parameter :: ntm=15,ntm_chem=15
      character*8, parameter :: trname(ntm)=(/
     *    'Ox      ','NOx     ','N2O5    ','HNO3    ','H2O2    ',
     *    'CH3OOH  ','HCHO    ','HO2NO2  ','CO      ','CH4     ',
     *    'PAN     ','Isoprene','AlkylNit','Alkenes ','Paraffin'/)
#endif  /* not SHINDELL_STRAT_CHEM */
#else
#if (defined TRACERS_AEROSOLS_Koch) && (defined TRACERS_HETCHEM) &&\
    (defined TRACERS_DUST)
      integer, parameter :: ntm=20,ntm_dust=4
      character*8, parameter :: trname(ntm)=(/
     *    'DMS     ','MSA     ','SO2     ','SO4     ','H2O2_s  ',
     *    'seasalt1','seasalt2','BCII    ','BCIA    ','BCB     ',
     *    'OCII    ','OCIA    ','OCB     ','Clay    ','Silt1   ',
     *    'Silt2   ','Silt3   ','SO4_d1  ','SO4_d2  ','SO4_d3  '/)
#else
#if (defined TRACERS_AEROSOLS_Koch) && (defined TRACERS_DUST)
      integer, parameter :: ntm=17,ntm_dust=4
      character*8, parameter :: trname(ntm)=(/
     *    'DMS     ','MSA     ','SO2     ','SO4     ','H2O2_s  ',
     *    'seasalt1','seasalt2','BCII    ','BCIA    ','BCB     ',
     *    'OCII    ','OCIA    ','OCB     ','Clay    ','Silt1   ',
     *    'Silt2   ','Silt3   '/)
#else
#if (defined TRACERS_AEROSOLS_Koch) && (defined TRACERS_COSMO) &&\
    (defined TRACERS_RADON)
      integer, parameter :: ntm=11
      character*8, parameter :: trname(ntm)=(/
     *    'DMS     ','MSA     ','SO2     ','SO4     ','H2O2_s  ',
     *    'seasalt1','seasalt2',
     *     'Pb210   ','Be7     ','Be10    ','Rn222   '/)
#else
#if (defined TRACERS_COSMO)
#if (defined TRACERS_RADON)
      integer, parameter :: ntm=4
      character*8, parameter :: trname(ntm)=(/
     *     'Pb210   ','Be7     ','Be10    ','Rn222   '/)
#else /* COSMO w/o RADON */
      integer, parameter :: ntm=2
      character*8, parameter :: trname(ntm)=(/
     *     'Be7     ','Be10    '/)
#endif
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
#ifdef TRACERS_AMP
#ifdef TRACERS_AMP_M1
! The order of the ntmAMP aerosols matters!!!
      integer, parameter :: ntm=56,ntmAMP=51,ntm_dust=4
      character*8, parameter :: trname(ntm)=(/
     *    'M_NO3   ','M_NH4   ','M_H2O   ','M_AKK_SU','N_AKK_1 ',!AKK
     *    'M_ACC_SU','N_ACC_1 ','M_DD1_SU','M_DD1_DU','N_DD1_1 ',!ACC,DD1
     *    'M_DS1_SU','M_DS1_DU','N_DS1_1 ','M_DD2_SU','M_DD2_DU',!DS1,DD2
     *    'N_DD2_1 ','M_DS2_SU','M_DS2_DU','N_DS2_1 ','M_SSA_SU',!DD2,DS2,SSA
     *    'M_SSA_SS','M_SSC_SS'                                 ,!SSA,SSC
     *    'M_OCC_SU','M_OCC_OC','N_OCC_1 ','M_BC1_SU','M_BC1_BC',!OCC,BC1
     *    'N_BC1_1 ','M_BC2_SU','M_BC2_BC','N_BC2_1 ','M_BC3_SU',!BC1,BC2,BC3
     *    'M_BC3_BC','N_BC3_1 ','M_DBC_SU','M_DBC_BC','M_DBC_DU',!BC3,DBC
     *    'N_DBC_1 ','M_BOC_SU','M_BOC_BC','M_BOC_OC','N_BOC_1 ',!DBC,BOC
     *    'M_BCS_SU','M_BCS_BC','N_BCS_1 ','M_MXX_SU','M_MXX_BC',!BCS,MXX
     *    'M_MXX_OC','M_MXX_DU','M_MXX_SS','N_MXX_1 ','H2SO4   ',
     *    'DMS     ','SO2     ','H2O2_s  ','NH3     '/)
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
#endif
#ifdef TRACERS_AMP_M2
      integer, parameter :: ntm=56,ntmAMP=51,ntm_dust=4
      character*8, parameter :: trname(ntm)=(/
     *    'M_NO3   ','M_NH4   ','M_H2O   ','M_AKK_SU','N_AKK_1 ',
     *    'M_ACC_SU','N_ACC_1 ','M_DD1_SU','M_DD1_DU','N_DD1_1 ',
     *    'M_DS1_SU','M_DS1_DU','N_DS1_1 ','M_DD2_SU','M_DD2_DU',
     *    'N_DD2_1 ','M_DS2_SU','M_DS2_DU','N_DS2_1 ','M_SSA_SU',
     *    'M_SSA_SS','M_SSC_SS',
     *    'M_OCC_SU','M_OCC_OC','N_OCC_1 ','M_BC1_SU','M_BC1_BC',
     *    'N_BC1_1 ','M_BC2_SU','M_BC2_BC','N_BC2_1 ','M_OCS_SU',
     *    'M_OCS_OC','N_OCS_1 ','M_DBC_SU','M_DBC_BC','M_DBC_DU',
     *    'N_DBC_1 ','M_BOC_SU','M_BOC_BC','M_BOC_OC','N_BOC_1 ',
     *    'M_BCS_SU','M_BCS_BC','N_BCS_1 ','M_MXX_SU','M_MXX_BC',
     *    'M_MXX_OC','M_MXX_DU','M_MXX_SS','N_MXX_1 ','H2SO4   ',
     *    'DMS     ','SO2     ','H2O2_s  ','NH3     '/)
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
     *    0 ,12,0 ,0 ,0 , !BC3,DBC
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
#endif
#ifdef TRACERS_AMP_M3
      integer, parameter :: ntm=46,ntmAMP=41,ntm_dust=4
      character*8, parameter :: trname(ntm)=(/
     *    'M_NO3   ','M_NH4   ','M_H2O   ','M_AKK_SU','N_AKK_1 ',
     *    'M_ACC_SU','N_ACC_1 ','M_DD1_SU','M_DD1_DU','N_DD1_1 ',
     *    'M_DS1_SU','M_DS1_DU','N_DS1_1 ','M_DD2_SU','M_DD2_DU',
     *    'N_DD2_1 ','M_DS2_SU','M_DS2_DU','N_DS2_1 ','M_SSA_SU',
     *    'M_SSA_SS','M_SSC_SS',
     *    'M_OCC_SU','M_OCC_OC','N_OCC_1 ','M_BC1_SU','M_BC1_BC',
     *    'N_BC1_1 ','M_BC2_SU','M_BC2_BC','N_BC2_1 ','M_BOC_SU',
     *    'M_BOC_BC','M_BOC_OC','N_BOC_1 ','M_MXX_SU','M_MXX_BC',
     *    'M_MXX_OC','M_MXX_DU','M_MXX_SS','N_MXX_1 ','H2SO4   ',
     *    'DMS     ','SO2     ','H2O2_s  ','NH3     '/)
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
#endif
#ifdef TRACERS_AMP_M4
      integer, parameter :: ntm=39,ntmAMP=34,ntm_dust=4
      character*8, parameter :: trname(ntm)=(/
     *    'M_NO3   ','M_NH4   ','M_H2O   ','M_ACC_SU','N_ACC_1 ',
     *    'M_DD1_SU','M_DD1_DU','N_DD1_1 ','M_DS1_SU','M_DS1_DU',
     *    'N_DS1_1 ','M_DD2_SU','M_DD2_DU','N_DD2_1 ','M_DS2_SU',
     *    'M_DS2_DU','N_DS2_1 ','M_SSS_SU','M_SSS_SS',
     *    'M_OCC_SU','M_OCC_OC','N_OCC_1 ','M_BC1_SU','M_BC1_BC',
     *    'N_BC1_1 ','M_BC2_SU','M_BC2_BC','N_BC2_1 ','M_MXX_SU',
     *    'M_MXX_BC','M_MXX_OC','M_MXX_DU','M_MXX_SS','N_MXX_1 ',
     *    'H2SO4   ','DMS     ','SO2     ','H2O2_s  ','NH3     '/)
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

#endif
#ifdef TRACERS_AMP_M5
      integer, parameter :: ntm=50,ntmAMP=45,ntm_dust=4
      character*8, parameter :: trname(ntm)=(/
     *    'M_NO3   ','M_NH4   ','M_H2O   ','M_AKK_SU','N_AKK_1 ',!AKK
     *    'M_ACC_SU','N_ACC_1 ','M_DD1_SU','M_DD1_DU','N_DD1_1 ',!ACC,DD1
     *    'M_DS1_SU','M_DS1_DU','N_DS1_1 ',                      !DS1
     *    'M_SSA_SU',                                            !SSA
     *    'M_SSA_SS','M_SSC_SS',                                 !SSA,SSC
     *    'M_OCC_SU','M_OCC_OC','N_OCC_1 ','M_BC1_SU','M_BC1_BC',!OCC,BC1
     *    'N_BC1_1 ','M_BC2_SU','M_BC2_BC','N_BC2_1 ','M_BC3_SU',!BC1,BC2,BC3
     *    'M_BC3_BC','N_BC3_1 ','M_DBC_SU','M_DBC_BC','M_DBC_DU',!BC3,DBC
     *    'N_DBC_1 ','M_BOC_SU','M_BOC_BC','M_BOC_OC','N_BOC_1 ',!DBC,BOC
     *    'M_BCS_SU','M_BCS_BC','N_BCS_1 ','M_MXX_SU','M_MXX_BC',!BCS,MXX
     *    'M_MXX_OC','M_MXX_DU','M_MXX_SS','N_MXX_1 ','H2SO4   ',
     *    'DMS     ','SO2     ','H2O2_s  ','NH3     '/)
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
      integer, parameter :: ntm=50,ntmAMP=45,ntm_dust=4
      character*8, parameter :: trname(ntm)=(/
     *    'M_NO3   ','M_NH4   ','M_H2O   ','M_AKK_SU','N_AKK_1 ',!AKK
     *    'M_ACC_SU','N_ACC_1 ','M_DD1_SU','M_DD1_DU','N_DD1_1 ',!ACC,DD1
     *    'M_DS1_SU','M_DS1_DU','N_DS1_1 ',                      !DS1
     *    'M_SSA_SU',                                            !SSA
     *    'M_SSA_SS','M_SSC_SS',                                 !SSA,SSC
     *    'M_OCC_SU','M_OCC_OC','N_OCC_1 ','M_BC1_SU','M_BC1_BC',!OCC,BC1
     *    'N_BC1_1 ','M_BC2_SU','M_BC2_BC','N_BC2_1 ','M_OCS_SU',!BC1,BC2,OCS
     *    'M_OCS_OC','N_OCS_1 ','M_DBC_SU','M_DBC_BC','M_DBC_DU',!OCS,DBC
     *    'N_DBC_1 ','M_BOC_SU','M_BOC_BC','M_BOC_OC','N_BOC_1 ',!DBC,BOC
     *    'M_BCS_SU','M_BCS_BC','N_BCS_1 ','M_MXX_SU','M_MXX_BC',!BCS,MXX
     *    'M_MXX_OC','M_MXX_DU','M_MXX_SS','N_MXX_1 ','H2SO4   ',
     *    'DMS     ','SO2     ','H2O2_s  ','NH3     '/)
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
      integer, parameter :: ntm=40,ntmAMP=35,ntm_dust=4
      character*8, parameter :: trname(ntm)=(/
     *    'M_NO3   ','M_NH4   ','M_H2O   ','M_AKK_SU','N_AKK_1 ',!AKK
     *    'M_ACC_SU','N_ACC_1 ','M_DD1_SU','M_DD1_DU','N_DD1_1 ',!ACC,DD1
     *    'M_DS1_SU','M_DS1_DU','N_DS1_1 ',                      !DS1
     *    'M_SSA_SU',                                            !SSA
     *    'M_SSA_SS','M_SSC_SS',                                 !SSA,SSC
     *    'M_OCC_SU','M_OCC_OC','N_OCC_1 ','M_BC1_SU','M_BC1_BC',!OCC,BC1
     *    'N_BC1_1 ','M_BC2_SU','M_BC2_BC','N_BC2_1 ',           !BC1,BC2
     *    'M_BOC_SU','M_BOC_BC','M_BOC_OC','N_BOC_1 ',           !BOC
     *    'M_MXX_SU','M_MXX_BC',                                 !MXX
     *    'M_MXX_OC','M_MXX_DU','M_MXX_SS','N_MXX_1 ','H2SO4   ',
     *    'DMS     ','SO2     ','H2O2_s  ','NH3     '/)
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
      integer, parameter :: ntm=33,ntmAMP=28,ntm_dust=4
      character*8, parameter :: trname(ntm)=(/
     *    'M_NO3   ','M_NH4   ','M_H2O   ',
     *    'M_ACC_SU','N_ACC_1 ','M_DD1_SU','M_DD1_DU','N_DD1_1 ',!ACC,DD1
     *    'M_DS1_SU','M_DS1_DU','N_DS1_1 ',                      !DS1
     *    'M_SSS_SU',                                            !SSS
     *    'M_SSS_SS',                                            !SSS
     *    'M_OCC_SU','M_OCC_OC','N_OCC_1 ','M_BC1_SU','M_BC1_BC',!OCC,BC1
     *    'N_BC1_1 ','M_BC2_SU','M_BC2_BC','N_BC2_1 ',           !BC1,BC2
     *    'M_MXX_SU','M_MXX_BC',                                 !MXX
     *    'M_MXX_OC','M_MXX_DU','M_MXX_SS','N_MXX_1 ','H2SO4   ',
     *    'DMS     ','SO2     ','H2O2_s  ','NH3     '/)
      integer, parameter :: AMP_MODES_MAP(ntmAMP)=(/
     *    0 ,0 ,0 ,       !
     *    1 ,1 ,2 ,2 ,2,  !ACC,DD1
     *    3 ,3 ,3 ,       !DS1
     *    4 ,             !SSS
     *    4 ,5 ,          !SSS
     *    5 ,5 ,5, 6, 6,  !OCC,BC1
     *    6, 7, 7, 7, 8,  !BC1,BC2
     *    8,8,8,8,8,/)    !MXX
      integer, parameter :: AMP_NUMB_MAP(ntmAMP)=(/
     *    0 ,0 ,0 ,       !
     *    0 ,1 ,0 ,0 ,2,  !ACC,DD1
     *    0 ,0 ,3 ,       !DS1
     *    0 ,             !SSS
     *    0 ,0 ,          !SSS
     *    0 ,0 ,5, 0, 0,  !OCC,BC1
     *    6, 0, 0, 7, 0,  !BC1,BC2
     *    0,0,0,0,8,/)    !MXX
      integer, parameter :: AMP_AERO_MAP(ntmAMP)=(/
     *    1 ,2 ,3 ,4 ,5 ,
     *    6 ,7 ,8 ,9 ,10,
     *    11,12,13,15,16,
     *    17,18,19,20,21,
     *    22,23,24,25,26,
     *    27,28,29   /)
#endif
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
      INTEGER,PARAMETER :: ntm=23,ntm_dust=23,ntm_min=20,ntm_quhe=3
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
#ifdef TRACERS_AGE_OCEAN
      character*8, parameter :: trname(ntm)=(/'Age     '/)
#else
      character*8, parameter :: trname(ntm)=(/'Water   '/)
#endif   /* TRACERS_AGE_OCEAN */
#else /* default for TRACERS_ON */
      integer, parameter :: ntm=1
      character*8, parameter :: trname(ntm)=(/'Air     '/)
#endif   /* TRACERS_OCEAN */
#endif   /* TRACERS_WATER */
#endif   /* TRACERS_QUARZHEM */
#endif   /* TRACERS_MINERALS */
#endif   /* TRACERS_MINERALS + QUARZHEM */
#endif   /* TRACERS_DUST */
#endif   /* TRACERS_AMP */
#endif   /* TRACERS_AEROSOLS_Koch */
#endif   /* TRACERS_COSMO */
#endif   /* TRACERS_AEROSOLS_Koch + COSMO + RADON */
#endif   /* TRACERS_AEROSOLS_Koch + DUST */
#endif   /* TRACERS_AEROSOLS_Koch + HETCHEM + DUST */
#endif   /* TRACERS_SPECIAL_Shindell */
#endif   /* TRACERS_SPECIAL_Shindell + AEROSOLS_Koch */
#endif   /* TRACERS_SPECIAL_Shindell + DUST */
#endif   /* TRACERS_SPECIAL_Shindell + DUST + AEROSOLS_Koch + NITRATE */
#endif   /* TRACERS_SPECIAL_Shindell + DUST + AEROSOLS_Koch */
#endif   /* TRACERS_AEROSOLS_Koch + NITRATE */
#endif   /* TRACERS_DUST + AEROSOLS_Koch + Nitrate */
#endif   /* TRACERS_DUST + SPECIAL_Shindell + AEROSOLS_Koch + NITRATE */
#endif   /* TRACERS_TRACERS_SPECIAL_Shindell + AEROSOLS_Koch + HETCHEM */
#endif   /* DUST + SPECIAL_Shind. + AEROSOLS_Koch + HETCHEM */
#endif   /* DUST + SPECIAL_Shind. + AEROSOLS_Koch + HETCHEM + NITRATE */
#endif   /* SHINDELL_STRAT_CHEM + TRACERS_AEROSOLS_Koch */
#endif   /* DUST + SH._STRAT_CHEM + AEROSOLS_Koch + NITRATE */
#endif   /* DUST + SH._STRAT_CHEM + AEROSOLS_Koch + HETCHEM + NITRATE */
#endif   /* TRACERS_SPECIAL_Lerner */
#endif   /* TRACERS_SPECIAL_Lerner + TRACERS_WATER*/
#endif   /* TRACERS_GASEXCH_Natassa */
#endif   /* TRACERS_SPECIAL_O18 */

!@var N_XXX: variable names of indices for tracers (init = 0)
      integer ::
     *     n_Air=0,    n_SF6=0,   n_Rn222=0, n_CO2=0,      n_N2O=0,
     *     n_CFC11=0,  n_14CO2=0, n_CH4=0,   n_O3=0,       n_water=0,
     *     n_H2O18=0,  n_HDO=0,   n_HTO=0,   n_Ox=0,       n_NOx=0,
     *     n_N2O5=0,   n_HNO3=0,  n_H2O2=0,  n_CH3OOH=0,   n_HCHO=0,
     *     n_HO2NO2=0, n_CO=0,    n_PAN=0,   n_H2O17=0,
     *     n_Isoprene=0, n_AlkylNit=0, n_Alkenes=0, n_Paraffin=0,
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
     *     n_OCI1=0,  n_OCI2=0,  n_OCI3=0,
     *     n_OCA1=0,  n_OCA2=0,  n_OCA3=0, n_OCA4,
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
     *     n_H2SO4=0

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
#if (defined TRACERS_SPECIAL_Shindell) || (defined TRACERS_AEROSOLS_Koch) ||\
    (defined TRACERS_AMP)
!@var MCHEM: timing index for chemistry
      integer mchem
!@var RSULF1, RSULF2, RSULF3, RSULF4: rate coefficients
c for gas phase sulfur chemistry used by aerosol and chemistry models
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:)::rsulf1,rsulf2,rsulf3,rsulf4
#endif
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_OM_SP) ||\
    (defined TRACERS_AMP)
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
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP)
!@dbparam imDUST is 1 for AEROCOM-prescribed simulations, 0 interactive
      INTEGER :: imDUST=0
!@param kim dimension 1 of lookup table for mean surface wind speed integration
!@param kjm dimension 2 of lookup table for mean surface wind speed integration
      INTEGER,PARAMETER :: kim=234,kjm=234
!@var table1 array for lookup table for calculation of mean surface wind speed
!@+          local to each grid box
      REAL*8,ALLOCATABLE,DIMENSION(:,:) :: table1
!@var x11 index of table1 for GCM surface wind speed from 0 to 50 m/s
!@var x21 index of table1 for sub grid scale velocity scale (sigma)
      REAL*8 :: x11(kim),x21(kjm)
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

!@var TRDN1: lowest level downdraft tracer concentration (kg/kg)
       REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: trdn1

!@var ntsurfsrcmax maximum number of surface 2D sources/sinks
      integer, parameter :: ntsurfsrcmax=15
!@var ntsurfsrc no. of non-interactive surface sources for each tracer
      integer, dimension(ntm) :: ntsurfsrc
!@var ntisurfsrc no. of interactive surface sources for each tracer
      integer, dimension(ntm) :: ntisurfsrc
!@var nt3Dsrcmax maximum number of 3D tracer sources/sinks
      integer, parameter :: nt3Dsrcmax=6
!@var sfc_src array holds tracer sources that go into trsource( )
!@+ maybe wasteful of memory, but works for now...
      real*8, allocatable, dimension(:,:,:,:) :: sfc_src
!@var ssname holds source name, read from file header, e.g. to be
!@+ placed into lname and sname arrays.
!@var freq frequency (annual? monthly?) read from emis file header
!@var res horiz. resolution (M? F?) read from emis file header
!@var nameT tracer name read from emis file header (should match trname)
!@var ty_start starting year/decade for a transient emissions file
!@var ty_end     ending year/decade for a transient emissions file
      character*30, dimension(ntm,ntsurfsrcmax) :: ssname ! some maybe
      character*10, dimension(ntm,ntsurfsrcmax) :: nameT  ! need not be
      character*1, dimension(ntm,ntsurfsrcmax) :: freq,res ! arrays
      character*9, dimension(ntm,ntsurfsrcmax) :: Tyears   ! here...
      integer, dimension(ntm,ntsurfsrcmax) :: ty_start,ty_end
!@param kstep number of years assumed between transient emissions slices
      integer, parameter :: kstep=10 ! assume decadal transient emissions
! ---- section for altering tracers sources by sector/region ----
!@param n_max_sect maximum number of sectors for emissions altering
!@param n_max_reg  maximum number of regions for emissions altering
      integer, parameter :: n_max_sect=10, n_max_reg=10
!@var num_tr_sectors number of sectors for a particular tracer and source
      integer, dimension(ntm,ntsurfsrcmax) :: num_tr_sectors
!@var num_tr_sectors3D number of sectors for a tracer's 3D source
      integer, dimension(ntm,nt3Dsrcmax) :: num_tr_sectors3D
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
      integer, dimension(ntm,ntsurfsrcmax,n_max_sect) :: tr_sect_index
!@var tr_sect_index3D holds 3d source sector index for given tracer/source
      integer, dimension(ntm,nt3Dsrcmax,n_max_sect) :: tr_sect_index3D
!@var tr_sect_name array hold the sector name for given tracer/source
      character*10,dimension(ntm,ntsurfsrcmax,n_max_sect):: tr_sect_name
!@var tr_sect_name3D holds 3d source sector name for given tracer/source
      character*10,dimension(ntm,nt3Dsrcmax,n_max_sect):: tr_sect_name3D
!@var sect_name array hold the sector names (all)
      character*10,dimension(n_max_sect):: sect_name
!@var ef_fact the actual factors that alter sources by region/sector
!@var ef_fact3d factors used to alter 3D sources (these are more
!@+ hard-coded for now...)
      real*8, dimension(n_max_sect,n_max_reg) :: ef_fact,ef_fact3D
! variables for outputting a map of the regions:
      real*8, allocatable, dimension(:,:) :: ef_REG_IJ
      real*8, dimension(IM,JM) :: ef_REG_IJ_glob
      real*4, dimension(IM,JM) :: ef_REG_IJ_glob4
! --- end of source-altering section ----------------------------
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
!@var rc_washt aerosol washout rate
      REAL*8 :: rc_washt(Ntm)=1.D-1
!@param nChemistry index for tracer chemistry 3D source
!@param nStratwrite index for tracer stratosphic overwrite 3D source
!@param nOther index for tracer misc. 3D source
!@param nAircraft index for tracer aircraft 3D source
!@param nVolcanic index for tracer volcano 3D source
! Must be a better way to do this, but for now, it is better than
! hardcoding indicies with a number like "3":
      INTEGER, PARAMETER :: nChemistry  = 1, nStratwrite = 2,
     &     nOther  = 3, nAircraft   = 4, nBiomass    = 5,
     &     nVolcanic = 6

#ifdef TRACERS_GASEXCH_Natassa
!@var ocmip_cfc: CFC-11 emissions estimated from OCMIP surf.conc.
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: ocmip_cfc
#endif

#ifdef TRACERS_WATER
!@param nWD_TYPES number of tracer types for wetdep purposes
      integer, parameter :: nWD_TYPES=3 !(gas,particle,water)
!@param tr_evap_fact fraction of re-evaporation of raindrops by tracer type
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
!@var iso_index indexing taking actual tracer number to isotope
!@+   fractionation number (1=water,2=h2o18,3=hdo,4=hto,5=h2o17)
      integer :: iso_index(5)
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
#if (defined TRACERS_HETCHEM) || (defined TRACERS_NITRATE)
      integer, parameter :: rhet=3
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: rxts,rxts1,rxts2,rxts3
     *                                         ,rxts4
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:,:) :: krate
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

      ALLOCATE(   ef_REG_IJ(IM,J_0H:J_1H) )
      ALLOCATE(     oh_live(IM,J_0H:J_1H,LM),
     *             no3_live(IM,J_0H:J_1H,LM),
     *                  trm(IM,J_0H:J_1H,LM,NTM),
     *                trmom(NMOM,IM,J_0H:J_1H,LM,NTM),
     *                trdn1(NTM,IM,J_0H:J_1H),
     *              sfc_src(IM,J_0H:J_1H,ntm,ntsurfsrcmax))

#ifdef TRACERS_GASEXCH_Natassa
      !60years (1939--1998) OCMIP surfc. concentr. converted to
      !global averaged emission rates
      !each value corresponds to the annual value
      ALLOCATE(ocmip_cfc(67,NTM))
#endif

#ifdef TRACERS_WATER
      ALLOCATE(        trwm(IM,J_0H:J_1H,LM,NTM) )
#endif
#ifdef TRACERS_HETCHEM
      ALLOCATE( rxts(IM,J_0H:J_1H,LM),rxts1(IM,J_0H:J_1H,LM),
     *         rxts2(IM,J_0H:J_1H,LM),rxts3(IM,J_0H:J_1H,LM),
     *         rxts4(IM,J_0H:J_1H,LM),krate(IM,J_0H:J_1H,LM,8,rhet))
#endif
#if (defined TRACERS_SPECIAL_Shindell) || (defined TRACERS_AEROSOLS_Koch) ||\
    (defined TRACERS_AMP)
      ALLOCATE(  rsulf1(IM,J_0H:J_1H,LM),rsulf2(IM,J_0H:J_1H,LM),
     *           rsulf3(IM,J_0H:J_1H,LM),rsulf4(IM,J_0H:J_1H,LM) )
#endif

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP)
      ALLOCATE(table1(Kim,Kjm),STAT=ier)
#endif

      END SUBROUTINE ALLOC_TRACER_COM

