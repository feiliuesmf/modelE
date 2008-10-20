#include "rundeck_opts.h"
module TRACERS_SOA
!-------------------------------------------------------------------------------
! SOA formation - created by Kostas Tsigaridis
!-------------------------------------------------------------------------------
!@sum module for the calculation of secondary organic aerosols (SOA). Requires
!@+   tropospheric chemistry and aerosols to be activated.
!@auth Kostas Tsigaridis (ktsigaridis@giss.nasa.gov)
use RESOLUTION, only: LM
use DOMAIN_DECOMP,only: write_parallel
use TRACER_COM, only: ntm,ntm_soa,tr_mm,&
                      n_isoprene,&
                      n_isopp1a,n_isopp2a
implicit none

!@var n_soa_i the first SOA-related species
!@var n_soa_e the last SOA-related species
integer                    :: n_soa_i,n_soa_e
!@param nsoa the total number of aerosol-phase SOA related species (=ntm_soa/2)
integer, parameter         :: nsoa=ntm_soa/2
!@var soagas the amount of gas-phase SOA related species just before chemistry (ug/m3)
!@var soachem the amount of newly (chemically) formed soa mass per timestep (ug/m3)
real*8, dimension(nsoa,LM) :: soagas,soachem
!@var mw the molecular weight of all tracers, in units of tracer mass per mole
!@+      in order to replace the tm_mm which in some exceptional cases is in units of
!@+      a certain atom (typically C or S) per mole.
real*8, dimension(ntm)     :: mw!,agas
!@var apartmass the mass-based yield of semivolatile species from chemistry
!@var apartmolar the molar-based yield of semivolatile species from chemistry
real*8, dimension(nsoa)    :: apartmass,apartmolar
!modelE!#ifdef SOA_LUMPED_NOX_DEP
!modelE!real, dimension(nsoa)    :: apartmass_nox,apartmolar_nox
!modelE!#endif
!@var  molec2ug converts molec/cm3 to ug/m3. 1/molec2ug converts ug/m3 to molec/cm3
real*8, dimension(ntm) :: molec2ug

!
! soa semivolatile products in aerosol phase
!
!@var whichsoa converts species index to soa index
integer, dimension(ntm)  :: whichsoa
!@var issoa converts soa index to species index
integer, dimension(nsoa) :: issoa

!@param soacomp number of different cases in Lambda calculations
integer, parameter                 :: soacomp=10
!@var Lambda empirical factor describing the affinity of species i with species j
!@+          for the activity coefficient calculations. A value close to unity means
!@+          high chemical similarity of the species. A value of -1.0 declares
!@+          symmetric treatment: Lambda(i,j)=Lambda(j,i).
real*8, dimension(soacomp,soacomp) :: Lambda
integer, parameter                 :: imfter=1
integer, parameter                 :: imfaro=2
integer, parameter                 :: imfisop=3
integer, parameter                 :: imfocii=4
integer, parameter                 :: imfocia=5
integer, parameter                 :: imfocb=6
integer, parameter                 :: imfbcii=7
integer, parameter                 :: imfbcia=8
integer, parameter                 :: imfbcb=9
integer, parameter                 :: imfinorg=10
!
!                                   ter     aro    isop    ocii    ocia     ocb    bcii    bcia     bcb   inorg
data Lambda(imfter,1:soacomp)   / 1.00d0, 1.00d0, 1.00d0, 0.80d0, 0.90d0, 0.85d0, 0.70d0, 0.80d0, 0.75d0, 0.70d0/
data Lambda(imfaro,1:soacomp)   /-1.00d0, 1.00d0, 1.00d0, 0.80d0, 0.90d0, 0.85d0, 0.70d0, 0.80d0, 0.75d0, 0.70d0/
data Lambda(imfisop,1:soacomp)  /-1.00d0,-1.00d0, 1.00d0, 0.80d0, 0.90d0, 0.85d0, 0.70d0, 0.80d0, 0.75d0, 0.70d0/
data Lambda(imfocii,1:soacomp)  /-1.00d0,-1.00d0,-1.00d0, 1.00d0, 0.90d0, 0.85d0, 0.70d0, 0.80d0, 0.75d0, 0.70d0/
data Lambda(imfocia,1:soacomp)  /-1.00d0,-1.00d0,-1.00d0,-1.00d0, 1.00d0, 0.85d0, 0.70d0, 0.80d0, 0.75d0, 0.70d0/
data Lambda(imfocb,1:soacomp)   /-1.00d0,-1.00d0,-1.00d0,-1.00d0,-1.00d0, 1.00d0, 0.70d0, 0.80d0, 0.75d0, 0.70d0/
data Lambda(imfbcii,1:soacomp)  /-1.00d0,-1.00d0,-1.00d0,-1.00d0,-1.00d0,-1.00d0, 1.00d0, 0.80d0, 0.75d0, 0.70d0/
data Lambda(imfbcia,1:soacomp)  /-1.00d0,-1.00d0,-1.00d0,-1.00d0,-1.00d0,-1.00d0,-1.00d0, 1.00d0, 0.75d0, 0.70d0/
data Lambda(imfbcb,1:soacomp)   /-1.00d0,-1.00d0,-1.00d0,-1.00d0,-1.00d0,-1.00d0,-1.00d0,-1.00d0, 1.00d0, 0.70d0/
data Lambda(imfinorg,1:soacomp) /-1.00d0,-1.00d0,-1.00d0,-1.00d0,-1.00d0,-1.00d0,-1.00d0,-1.00d0,-1.00d0, 1.00d0/
!
! enthalpies of vaporization (dH, in kJ/mol):
!    0.     No temperature dependence on vapor pressure
!   42.     Chung and Seinfeld, 2002; Kleindienst et al., GRL, 2007 (for isoprene products)
!   72.4    Pun et al., 2003
!   79.     Andersson-Skold and Simpson, 2001
!  109.     for pinic acid; Bilde and Pandis, EST, 2001
!  156.     Strader et al., 1999
!
!@param dH_isoprene enthalpy of vaporization for the isoprene-produced SOA species (KJ/mol)
real*8, parameter :: dH_isoprene=42.d0 ! Chung and Seinfeld, JGR, 2002; Henze and Seinfeld, GRL, 2006; Kleindienst et al., GRL, 2007
!real*8, parameter :: dH_terpenes=72.9d0
!real*8, parameter :: dH_terpenes_pinic=109.d0
!real*8, parameter :: dH_aromatics=72.9d0
!real*8, parameter :: dH_terpenes=42.d0
!real*8, parameter :: dH_terpenes_pinic=42.d0
!real*8, parameter :: dH_aromatics=15.d0

!@var kpart partitioning coefficient of SOA species (m3/ug)
real*8, dimension(LM,nsoa)         :: kpart
!@param kpart_ref partitioning coefficient of SOA species (m3/ug) at the reference temperature kpart_temp_ref
real*8, dimension(nsoa), parameter :: kpart_ref=(/&
! isoprene + OH low NOx SOAb formation, Henze and Seinfeld, GRL 2006
!modelE!#ifdef SOA_FULL
!modelE!0.00862, 1.62,  & ! iisopp1a_nox, iisopp2a_nox
!modelE!#endif
0.00862d0, 1.62d0  & ! iisopp1a_hox, iisopp2a_hox
!modelE!#ifdef SOA_FULL
!modelE!1./15.7, 1./385,& ! iapinp1a_nox, iapinp2a_nox: Presto et al., EST, 2005
!modelE!#endif
!modelE!1./15.7, 1./385,& ! iapinp1a_hox, iapinp2a_hox: Presto et al., EST, 2005
!modelE!#ifndef SOA_MINIMUM
!modelE!#  ifdef SOA_FULL
!modelE!0.195,   0.003, & ! ibpinp1a_nox, ibpinp2a_nox: Griffin et al., JGR, 1999
!modelE!#  endif
!modelE!0.195,   0.003, & ! ibpinp1a_hox, ibpinp2a_hox: Griffin et al., JGR, 1999
!modelE!#endif
!modelE!#ifdef SOA_FULL
!modelE!0.053,   0.0019,& ! itolp1a_nox,  itolp2a_nox: Odum et al., Science, 1997
!modelE!#endif
!modelE!0.053,   0.0019 & ! itolp1a_hox,  itolp2a_hox: Odum et al., Science, 1997
!modelE!#ifndef SOA_MINIMUM
!modelE!#  ifdef SOA_FULL
!modelE!,0.301,   0.008  & ! ixylp1a_nox,  ixylp2a_nox: Song et al., EST, 2005 !!!!!!!!!!!!!!CHANGE ALSO HARDCODED PART IN sources_sinks.f90!!!!!!!!!!!!!!!
!modelE!#  endif
!modelE!,0.229,   0.004  & ! ixylp1a_hox,  ixylp2a_hox: Song et al., EST, 2005 !!!!!!!!!!!!!!CHANGE ALSO HARDCODED PART IN sources_sinks.f90!!!!!!!!!!!!!!!
!modelE!#endif
!modelE!!#ifdef SOA_FULL
!modelE!!0.430,   0.047, & ! itolp1a_nox,  itolp2a_nox: Ng et al., ACP, 2007
!modelE!!#endif
!modelE!!1.e10,   1.e10  & ! itolp1a_hox,  itolp2a_hox: Ng et al., ACP, 2007
!modelE!!#ifndef SOA_MINIMUM
!modelE!!#  ifdef SOA_FULL
!modelE!!,0.761,   0.029  & ! ixylp1a_nox,  ixylp2a_nox: Ng et al., ACP, 2007
!modelE!!#  endif
!modelE!!,1.e10,   1.e10  & ! ixylp1a_hox,  ixylp2a_hox: Ng et al., ACP, 2007
!modelE!!#endif
/)
!modelE!#  ifdef SOA_LUMPED_NOX_DEP
!modelE!real, dimension(nsoa), parameter :: kpart_nox_ref=(/&
!modelE!0.00862, 1.62,  & ! iisopp1a_nox, iisopp2a_nox
!modelE!1./15.7, 1./385,& ! iapinp1a_nox, iapinp2a_nox: Presto et al., EST, 2005
!modelE!0.195,   0.003, & ! ibpinp1a_nox, ibpinp2a_nox: Griffin et al., JGR, 1999
!modelE!0.053,   0.0019,& ! itolp1a_nox,  itolp2a_nox: Odum et al., Science, 1997
!modelE!0.301,   0.008  & ! ixylp1a_nox,  ixylp2a_nox: Song et al., EST, 2005
!modelE!!1.e10,   1.e10, & ! itolp1a_nox,  itolp2a_nox: Ng et al., ACP, 2007
!modelE!!1.e10,   1.e10  & ! ixylp1a_nox,  ixylp2a_nox: Ng et al., ACP, 2007
!modelE!/)
!modelE!#endif

!@param kpart_temp_ref reference temperature where kpart_ref was derived
real*8, dimension(nsoa), parameter :: kpart_temp_ref=(/&
! isoprene + OH low NOx SOAb formation, Henze and Seinfeld, GRL 2006
!modelE!#ifdef SOA_FULL
!modelE!295.,295., & ! iisopp1a_nox, iisopp2a_nox
!modelE!#endif
295.d0,295.d0 & ! iisopp1a_hox, iisopp2a_hox
!modelE!#ifdef SOA_FULL
!modelE!295.,295., & ! iapinp1a_nox, iapinp2a_nox: Presto et al., EST, 2005
!modelE!#endif
!modelE!295.,295., & ! iapinp1a_hox, iapinp2a_hox: Presto et al., EST, 2005
!modelE!#ifndef SOA_MINIMUM
!modelE!#  ifdef SOA_FULL
!modelE!298.,298., & ! ibpinp1a_nox, ibpinp2a_nox: Griffin et al., JGR, 1999
!modelE!#  endif
!modelE!298.,298., & ! ibpinp1a_hox, ibpinp2a_hox: Griffin et al., JGR, 1999
!modelE!#endif
!modelE!#ifdef SOA_FULL
!modelE!298.,298., & ! itolp1a_nox,  itolp2a_nox: Odum et al., Science, 1997
!modelE!#endif
!modelE!298.,298.  & ! itolp1a_hox,  itolp2a_hox: Odum et al., Science, 1997
!modelE!#ifndef SOA_MINIMUM
!modelE!#  ifdef SOA_FULL
!modelE!,300.,300.  & ! ixylp1a_nox,  ixylp2a_nox: Song et al., EST, 2005
!modelE!#  endif
!modelE!,300.,300.  & ! ixylp1a_hox,  ixylp2a_hox: Song et al., EST, 2005
!modelE!#endif
!modelE!!#ifdef SOA_FULL
!modelE!!295.,295., & ! itolp1a_nox,  itolp2a_nox: Ng et al., ACP, 2007
!modelE!!#endif
!modelE!!295.,295.  & ! itolp1a_hox,  itolp2a_hox: Ng et al., ACP, 2007
!modelE!!#ifndef SOA_MINIMUM
!modelE!!#  ifdef SOA_FULL
!modelE!!,295.,295.  & ! ixylp1a_nox,  ixylp2a_nox: Ng et al., ACP, 2007
!modelE!!#  endif
!modelE!!,295.,295.  & ! ixylp1a_hox,  ixylp2a_hox: Ng et al., ACP, 2007
!modelE!!#endif
/)
!modelE!#  ifdef SOA_LUMPED_NOX_DEP
!modelE!real, dimension(nsoa), parameter :: kpart_temp_nox_ref=(/&
!modelE!295.,295., & ! iisopp1a_nox, iisopp2a_nox
!modelE!295.,295., & ! iapinp1a_nox, iapinp2a_nox: Presto et al., EST, 2005
!modelE!298.,298., & ! ibpinp1a_nox, ibpinp2a_nox: Griffin et al., JGR, 1999
!modelE!298.,298., & ! itolp1a_nox,  itolp2a_nox: Odum et al., Science, 1997
!modelE!300.,300.  & ! ixylp1a_nox,  ixylp2a_nox: Song et al., EST, 2005
!modelE!!295.,295., & ! itolp1a_nox,  itolp2a_nox: Ng et al., ACP, 2007
!modelE!!295.,295.  & ! ixylp1a_nox,  ixylp2a_nox: Ng et al., ACP, 2007
!modelE!/)
!modelE!#endif

character(len=300) :: out_line


contains

subroutine soa_init

use TRACER_COM, only: n_bcii,n_bcia,n_bcb,n_ocii,n_ocia,n_ocb
use CONSTANT, only: avog
implicit none

integer   :: i,j
!modelE!real*8    :: apartmolar_tot

!
! define soa species
!
issoa(1)=n_isopp1a
issoa(2)=n_isopp2a

!
! create whichsoa from issoa, in order to correlate the two variables
!
whichsoa=0
do i=1,nsoa
  do j=1,ntm
    if (issoa(i)==j) then
      whichsoa(j)=i
    endif
  enddo
enddo

!
! mass based stoicheiometric coefficients
!
apartmass=0.d0
!modelE!! a-pinene SOAb formation, Presto et al., EST, 2005
!modelE!#ifdef SOA_FULL
!modelE!apartmass(whichsoa(iapinp1a_nox))=0.0138
!modelE!apartmass(whichsoa(iapinp2a_nox))=0.461
!modelE!#endif
!modelE!#ifdef SOA_LUMPED_NOX_DEP
!modelE!apartmass_nox(whichsoa(iapinp1a_hox))=0.0138
!modelE!apartmass_nox(whichsoa(iapinp2a_hox))=0.461
!modelE!#endif
!modelE!apartmass(whichsoa(iapinp1a_hox))=0.192
!modelE!apartmass(whichsoa(iapinp2a_hox))=0.215
!modelE!! xylene SOAa formation, Song et al., EST, 2005 (first two lines) or Ng et al., ACP, 2007 (latter two lines)
!modelE!#ifndef SOA_MINIMUM
!modelE!#  ifdef SOA_FULL
!modelE!apartmass(whichsoa(ixylp1a_nox))=0.049
!modelE!apartmass(whichsoa(ixylp2a_nox))=0.178
!modelE!!apartmass(whichsoa(ixylp1a_nox))=0.031
!modelE!!apartmass(whichsoa(ixylp2a_nox))=0.090
!modelE!#  endif
!modelE!#  ifdef SOA_LUMPED_NOX_DEP
!modelE!apartmass_nox(whichsoa(ixylp1a_hox))=0.049
!modelE!apartmass_nox(whichsoa(ixylp2a_hox))=0.178
!modelE!!apartmass_nox(whichsoa(ixylp1a_hox))=0.031
!modelE!!apartmass_nox(whichsoa(ixylp2a_hox))=0.090
!modelE!#  endif
!modelE!apartmass(whichsoa(ixylp1a_hox))=0.024
!modelE!apartmass(whichsoa(ixylp2a_hox))=0.152
!modelE!!apartmass(whichsoa(ixylp1a_hox))=0.30
!modelE!!apartmass(whichsoa(ixylp2a_hox))=0.
!modelE!#endif
! isoprene + OH low NOx SOAb formation, Henze and Seinfeld, GRL 2006; scaled, based on a-pinene
!modelE!#ifdef SOA_FULL
apartmass(whichsoa(n_isopp1a))=0.232d0!/0.125*apartmass(whichsoa(iapinp1a_nox))
apartmass(whichsoa(n_isopp2a))=0.0288d0!/0.102*apartmass(whichsoa(iapinp2a_nox))
!modelE!apartmass(whichsoa(iisopp1a_nox))=0.232/0.125*apartmass(whichsoa(iapinp1a_nox))
!modelE!apartmass(whichsoa(iisopp2a_nox))=0.0288/0.102*apartmass(whichsoa(iapinp2a_nox))
!modelE!#endif
!modelE!#ifdef SOA_LUMPED_NOX_DEP
!modelE!apartmass_nox(whichsoa(iisopp1a_hox))=0.232/0.125*apartmass_nox(whichsoa(iapinp1a_hox))
!modelE!apartmass_nox(whichsoa(iisopp2a_hox))=0.0288/0.102*apartmass_nox(whichsoa(iapinp2a_hox))
!modelE!#endif
!modelE!apartmass(whichsoa(iisopp1a_hox))=0.232/0.125*apartmass(whichsoa(iapinp1a_hox))
!modelE!apartmass(whichsoa(iisopp2a_hox))=0.0288/0.102*apartmass(whichsoa(iapinp2a_hox))
!modelE!! b-pinene SOAb formation, Griffin et al., JGR, 1999; scaled, based on a-pinene
!modelE!#ifndef SOA_MINIMUM
!modelE!#  ifdef SOA_FULL
!modelE!apartmass(whichsoa(ibpinp1a_nox))=0.026/0.125*apartmass(whichsoa(iapinp1a_nox))
!modelE!apartmass(whichsoa(ibpinp2a_nox))=0.485/0.102*apartmass(whichsoa(iapinp2a_nox))
!modelE!#  endif
!modelE!#  ifdef SOA_LUMPED_NOX_DEP
!modelE!apartmass_nox(whichsoa(ibpinp1a_hox))=0.026/0.125*apartmass_nox(whichsoa(iapinp1a_hox))
!modelE!apartmass_nox(whichsoa(ibpinp2a_hox))=0.485/0.102*apartmass_nox(whichsoa(iapinp2a_hox))
!modelE!#  endif
!modelE!apartmass(whichsoa(ibpinp1a_hox))=0.026/0.125*apartmass(whichsoa(iapinp1a_hox))
!modelE!apartmass(whichsoa(ibpinp2a_hox))=0.485/0.102*apartmass(whichsoa(iapinp2a_hox))
!modelE!#endif
!modelE!! toluene SOAa formation, Odum et al., Science, 1997; scaled, based on xylene (first two lines) or Ng et al., ACP, 2007, no scaling (latter two lines)
!modelE!#ifdef SOA_FULL
!modelE!apartmass(whichsoa(itolp1a_nox))=0.071/0.038*apartmass(whichsoa(ixylp1a_nox))
!modelE!apartmass(whichsoa(itolp2a_nox))=0.138/0.167*apartmass(whichsoa(ixylp2a_nox))
!modelE!!apartmass(whichsoa(itolp1a_nox))=0.058
!modelE!!apartmass(whichsoa(itolp2a_nox))=0.113
!modelE!#endif
!modelE!#ifdef SOA_LUMPED_NOX_DEP
!modelE!apartmass_nox(whichsoa(itolp1a_hox))=0.071/0.038*0.049!apartmass_nox(whichsoa(ixylp1a_hox)) ! HARDCODED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!modelE!apartmass_nox(whichsoa(itolp2a_hox))=0.138/0.167*0.178!apartmass_nox(whichsoa(ixylp2a_hox)) ! HARDCODED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!modelE!!apartmass_nox(whichsoa(itolp1a_hox))=0.058
!modelE!!apartmass_nox(whichsoa(itolp2a_hox))=0.113
!modelE!#endif
!modelE!apartmass(whichsoa(itolp1a_hox))=0.071/0.038*0.024!apartmass(whichsoa(ixylp1a_hox)) ! HARDCODED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!modelE!apartmass(whichsoa(itolp2a_hox))=0.138/0.167*0.152!apartmass(whichsoa(ixylp2a_hox)) ! HARDCODED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!modelE!!apartmass(whichsoa(itolp1a_hox))=0.36
!modelE!!apartmass(whichsoa(itolp2a_hox))=0.

!
! molar based stoicheiometric coefficients
!
apartmolar=0.d0
!modelE!#ifdef SOA_FULL
apartmolar(whichsoa(n_isopp1a))=apartmass(whichsoa(n_isopp1a))*tr_mm(n_isoprene)/tr_mm(n_isopp1a)
apartmolar(whichsoa(n_isopp2a))=apartmass(whichsoa(n_isopp2a))*tr_mm(n_isoprene)/tr_mm(n_isopp2a)
!modelE!apartmolar(whichsoa(iisopp1a_nox))=apartmass(whichsoa(iisopp1a_nox))*tr_mm(iisop)/tr_mm(iisopp1a_nox)
!modelE!apartmolar(whichsoa(iisopp2a_nox))=apartmass(whichsoa(iisopp2a_nox))*tr_mm(iisop)/tr_mm(iisopp2a_nox)
!modelE!#endif
!modelE!#ifdef SOA_LUMPED_NOX_DEP
!modelE!apartmolar_nox(whichsoa(iisopp1a_hox))=apartmass_nox(whichsoa(iisopp1a_hox))*tr_mm(iisop)/tr_mm(iisopp1a_hox)
!modelE!apartmolar_nox(whichsoa(iisopp2a_hox))=apartmass_nox(whichsoa(iisopp2a_hox))*tr_mm(iisop)/tr_mm(iisopp2a_hox)
!modelE!#endif
!modelE!apartmolar(whichsoa(iisopp1a_hox))=apartmass(whichsoa(iisopp1a_hox))*tr_mm(iisop)/tr_mm(iisopp1a_hox)
!modelE!apartmolar(whichsoa(iisopp2a_hox))=apartmass(whichsoa(iisopp2a_hox))*tr_mm(iisop)/tr_mm(iisopp2a_hox)
!modelE!#ifdef SOA_FULL
!modelE!apartmolar(whichsoa(iapinp1a_nox))=apartmass(whichsoa(iapinp1a_nox))*tr_mm(iapin)/tr_mm(iapinp1a_nox)
!modelE!apartmolar(whichsoa(iapinp2a_nox))=apartmass(whichsoa(iapinp2a_nox))*tr_mm(iapin)/tr_mm(iapinp2a_nox)
!modelE!#endif
!modelE!#ifdef SOA_LUMPED_NOX_DEP
!modelE!apartmolar_nox(whichsoa(iapinp1a_hox))=apartmass_nox(whichsoa(iapinp1a_hox))*tr_mm(iapin)/tr_mm(iapinp1a_hox)
!modelE!apartmolar_nox(whichsoa(iapinp2a_hox))=apartmass_nox(whichsoa(iapinp2a_hox))*tr_mm(iapin)/tr_mm(iapinp2a_hox)
!modelE!#endif
!modelE!apartmolar(whichsoa(iapinp1a_hox))=apartmass(whichsoa(iapinp1a_hox))*tr_mm(iapin)/tr_mm(iapinp1a_hox)
!modelE!apartmolar(whichsoa(iapinp2a_hox))=apartmass(whichsoa(iapinp2a_hox))*tr_mm(iapin)/tr_mm(iapinp2a_hox)
!modelE!#ifndef SOA_MINIMUM
!modelE!#  ifdef SOA_FULL
!modelE!apartmolar(whichsoa(ibpinp1a_nox))=apartmass(whichsoa(ibpinp1a_nox))*tr_mm(ibpin)/tr_mm(ibpinp1a_nox)
!modelE!apartmolar(whichsoa(ibpinp2a_nox))=apartmass(whichsoa(ibpinp2a_nox))*tr_mm(ibpin)/tr_mm(ibpinp2a_nox)
!modelE!apartmolar_tot=apartmolar(whichsoa(ibpinp1a_nox))+apartmolar(whichsoa(ibpinp2a_nox))
!modelE!if (apartmolar_tot > 1.) then
!modelE!  apartmolar(whichsoa(ibpinp1a_nox)) = apartmolar(whichsoa(ibpinp1a_nox))/apartmolar_tot
!modelE!  apartmolar(whichsoa(ibpinp2a_nox)) = apartmolar(whichsoa(ibpinp2a_nox))/apartmolar_tot
!modelE!  apartmass(whichsoa(ibpinp1a_nox)) = apartmolar(whichsoa(ibpinp1a_nox))*tr_mm(ibpinp1a_nox)/tr_mm(ibpin)
!modelE!  apartmass(whichsoa(ibpinp2a_nox)) = apartmolar(whichsoa(ibpinp2a_nox))*tr_mm(ibpinp2a_nox)/tr_mm(ibpin)
!modelE!endif
!modelE!#  endif
!modelE!#  ifdef SOA_LUMPED_NOX_DEP
!modelE!apartmolar_nox(whichsoa(ibpinp1a_hox))=apartmass_nox(whichsoa(ibpinp1a_hox))*tr_mm(ibpin)/tr_mm(ibpinp1a_hox)
!modelE!apartmolar_nox(whichsoa(ibpinp2a_hox))=apartmass_nox(whichsoa(ibpinp2a_hox))*tr_mm(ibpin)/tr_mm(ibpinp2a_hox)
!modelE!apartmolar_tot=apartmolar_nox(whichsoa(ibpinp1a_hox))+apartmolar_nox(whichsoa(ibpinp2a_hox))
!modelE!if (apartmolar_tot > 1.) then
!modelE!  apartmolar_nox(whichsoa(ibpinp1a_hox)) = apartmolar_nox(whichsoa(ibpinp1a_hox))/apartmolar_tot
!modelE!  apartmolar_nox(whichsoa(ibpinp2a_hox)) = apartmolar_nox(whichsoa(ibpinp2a_hox))/apartmolar_tot
!modelE!  apartmass_nox(whichsoa(ibpinp1a_hox)) = apartmolar_nox(whichsoa(ibpinp1a_hox))*tr_mm(ibpinp1a_hox)/tr_mm(ibpin)
!modelE!  apartmass_nox(whichsoa(ibpinp2a_hox)) = apartmolar_nox(whichsoa(ibpinp2a_hox))*tr_mm(ibpinp2a_hox)/tr_mm(ibpin)
!modelE!endif
!modelE!#  endif
!modelE!apartmolar(whichsoa(ibpinp1a_hox))=apartmass(whichsoa(ibpinp1a_hox))*tr_mm(ibpin)/tr_mm(ibpinp1a_hox)
!modelE!apartmolar(whichsoa(ibpinp2a_hox))=apartmass(whichsoa(ibpinp2a_hox))*tr_mm(ibpin)/tr_mm(ibpinp2a_hox)
!modelE!apartmolar_tot=apartmolar(whichsoa(ibpinp1a_hox))+apartmolar(whichsoa(ibpinp2a_hox))
!modelE!if (apartmolar_tot > 1.) then
!modelE!  apartmolar(whichsoa(ibpinp1a_hox)) = apartmolar(whichsoa(ibpinp1a_hox))/apartmolar_tot
!modelE!  apartmolar(whichsoa(ibpinp2a_hox)) = apartmolar(whichsoa(ibpinp2a_hox))/apartmolar_tot
!modelE!  apartmass(whichsoa(ibpinp1a_hox)) = apartmolar(whichsoa(ibpinp1a_hox))*tr_mm(ibpinp1a_hox)/tr_mm(ibpin)
!modelE!  apartmass(whichsoa(ibpinp2a_hox)) = apartmolar(whichsoa(ibpinp2a_hox))*tr_mm(ibpinp2a_hox)/tr_mm(ibpin)
!modelE!endif
!modelE!#endif
!modelE!#ifdef SOA_FULL
!modelE!apartmolar(whichsoa(itolp1a_nox))=apartmass(whichsoa(itolp1a_nox))*tr_mm(itol)/tr_mm(itolp1a_nox)
!modelE!apartmolar(whichsoa(itolp2a_nox))=apartmass(whichsoa(itolp2a_nox))*tr_mm(itol)/tr_mm(itolp2a_nox)
!modelE!#endif
!modelE!#ifdef SOA_LUMPED_NOX_DEP
!modelE!apartmolar_nox(whichsoa(itolp1a_hox))=apartmass_nox(whichsoa(itolp1a_hox))*tr_mm(itol)/tr_mm(itolp1a_hox)
!modelE!apartmolar_nox(whichsoa(itolp2a_hox))=apartmass_nox(whichsoa(itolp2a_hox))*tr_mm(itol)/tr_mm(itolp2a_hox)
!modelE!#endif
!modelE!apartmolar(whichsoa(itolp1a_hox))=apartmass(whichsoa(itolp1a_hox))*tr_mm(itol)/tr_mm(itolp1a_hox)
!modelE!apartmolar(whichsoa(itolp2a_hox))=apartmass(whichsoa(itolp2a_hox))*tr_mm(itol)/tr_mm(itolp2a_hox)
!modelE!#ifndef SOA_MINIMUM
!modelE!#  ifdef SOA_FULL
!modelE!apartmolar(whichsoa(ixylp1a_nox))=apartmass(whichsoa(ixylp1a_nox))*tr_mm(ixyl)/tr_mm(ixylp1a_nox)
!modelE!apartmolar(whichsoa(ixylp2a_nox))=apartmass(whichsoa(ixylp2a_nox))*tr_mm(ixyl)/tr_mm(ixylp2a_nox)
!modelE!#  endif
!modelE!#  ifdef SOA_LUMPED_NOX_DEP
!modelE!apartmolar_nox(whichsoa(ixylp1a_hox))=apartmass_nox(whichsoa(ixylp1a_hox))*tr_mm(ixyl)/tr_mm(ixylp1a_hox)
!modelE!apartmolar_nox(whichsoa(ixylp2a_hox))=apartmass_nox(whichsoa(ixylp2a_hox))*tr_mm(ixyl)/tr_mm(ixylp2a_hox)
!modelE!#  endif
!modelE!apartmolar(whichsoa(ixylp1a_hox))=apartmass(whichsoa(ixylp1a_hox))*tr_mm(ixyl)/tr_mm(ixylp1a_hox)
!modelE!apartmolar(whichsoa(ixylp2a_hox))=apartmass(whichsoa(ixylp2a_hox))*tr_mm(ixyl)/tr_mm(ixylp2a_hox)
!modelE!#endif

!modelE!agas=0.d0
!modelE!agas(iisop)=1.-(&
!modelE!#ifdef SOA_FULL
!modelE!                0.5*(apartmolar(whichsoa(iisopp1a_nox))+apartmolar(whichsoa(iisopp2a_nox)))+0.5*&
!modelE!#endif
!modelE!#ifdef SOA_LUMPED_NOX_DEP
!modelE!                0.5*(apartmolar_nox(whichsoa(iisopp1a_hox))+apartmolar_nox(whichsoa(iisopp2a_hox)))+0.5*&
!modelE!#endif
!modelE!                (apartmolar(whichsoa(iisopp1a_hox))+apartmolar(whichsoa(iisopp2a_hox)))&
!modelE!               )
!modelE!agas(iapin)=1.-(&
!modelE!#ifdef SOA_FULL
!modelE!                0.5*(apartmolar(whichsoa(iapinp1a_nox))+apartmolar(whichsoa(iapinp2a_nox)))+0.5*&
!modelE!#endif
!modelE!#ifdef SOA_LUMPED_NOX_DEP
!modelE!                0.5*(apartmolar_nox(whichsoa(iapinp1a_hox))+apartmolar_nox(whichsoa(iapinp2a_hox)))+0.5*&
!modelE!#endif
!modelE!                (apartmolar(whichsoa(iapinp1a_hox))+apartmolar(whichsoa(iapinp2a_hox)))&
!modelE!               )
!modelE!#ifndef SOA_MINIMUM
!modelE!agas(ibpin)=1.-(&
!modelE!#  ifdef SOA_FULL
!modelE!                0.5*(apartmolar(whichsoa(ibpinp1a_nox))+apartmolar(whichsoa(ibpinp2a_nox)))+0.5*&
!modelE!#  endif
!modelE!#  ifdef SOA_LUMPED_NOX_DEP
!modelE!                0.5*(apartmolar_nox(whichsoa(ibpinp1a_hox))+apartmolar_nox(whichsoa(ibpinp2a_hox)))+0.5*&
!modelE!#  endif
!modelE!                (apartmolar(whichsoa(ibpinp1a_hox))+apartmolar(whichsoa(ibpinp2a_hox)))&
!modelE!               )
!modelE!#endif
!modelE!agas(itol)=1.-(&
!modelE!#ifdef SOA_FULL
!modelE!               0.5*(apartmolar(whichsoa(itolp1a_nox))+apartmolar(whichsoa(itolp2a_nox)))+0.5*&
!modelE!#endif
!modelE!#ifdef SOA_LUMPED_NOX_DEP
!modelE!               0.5*(apartmolar_nox(whichsoa(itolp1a_hox))+apartmolar_nox(whichsoa(itolp2a_hox)))+0.5*&
!modelE!#endif
!modelE!               (apartmolar(whichsoa(itolp1a_hox))+apartmolar(whichsoa(itolp2a_hox)))&
!modelE!              )
!modelE!#ifndef SOA_MINIMUM
!modelE!agas(ixyl)=1.-(&
!modelE!#  ifdef SOA_FULL
!modelE!               0.5*(apartmolar(whichsoa(ixylp1a_nox))+apartmolar(whichsoa(ixylp2a_nox)))+0.5*&
!modelE!#  endif
!modelE!#  ifdef SOA_LUMPED_NOX_DEP
!modelE!               0.5*(apartmolar_nox(whichsoa(ixylp1a_hox))+apartmolar_nox(whichsoa(ixylp2a_hox)))+0.5*&
!modelE!#  endif
!modelE!               (apartmolar(whichsoa(ixylp1a_hox))+apartmolar(whichsoa(ixylp2a_hox)))&
!modelE!               )
!modelE!#endif

do i=1,soacomp
  do j=1,soacomp
    if (Lambda(i,j).eq.-1.0d0) Lambda(i,j)=Lambda(j,i) ! symmetric interactions
  enddo
enddo

! Correct tr_mm in order to use the real MW of species in meanmw and activity coefficient calculations only
mw=tr_mm
mw(n_bcii)=170.d0
mw(n_bcia)=170.d0
mw(n_bcb)=170.d0
mw(n_ocii)=170.d0
mw(n_ocia)=170.d0
mw(n_ocb)=170.d0

do i=1,ntm
  molec2ug(i)=mw(i)*1.d12/avog
enddo

out_line="Initialization of SOA formation completed"
call write_parallel(trim(out_line),crit=.true.)

end subroutine soa_init


subroutine soa_aerosolphase(L,changeL,bypfactor)

use TRACER_COM, only: n_bcii,n_bcia,n_bcb,n_ocii,n_ocia,n_ocb,&
#ifdef TRACERS_NITRATE
                      n_nh4,n_no3p,&
#endif
                      n_msa,n_so4,&
                      mass2vol
use TRCHEM_Shindell_COM, only: y
implicit none

!@var y0_ug concentration of species (in molecules/cm3) before chemistry and partitioning
!@var y_ug concentration of species (in molecules/cm3) after chemistry and partitioning
!@var y_mw =y/tr_mm*mw
real*8,dimension(ntm)                    :: y0_ug,y_ug,y_mw
real*8, intent(in)                       :: bypfactor
real*8, dimension(LM,ntm), intent(inout) :: changeL
!@var jl looping index
integer                                  :: jl
integer, intent(in)                      :: L
!
! For correction of Kp due to composition
!
!@var AEROtot the total aerosol concentration, which includes both the semivolatile aerosol phase and
!@+           the non-volatile aerosol phase that affects condensation, BEFORE new condensation
!@var sigma-ij,sigma_ik,sigma_jk used for the activity coefficient calculations
!@var zc the activity coefficient BEFORE new condensation and BEFORE the validity check (needed only for diagnistic output)
!@var meanmw the mean molecular weight of the total condensation affected aerosol phase, BEFORE new condensation
real*8                      :: AEROtot,sigma_ij,sigma_ik,sigma_kj,zc,meanmw
!@var imf,jmf,kmf indices for the activity coefficient calculation
integer                     :: imf,jmf,kmf
!@var xmf the mole fraction of species in the aerosol phase BEFORE new condensation
!@var zcoef final value of activity coefficient BEFORE new condensation (this is the value that is
!@+         going to be used). The allowed range is 0.3 < zcoef < 5.0
real*8, dimension(soacomp)  :: xmf,zcoef
!@var kp final value of partitioning coefficient after applying all relevant parameters to the reference case kpart_ref
real*8,dimension(n_soa_i:n_soa_e) :: kp !!!!!!!!!!!!!!!!!!!!!! CONTINUE FROM HERE !!!!!!!!!!!!!!!!!!!!!!!!! the documentation
! Other SOA parameters
real*8, parameter           :: M0err=1.d-10
real*8                      :: M0,PCP,M0temp,M0a,M0b,M0err_curr
real*8, dimension(nsoa)     :: soamass,partfact
integer                     :: i,iternum
!modelE!character(len=11), parameter :: itermethod="chord" ! valid values are "bisectional" and "chord"
real*8                      :: x1,x2,y1,y2,a,b

logical, parameter :: SO4part=.true.
#ifdef TRACERS_NITRATE
logical, parameter :: NH4part=.true.
#endif
logical, parameter :: SOAevap=.true.

!DO JL=1,LM
DO JL=L,L

!
! Change concentrations to ug/m3
! Note that we want the concentrations AFTER the chemistry, 
! thus we apply changeL artificially (but not hardcoded)
!
  do i=1,ntm
    y0_ug(i)=(y(i,jl)+changeL(jl,i)*bypfactor*mass2vol(i))*molec2ug(i)
  enddo
  y_ug=y0_ug

!
! Calculate the primary aerosol able to absorb semivolatile compounds
! Partitioning in both OC and BC always happens
! Partitioning in (SO4+MSA) and/or NH4/NO3 is also an option
!
  PCP=y_ug(n_bcii)+y_ug(n_bcia)+y_ug(n_bcb)+&
      y_ug(n_ocii)+y_ug(n_ocia)+y_ug(n_ocb)
  if(SO4part) PCP=PCP+y_ug(n_msa)+y_ug(n_so4)
#ifdef TRACERS_NITRATE
  if(NH4part) PCP=PCP+y_ug(n_nh4)+y_ug(n_no3p)
#endif
!
! Correct the Kp to take into account the change of activity coefficient due to change in composition
! using the Wilson equation. Method described in Bowman and Karamalegos, EST, 2002, 36, 2701-2707.
!
  do i=1,ntm
    y_mw(i)=y(i,jl)/tr_mm(i)*mw(i)
  enddo
  AEROtot=y_mw(n_bcii)+y_mw(n_bcia)+y_mw(n_bcb)+&
          y_mw(n_ocii)+y_mw(n_ocia)+y_mw(n_ocb)
  if(SO4part) AEROtot=AEROtot+y_mw(n_msa)+y_mw(n_so4)
#ifdef TRACERS_NITRATE
  if(NH4part) AEROtot=AEROtot+y_mw(n_nh4)+y_mw(n_no3p)
#endif
  do i=1,nsoa
    AEROtot=AEROtot+y_mw(issoa(i))
  enddo
  if (AEROtot.le.0.d0) goto 60
!
! Calculate mole fraction of individual species.
!
  xmf=0.d0
  xmf(imfisop)=(&
!modelE!#ifdef SOA_FULL
                y_mw(n_isopp1a)+y_mw(n_isopp2a)&
!modelE!#endif
!modelE!               +y(jl,iisopp1a_hox)/tr_mm(iisopp1a_hox)*mw(iisopp1a_hox)+y(jl,iisopp2a_hox)/tr_mm(iisopp2a_hox)*mw(iisopp2a_hox)&
               )/AEROtot
!modelE!  xmf(imfter)=(&
!modelE!#ifdef SOA_FULL
!modelE!                y(jl,iapinp1a_nox)/tr_mm(iapinp1a_nox)*mw(iapinp1a_nox)+y(jl,iapinp2a_nox)/tr_mm(iapinp2a_nox)*mw(iapinp2a_nox)&
!modelE!#endif
!modelE!              +y(jl,iapinp1a_hox)/tr_mm(iapinp1a_hox)*mw(iapinp1a_hox)+y(jl,iapinp2a_hox)/tr_mm(iapinp2a_hox)*mw(iapinp2a_hox)&
!modelE!#ifndef SOA_MINIMUM
!modelE!#  ifdef SOA_FULL
!modelE!              +y(jl,ibpinp1a_nox)/tr_mm(ibpinp1a_nox)*mw(ibpinp1a_nox)+y(jl,ibpinp2a_nox)/tr_mm(ibpinp2a_nox)*mw(ibpinp2a_nox)&
!modelE!#  endif
!modelE!              +y(jl,ibpinp1a_hox)/tr_mm(ibpinp1a_hox)*mw(ibpinp1a_hox)+y(jl,ibpinp2a_hox)/tr_mm(ibpinp2a_hox)*mw(ibpinp2a_hox)&
!modelE!#endif
!modelE!              )/AEROtot
!modelE!  xmf(imfaro)=(&
!modelE!#ifdef SOA_FULL
!modelE!               y(jl,itolp1a_nox)/tr_mm(itolp1a_nox)*mw(itolp1a_nox)+y(jl,itolp2a_nox)/tr_mm(itolp2a_nox)*mw(itolp2a_nox)&
!modelE!#endif
!modelE!              +y(jl,itolp1a_hox)/tr_mm(itolp1a_hox)*mw(itolp1a_hox)+y(jl,itolp2a_hox)/tr_mm(itolp2a_hox)*mw(itolp2a_hox)&
!modelE!#ifndef SOA_MINIMUM
!modelE!#  ifdef SOA_FULL
!modelE!              +y(jl,ixylp1a_nox)/tr_mm(ixylp1a_nox)*mw(ixylp1a_nox)+y(jl,ixylp2a_nox)/tr_mm(ixylp2a_nox)*mw(ixylp2a_nox)&
!modelE!#  endif
!modelE!              +y(jl,ixylp1a_hox)/tr_mm(ixylp1a_hox)*mw(ixylp1a_hox)+y(jl,ixylp2a_hox)/tr_mm(ixylp2a_hox)*mw(ixylp2a_hox)&
!modelE!#endif
!modelE!              )/AEROtot
  xmf(imfbcii)=y_mw(n_bcii)/AEROtot
  xmf(imfbcia)=y_mw(n_bcia)/AEROtot
  xmf(imfbcb)=y_mw(n_bcb)/AEROtot
  xmf(imfocii)=y_mw(n_ocii)/AEROtot
  xmf(imfocia)=y_mw(n_ocia)/AEROtot
  xmf(imfocb)=y_mw(n_ocb)/AEROtot
!
! If partitioning does not occur on some aerosol species, it's mole fraction equals to zero, no matter the real
! concentration is, in order not to affect the calculation of the activity coefficient.
!
  if(SO4part) then
    xmf(imfinorg)=xmf(imfinorg)+y_mw(n_msa)/AEROtot
    xmf(imfinorg)=xmf(imfinorg)+y_mw(n_so4)/AEROtot
  endif
#ifdef TRACERS_NITRATE
  if(NH4part) then
    xmf(imfinorg)=xmf(imfinorg)+y_mw(n_nh4)/AEROtot
    xmf(imfinorg)=xmf(imfinorg)+y_mw(n_no3p)/AEROtot
  endif
#endif
!
! Calculate activity coefficient zcoef
!
  zcoef=1.d0
!  do kmf=1,soacomp
  do kmf=imfter,imfisop
    sigma_kj=0.d0
    do jmf=1,soacomp
      sigma_kj=sigma_kj+xmf(jmf)*Lambda(kmf,jmf)
    enddo
    sigma_ik=0.d0
    do imf=1,soacomp
      sigma_ij=0.d0
      do jmf=1,soacomp
        sigma_ij=sigma_ij+xmf(jmf)*Lambda(imf,jmf)
      enddo
      sigma_ik=sigma_ik+xmf(imf)*Lambda(imf,kmf)/sigma_ij
    enddo
    zc=exp(1.d0-log(sigma_kj)-sigma_ik)
    zcoef(kmf)=max(0.3d0,min(5.0d0,zc))
    if((zc.lt.0.3d0).or.(zc.gt.5.0d0)) then
      write(out_line,*)'WARNING: zcoef set to',zcoef(kmf),' (was',zc,')'
      call write_parallel(trim(out_line),crit=.true.)
    endif
  enddo
!
! Calculate mean molecular weight
!
  meanmw=y_mw(n_ocii)*mw(n_ocii)/AEROtot
  meanmw=meanmw+y_mw(n_ocia)*mw(n_ocia)/AEROtot
  meanmw=meanmw+y_mw(n_ocb)*mw(n_ocb)/AEROtot
  meanmw=meanmw+y_mw(n_bcii)*mw(n_bcii)/AEROtot
  meanmw=meanmw+y_mw(n_bcia)*mw(n_bcia)/AEROtot
  meanmw=meanmw+y_mw(n_bcb)*mw(n_bcb)/AEROtot
  if(SO4part) then
    meanmw=meanmw+y_mw(n_msa)*mw(n_msa)/AEROtot
    meanmw=meanmw+y_mw(n_so4)*mw(n_so4)/AEROtot
  endif
#ifdef TRACERS_NITRATE
  if(NH4part) then
    meanmw=meanmw+y_mw(n_nh4)*mw(n_nh4)/AEROtot
    meanmw=meanmw+y_mw(n_no3p)*mw(n_no3p)/AEROtot
  endif
#endif
  do i=1,nsoa
    meanmw=meanmw+y_mw(issoa(i))*mw(issoa(i))/AEROtot
  enddo
!
! Calculate final value of partitioning coefficient
!
!modelE!#ifdef SOA_FULL
  kp(n_isopp1a)=kpart(jl,whichsoa(n_isopp1a))/zcoef(imfisop)*mw(n_isopp1a)/meanmw
  kp(n_isopp2a)=kpart(jl,whichsoa(n_isopp2a))/zcoef(imfisop)*mw(n_isopp2a)/meanmw
!modelE!#endif
!modelE!  kp(iisopp1a_hox)=kpart(jl,whichsoa(iisopp1a_hox))/zcoef(imfisop)*mw(iisopp1a_hox)/meanmw
!modelE!  kp(iisopp2a_hox)=kpart(jl,whichsoa(iisopp2a_hox))/zcoef(imfisop)*mw(iisopp2a_hox)/meanmw
!modelE!#ifdef SOA_FULL
!modelE!  kp(iapinp1a_nox)=kpart(jl,whichsoa(iapinp1a_nox))/zcoef(imfter)*mw(iapinp1a_nox)/meanmw
!modelE!  kp(iapinp2a_nox)=kpart(jl,whichsoa(iapinp2a_nox))/zcoef(imfter)*mw(iapinp2a_nox)/meanmw
!modelE!#endif
!modelE!  kp(iapinp1a_hox)=kpart(jl,whichsoa(iapinp1a_hox))/zcoef(imfter)*mw(iapinp1a_hox)/meanmw
!modelE!  kp(iapinp2a_hox)=kpart(jl,whichsoa(iapinp2a_hox))/zcoef(imfter)*mw(iapinp2a_hox)/meanmw
!modelE!#ifndef SOA_MINIMUM
!modelE!#  ifdef SOA_FULL
!modelE!  kp(ibpinp1a_nox)=kpart(jl,whichsoa(ibpinp1a_nox))/zcoef(imfter)*mw(ibpinp1a_nox)/meanmw
!modelE!  kp(ibpinp2a_nox)=kpart(jl,whichsoa(ibpinp2a_nox))/zcoef(imfter)*mw(ibpinp2a_nox)/meanmw
!modelE!#  endif
!modelE!  kp(ibpinp1a_hox)=kpart(jl,whichsoa(ibpinp1a_hox))/zcoef(imfter)*mw(ibpinp1a_hox)/meanmw
!modelE!  kp(ibpinp2a_hox)=kpart(jl,whichsoa(ibpinp2a_hox))/zcoef(imfter)*mw(ibpinp2a_hox)/meanmw
!modelE!#endif
!modelE!#ifdef SOA_FULL
!modelE!  kp(itolp1a_nox)=kpart(jl,whichsoa(itolp1a_nox))/zcoef(imfaro)*mw(itolp1a_nox)/meanmw
!modelE!  kp(itolp2a_nox)=kpart(jl,whichsoa(itolp2a_nox))/zcoef(imfaro)*mw(itolp2a_nox)/meanmw
!modelE!#endif
!modelE!  kp(itolp1a_hox)=kpart(jl,whichsoa(itolp1a_hox))/zcoef(imfaro)*mw(itolp1a_hox)/meanmw
!modelE!  kp(itolp2a_hox)=kpart(jl,whichsoa(itolp2a_hox))/zcoef(imfaro)*mw(itolp2a_hox)/meanmw
!modelE!#ifndef SOA_MINIMUM
!modelE!#  ifdef SOA_FULL
!modelE!  kp(ixylp1a_nox)=kpart(jl,whichsoa(ixylp1a_nox))/zcoef(imfaro)*mw(ixylp1a_nox)/meanmw
!modelE!  kp(ixylp2a_nox)=kpart(jl,whichsoa(ixylp2a_nox))/zcoef(imfaro)*mw(ixylp2a_nox)/meanmw
!modelE!#  endif
!modelE!  kp(ixylp1a_hox)=kpart(jl,whichsoa(ixylp1a_hox))/zcoef(imfaro)*mw(ixylp1a_hox)/meanmw
!modelE!  kp(ixylp2a_hox)=kpart(jl,whichsoa(ixylp2a_hox))/zcoef(imfaro)*mw(ixylp2a_hox)/meanmw
!modelE!#endif
!
! Add previous gas phase (soagas - after destruction), previous particulate phase and soachem in variable soamass (ug m-3)
! If no evaporation, add only soagas and soachem
!
  do i=1,nsoa
!    soamass(i)=soagas(jl,i)+soachem(jl,i)
    soamass(i)=y_ug(issoa(i)-1)
    if (SOAevap) soamass(i)=soamass(i)+y_ug(issoa(i))
  enddo
  if (sum(soamass)==0.d0) goto 60
  M0a=PCP
  M0b=PCP
  do i=1,nsoa
    M0b=M0b+soamass(i)
    if(.not.SOAevap) then
      M0a=M0a+y_ug(issoa(i))
      M0b=M0b+y_ug(issoa(i))
    endif
  enddo
  iternum=0
  M0err_curr=max(tiny(M0b),M0b*M0err)!min(M0err,max(1.d-15,(M0b-M0a)*M0err))
!
! Iteration - chord method
!
  M0=M0b
11  continue ! Try with another M0
  if (iternum==0) then
    x2=M0
    x1=x2-M0err_curr
  else
    x1=M0
    x2=x1+M0err_curr
  endif
  y1=PCP-x1
  y2=PCP-x2
  if(.not.SOAevap) then
    do i=1,nsoa
      y1=y1+y_ug(issoa(i))
      y2=y2+y_ug(issoa(i))
    enddo
  endif
  do i=1,nsoa
    partfact(i)=kp(issoa(i))*x2/(1.d0+kp(issoa(i))*x2)
    y2=y2+soamass(i)*partfact(i)
    partfact(i)=kp(issoa(i))*x1/(1.d0+kp(issoa(i))*x1)
    y1=y1+soamass(i)*partfact(i)
  enddo
  iternum=iternum+1
!  if(abs(y1).lt.M0err_curr/2.d0)goto 20
  if(abs(y1).lt.M0err)goto 20
!
! Stop iteration after 100 iterations (no solution found)
!
  if(iternum.ge.100) goto 30
!
! If -M0err/2<M0temp<M0err/2, M0 is the solution (goto 20)
! Otherwise change the limits of M0a and M0b and restart (goto 11)
!
  a=(y2-y1)/(x2-x1)
  b=y1-a*x1
  M0=-b/a
  if (((iternum==1).and.(M0==x2)).or.((iternum>1).and.(M0==x1))) goto 20
  if (M0.lt.0.d0) then
    write(out_line,*)'WARNING: M0 set to zero (was ',M0,' after ',iternum,' iterations)'
    call write_parallel(trim(out_line),crit=.true.)
    M0=0.d0
  endif
  goto 11
!
! No solution found, aerosol phase set to previous values, chemistry still applies to gas-phase species
!
30 continue
  write(*,"('WARNING: Too many iteration steps. SOA forced to the previous values. M0=',e,' PCP=',e)") M0,PCP
  goto 60
!
! Found solution for M0
!
20 continue
  if (M0.ne.0.d0) then
    do i=1,nsoa
      partfact(i)=max(0.d0,partfact(i))
      if(SOAevap) then
        y_ug(issoa(i))=soamass(i)*partfact(i)             ! Calculate new aerosol-phase concentration
        y_ug(issoa(i)-1)=y_ug(issoa(i))/(kp(issoa(i))*M0) ! Calculate new gas-phase concentration
      else
        y_ug(issoa(i))=y_ug(issoa(i))+soamass(i)*partfact(i)
        y_ug(issoa(i)-1)=max(0.d0,y_ug(issoa(i)-1)-soamass(i)*partfact(i))
      endif
    enddo
  endif

60 continue ! everytime aerosols are zero, goto 60 is called
!
! Save results to changeL (kg)
!
  do i=1,nsoa
    changeL(jl,issoa(i)-1)=changeL(jl,issoa(i)-1)+(y_ug(issoa(i)-1)-y0_ug(issoa(i)-1))/&
                           (molec2ug(issoa(i)-1)*bypfactor*mass2vol(issoa(i)-1))
    changeL(jl,issoa(i))=changeL(jl,issoa(i))+(y_ug(issoa(i))-y0_ug(issoa(i)))/&
                         (molec2ug(issoa(i))*bypfactor*mass2vol(issoa(i)))
  enddo

ENDDO !JL

end subroutine soa_aerosolphase


real*8 function KpCALC(dH,Ksc,temp,Tsc)
!-------------------------------------------------------------------------------
!KpCALC calculates the temperature dependence of the partitioning coefficients
!      interface
!      ---------
!      call KpCALC(dH,Ksc,temp)
!        where: dH is the enthalpy of vaporization of the selected species, in KJ/mol
!               Ksc is the partitioning coefficient of the selected species at temperature Tsc
!               temp is the temperature in current box
!      reference
!      ---------
!      Chung and Seinfeld, JGR, 2002
!-------------------------------------------------------------------------------
use CONSTANT, only: gasc
implicit none
real*8, intent(in):: dH         ! KJ/mol
real*8, intent(in):: Ksc,temp,Tsc

KpCALC=Ksc*(temp/Tsc)*exp(dH*1.d3/gasc*(1.d0/temp-1.d0/Tsc))

end function KpCALC


end module TRACERS_SOA
