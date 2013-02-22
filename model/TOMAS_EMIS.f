#include "rundeck_opts.h"

!@sum  TOMAS_EMIS: Emission size assumptioned used for TOMAS model 
      MODULE TOMAS_EMIS


      USE TRACER_COM, only : nbins
      IMPLICIT NONE 

!@param scalesizeSO4 : SO4 emission size distribution (mass fraction at each bin)
!@param scaleCARBO30/scaleCARBO100 : EC/OC emission size distribution 
!@+   CARBO30 is for fossil fuel and uses bimodal distribution [Ban-Weiss et al, 2010] 
!@+          (NMD: 17.5 nm and GSD:1.6 ; NMD: 60 nm and GSD: 1.9)  
!@+   CARBO100 is for Biofuel and Biomass burning emissions
!@_          (NMD of 100 nm and GSD of 2)
!@param scalesizeClay :Clay emissions assuming a lognormal with NMD=0.14 um and Sigma=2
!@+      sum of scalesizeClay = ~1 (~5% of total clay emission will be in Dp>2um) 
!@param scalesizeSilt :Silt (Dp>2um) emissions assuming a lognormal with NMD=1.14 um and Sigma=2
!@+      sum of scalesizeSilt = 0.8415 (~15% of total silt emission will be missing 
!@+     due to upper size limit, ~10um, in TOMAS. ~8% will be in clay size range)

!@param scalesizeSalt :Sea-salt emissions based on Gong et al. (2004) 
!@+      Unlike other scalesize, sum of scalesizeSalt do not added up to 1. 
!@+      These numbers are from Kostas. 
!@+      Seasalt emission will be computed in TRACERS_AEROSOL_KOCH_e4.f 

#if (defined TOMAS_12_10NM) 
      real*8, parameter :: scalesizeSO4(nbins)=(/
     &     4.3760E-02, 6.2140E-02, 3.6990E-02, 1.8270E-02,
     &     4.2720E-02, 1.1251E-01, 1.9552E-01, 2.2060E-01,
     &     1.6158E-01, 7.6810E-02, 2.8884E-02, 2.0027E-04/)

      real*8, parameter :: scalesizeCARBO30(nbins)=(/
     &     1.1291E-03,4.9302E-03,1.2714E-02,3.6431E-02,
     &     1.0846E-01,2.1994E-01,2.7402E-01,2.0750E-01,
     &     9.5304E-02,2.6504E-02,1.2925E-02,1.6069E-05/) 
      real*8, parameter :: scalesizeCARBO100(nbins)=(/
     &     1.9827E-06,3.9249E-05,5.0202E-04,4.1538E-03,
     &     2.2253E-02,7.7269E-02,1.7402E-01,2.5432E-01,
     &     2.4126E-01,1.4856E-01,7.6641E-02,9.8120E-04/) 

      real*8, parameter :: scalesizeClay(nbins)=(/
     *    3.883E-08,1.246E-06,2.591E-05,3.493E-04,
     *    3.059E-03,1.741E-02,6.444E-02,1.553E-01,
     *    2.439E-01,2.495E-01,2.530E-01,1.292E-02/)
      real*8, parameter :: scalesizeSilt(nbins)=(/
     *    2.310E-17,5.376E-15,8.075E-13,7.832E-11,
     *    4.910E-09,1.991E-07,5.229E-06,8.900E-05,
     *    9.831E-04,7.054E-03,2.183E-01,6.150E-01/)

      real*8, parameter :: scalesizeSalt(nbins)=(/
     *     3.7028916E-23,5.7123891E-22,1.1839277E-20,
     *     2.8514155E-19,4.5906335E-18,4.3888240E-17,
     *     1.8982601E-16,4.6295965E-16,7.4091490E-16,
     *     1.4321188E-15,3.6773131E-14,7.6566702E-14/)


#elif (defined TOMAS_12_3NM) 
      real*8, parameter :: scalesizeSO4(nbins)=(/0.,0.,0.,
     &     4.3760E-02, 6.2140E-02, 3.6990E-02, 1.8270E-02,
     &     4.2720E-02, 1.1251E-01, 1.9552E-01, 2.2060E-01,
     &     1.6158E-01, 7.6810E-02, 2.8884E-02, 2.0027E-04/)

      real*8, parameter :: scalesizeCARBO30(nbins)=(/0.,0.,0.,
     &     1.1291E-03,4.9302E-03,1.2714E-02,3.6431E-02,
     &     1.0846E-01,2.1994E-01,2.7402E-01,2.0750E-01,
     &     9.5304E-02,2.6504E-02,1.2925E-02,1.6069E-05/) 
      real*8, parameter :: scalesizeCARBO100(nbins)=(/0.,0.,0.,
     &     1.9827E-06,3.9249E-05,5.0202E-04,4.1538E-03,
     &     2.2253E-02,7.7269E-02,1.7402E-01,2.5432E-01,
     &     2.4126E-01,1.4856E-01,7.6641E-02,9.8120E-04/) 

      real*8, parameter :: scalesizeClay(nbins)=(/0.,0.,0.,
     *    3.883E-08,1.246E-06,2.591E-05,3.493E-04,
     *    3.059E-03,1.741E-02,6.444E-02,1.553E-01,
     *    2.439E-01,2.495E-01,2.530E-01,1.292E-02/)
      real*8, parameter :: scalesizeSilt(nbins)=(/0.,0.,0.,
     *    2.310E-17,5.376E-15,8.075E-13,7.832E-11,
     *    4.910E-09,1.991E-07,5.229E-06,8.900E-05,
     *    9.831E-04,7.054E-03,2.183E-01,6.150E-01/)

      real*8, parameter :: scalesizeSalt(nbins)=(/0.,0.,0.,
     *     3.7028916E-23,5.7123891E-22,1.1839277E-20,
     *     2.8514155E-19,4.5906335E-18,4.3888240E-17,
     *     1.8982601E-16,4.6295965E-16,7.4091490E-16,
     *     1.4321188E-15,3.6773131E-14,7.6566702E-14/)

#endif

#if (defined TOMAS_15_10NM)
      real*8, parameter :: scalesizeSO4(nbins)=(/
     &     4.3760E-02, 6.2140E-02, 3.6990E-02, 1.8270E-02,
     &     4.2720E-02, 1.1251E-01, 1.9552E-01, 2.2060E-01,
     &     1.6158E-01, 7.6810E-02, 2.8884E-02, 2.0027E-04/)

      real*8, parameter :: scalesizeCARBO30(nbins)=(/
     &     3.8100E-03,2.0700E-02,7.2900E-02,1.6750E-01,
     &     2.5000E-01,2.4500E-01,1.5510E-01,6.4200E-02,1.7320E-02,
     &     3.0320E-03,3.7079E-04,2.3265E-07/) ! use for fossil fuel

      real*8, parameter :: scalesizeCARBO100(nbins)=(/
     &     1.9827E-06,3.9249E-05,5.0202E-04,4.1538E-03,
     &     2.2253E-02,7.7269E-02,1.7402E-01,2.5432E-01,2.4126E-01,
     &     1.4856E-01,7.6641E-02,9.8120E-04/) ! use for biomass burning
#elif (defined TOMAS_15_3NM)
      real*8, parameter :: scalesizeSO4(nbins)=(/
     &     0.,0.,0.,
     &     4.3760E-02, 6.2140E-02, 3.6990E-02, 1.8270E-02,
     &     4.2720E-02, 1.1251E-01, 1.9552E-01, 2.2060E-01,
     &     1.6158E-01, 7.6810E-02, 2.8884E-02, 2.0027E-04/)

      real*8, parameter :: scalesizeCARBO30(nbins)=(/
     &     0.,0.,0.,
     &     3.8100E-03,2.0700E-02,7.2900E-02,1.6750E-01,
     &     2.5000E-01,2.4500E-01,1.5510E-01,6.4200E-02,1.7320E-02,
     &     3.0320E-03,3.7079E-04,2.3265E-07/) 

      real*8, parameter :: scalesizeCARBO100(nbins)=(/
     &     0.,0.,0.,
     &     1.9827E-06,3.9249E-05,5.0202E-04,4.1538E-03,
     &     2.2253E-02,7.7269E-02,1.7402E-01,2.5432E-01,2.4126E-01,
     &     1.4856E-01,7.6641E-02,9.8120E-04/) 
#endif

!#if (defined TOMAS_30_10NM)

!#elif (defined TOMAS_30_3NM)
           
!#endif

      END MODULE TOMAS_EMIS
