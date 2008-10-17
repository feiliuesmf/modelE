#include "rundeck_opts.h"

      MODULE TRCHEM_Shindell_COM
!@sum  TRCHEM_Shindell_COM declares variables for tracer chemistry
!@+    and sources.
!@auth Drew Shindell (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on various chemistry modules of B436Tds3YM23 model)
c
#ifdef TRACERS_ON
      USE MODEL_COM, only  : im,jm,lm,psf,ptop,sig,sige,dsig,bydsig,
     &                       dtsrc,Itime,ItimeI,T,JEQ
      USE CONSTANT, only   : pi, mair, mwat, radian,avog
      USE DYNAMICS, only   : am, byam, PMID, PK
      USE GEOM, only       : BYDXYP,dxyp
      USE RAD_COM, only    : rcloudfj=>rcld !!! ,salbfj=>salb
      USE TRACER_COM, only : ntm, trm, TR_MM, ntm_soa

      IMPLICIT NONE
      SAVE

C**************  P  A  R  A  M  E  T  E  R  S  *******************
!@param p_1 number of reactants per reaction
!@param p_2 number of rxns in assembled lists (check with print rxn list)
!@param p_3 number of rxns in assembled lists (check with print rxn list)
!@param p_4 number of rxns in assembled lists (check with print rxn list)
!@param p_5 number of levels from top down with SRB flux
!@param n_rx maximum number of chemical reactions
!@param n_bi maximum number of bimolecular reactions
!@param n_tri maximum number of trimolecular reactions
!@param n_bnd1 maximum number of spectral bands 1
!@param n_bnd2 maximum number of spectral bands 2
!@param n_bnd3 maximum number of spectral bands 3
!@param n_nst maximum number of monomolecular decompositions
!@param n_fam maximum number of chemical families
!@param n_oig max number of optically-important gases
!@param n_srb max number of schumann-runge bands
!@param N__ Number of levels in Mie grid: 2*(2*lpar+2+jaddto(1))+3
!@param M__ Number of Gauss points used
!@param NLFASTJ maximum number levels after inserting extra Mie levels
!@param NS maximum number of species which require J-values calculating
!@param NWFASTJ maximum number of wavelength bins that can be used
!@param JPPJ number of photolysis reactions
!@param JPNL number of photolysis levels
!@param NJVAL Number of species for which to calculate J-values
!@param szamax max Zenith Angle(98 deg at 63 km;99 degrees at 80 km)
!@param dtaumax max optical depth above which must instert new level
!@param ZZHT Scale height above top of atmosphere (cm)
!@param odmax Maximum allowed optical depth, above which they're scaled
!@param luselb Use reflective photolysis boundary treatment
!@param zlbatm Optical depth above which to set lower boundary
!@param NFASTJ number of quadrature points in OPMIE
!@param MFIT expansion of phase function in OPMIE
!@param MFASTJ ?
!@param CMEQ1 ?
!@param nc total number of molecules included (incl. O2 and N2)
!@param ny number of chemically calculated gases (no O2 or N2)
!@param numfam number of chemical families
!@param n_phot how often to do photolysis (in increments of DTsrc)
!@param O3MULT =2.14d-2 This is the conversion from (atm*cm) units
!@+     (i.e. 1000 Dobson Units) to KG/m2. It is: 
!@+     1.E4*2.69E19*48./6.02E26 where 1.E4 is cm2/m2, 2.69E19 is 
!@+     molecules/cm3 at 1 atm pressure, 48. is molecular wt of O3,
!@+     and 6.02E26 is Avogadro's number in molecules/Kmol.
!@param cpd conversion from molecules/cm3 to mole/m
!@param BYO3MULT = 1/O3MULT
!@param pfix_O2 fixed ratio of O2/M
!@param pfix_H2 fixed ratio of H2/M
!@param pfix_Aldehyde fixed ratio of Aldehyde/M for initial conditions
!@param checktracer_on integer to turn on the checktracer call
!@param correct_strat_Ox logical tourn on lower strat Ox corrections
!@param MWabyMWw ratio of molecular weights of air/water
!@param O3_1_fact factor to alter surface O3 that is passed to FASTJ
!@+     this is fastj level 1, not model level 1.  Currently, it is 
!@+     decreased by a factor of (972/1000)mb
!@param RGAMMASULF unknown
!@param RKBYPIM=8.*RBOLTZ/pi/MASSN2O55=8.*1.38062D-23/3.14159/1.793D-25
!@param cboltz Boltzman's Constant = 1.3806d-19
!@param dlogp 10.d0**(-2./16.)
!@param dlogp2 10.d0**(-1./16.)
!@param byradian 1/radian = conversion from radians to degrees
!@param LCOalt number of levels in the COaltIN array 
!@param LCH4alt number of levels in the CH4altIN array
!@param LcorrOx number of levels in the corrOxIN array
!@param JCOlat number of latitudes in the COlat array 
!@param PCOalt pressures at LCOalt levels
!@param PCH4alt pressures at LCH4alt levels
!@param PcorrOx pressures at LcorrOx levels
!@param NCFASTJ2 number of levels in the fastj2 atmosphere
!@param NBFASTJ for fastj2 (=LM+1)
!@param MXFASTJ "Number of aerosol/cloud types supplied from CTM"
!@param dtausub # optic. depths at top of cloud requiring subdivision
!@param dsubdiv additional levels in first dtausub of cloud (fastj2) 
!@param masfac Conversion factor, pressure to column density (fastj2)
!@param NP maximum aerosol phase functions
!@param MIEDX2 choice of aerosol types for fastj2
!@param T_thresh threshold temperature used in master chem
!@param n2o_pppv default N2O L=1 overwriting in pppv
!@param cfc_pppv default CFC L=1 overwriting in pppv
!@param cfc_rad95 the average L=1 radiation code CFC11+CFC12 value
!@+     for 1995 (pppv). 
!@param fact_cfc ratio of our default CFC L=1 overwriting to the 
!@+     radiation's 1995 L=1 CFC11+CFC12 value.
      INTEGER, PARAMETER ::
     & LCOalt =   23,
     & JCOlat =   19,
     & LCH4alt=    6,
     & LcorrOX=    7,
     & p_1   =     2 
#ifdef SHINDELL_STRAT_CHEM
      INTEGER, PARAMETER ::
     & p_2   =   209,
     & p_3   =   500,
     & p_4   =   209,
     & n_rx  =   110,
     & n_bi  =    94,
     & n_tri =    11,
     & n_nst =     3,
     & nc     =   53+ntm_soa,     !formerly in param sub
     & ny     =   51+ntm_soa,     !formerly in param sub  
     & numfam =    4,     !formerly in param sub  
     & nC2O3=     26+ntm_soa,
     & nXO2=      27+ntm_soa,
     & nXO2N=     28+ntm_soa,
     & nRXPAR=    29+ntm_soa,
     & nROR=      30+ntm_soa,
     & nAldehyde= 31+ntm_soa,
     & nH2O=      32+ntm_soa,
     & nCH3O2=    33+ntm_soa,
     & nH2=       34+ntm_soa,
     & nOH=       35+ntm_soa,
     & nHO2=      36+ntm_soa,
     & nO3=       37+ntm_soa,
     & nO=        38+ntm_soa,
     & nO1D=      39+ntm_soa,
     & nNO=       40+ntm_soa,
     & nNO2=      41+ntm_soa,
     & nNO3=      42+ntm_soa,
     & nHONO=     43+ntm_soa,
     & nCl2O2=    44+ntm_soa,
     & nClO=      45+ntm_soa,
     & nOClO=     46+ntm_soa,
     & nCl2=      47+ntm_soa,
     & nCl=       48+ntm_soa,
     & nBrCl=     49+ntm_soa,
     & nBrO=      50+ntm_soa,
     & nBr=       51+ntm_soa,
     & nO2=       52+ntm_soa,
     & nM=        53+ntm_soa,     !you must always put nM last (highest number)
     & JPPJ   =   28,
     & NJVAL  =   27,     !formerly read in from jv_spec00_15.dat
     & NLFASTJ=  400,     !400 is arbitrary for now
     & NWFASTJ=   18, 
     & JPNL   =   LM,      ! OK? used to be set to 23
     & NCFASTJ2 = 2*LM+2,  ! fastj2
     & NBFASTJ  = LM+1,    ! fastj2
     & MXFASTJ  =  3,      ! fastj2
     & n_fam =     5,      ! fastj2
     & NP       = 21       ! fastj2
      INTEGER, PARAMETER, DIMENSION(MXFASTJ) :: MIEDX2 = (/3,10,14/)
!     (/black carbon absorber, water cloud (Deirmenjian 8 micron), 
!      irregular ice cloud (Mishchenko)/)
#else
      INTEGER, PARAMETER ::
     & p_2   =   111,
     & p_3   =   200,
     & p_4   =    70,
     & n_rx  =   101,
     & n_bi  =    96,
     & n_tri =    14,
     & n_nst =     5,
     & nc     =   35+ntm_soa,     !formerly in param sub
     & ny     =   33+ntm_soa,     !formerly in param sub  
     & numfam =    2,     !formerly in param sub  
     & nC2O3=     16+ntm_soa,
     & nXO2=      17+ntm_soa,
     & nXO2N=     18+ntm_soa,
     & nRXPAR=    19+ntm_soa,
     & nROR=      20+ntm_soa,
     & nAldehyde= 21+ntm_soa,
     & nH2O=      22+ntm_soa,
     & nCH3O2=    23+ntm_soa,
     & nH2=       24+ntm_soa,
     & nOH=       25+ntm_soa,
     & nHO2=      26+ntm_soa,
     & nO3=       27+ntm_soa,
     & nO=        28+ntm_soa,
     & nO1D=      29+ntm_soa,
     & nNO=       30+ntm_soa,
     & nNO2=      31+ntm_soa,
     & nNO3=      32+ntm_soa,
     & nHONO=     33+ntm_soa,
     & nO2=       34+ntm_soa,
     & nM=        35+ntm_soa,     !you must always put nM last (highest number)
     & JPPJ   =   16,
     & NJVAL  =   16,     !formerly read in from jv_spec00_15.dat
     & NLFASTJ=  350,     !350 is arbitrary for now
     & NWFASTJ=   15, 
     & JPNL   =   12,     ! change to LS1???
     & n_fam =     4,
     & NP     =    9
#endif
      INTEGER, PARAMETER ::
     & p_5   =    14,
     & n_bnd1=    31,
     & n_bnd2=    87,
     & n_bnd3=   107,
     & n_oig =     3,
     & n_srb =    18,
C ----------------------------------------------     
c     & n_Ox=        1,    ! note, these
c     & n_NOx=       2,    ! first 15 species are
c     & n_N2O5=      3,    ! tracers, and therefore
c     & n_HNO3=      4,    ! these parameters are
c     & n_H2O2=      5,    ! to be defined in 
c     & n_CH3OOH=    6,    ! TRACER_COM.f.
c     & n_HCHO=      7,    ! Note the UNDERSCORE!
c     & n_HO2NO2=    8,    !  T
c     & n_CO=        9,    !  R
c     & n_CH4=      10,    !  A
c     & n_PAN=      11,    !  C
c     & n_Isoprene= 12,    !  E
c     & n_AlkylNit= 13,    !  R
c     & n_Alkenes=  14,    !  S
c     & n_Paraffin= 15,    ! ---------------
C ----------------------------------------------   
     & N__=     1800,     !jan00, was 450, then 900 in Nov99
     & M__=        4,
     & NS     =   51,
     & MFIT   =    2*M__,
     & NFASTJ =    4,
     & MFASTJ =    1,
     & n_phot=     2  
      INTEGER, PARAMETER, DIMENSION(12) :: MDOFM =
     & (/31,59,90,120,151,181,212,243,273,304,334,365/)
     
      REAL*8, PARAMETER ::  O3MULT       = 2.14d-2,
     &                      BYO3MULT     = 1./O3MULT,
     &                      T_thresh     = 203.d0,
     &                      pfix_O2      = 0.209476d0,
     &                      pfix_H2      = 560.d-9,
     &                      pfix_Aldehyde= 2.d-9,
     &                      MWabyMWw     = mair/mwat,
     &                      O3_1_fact    = 0.972d0,
     &                      RGAMMASULF   = 0.1d0,
     &                      RKBYPIM      = 1.961d2,
     &                      cboltz       = 1.3806d-19,
     &                      dlogp        = 7.49894209d-1, !=10^(-.125)
     &                      szamax       = 98.0d0,
     &                      dtaumax      = 1.0d0,
     &                      ZZHT         = 5.d5,
     &                      odmax        = 200.d0,
     &                      zlbatm       = 4.d0,
     &                      CMEQ1        = 0.25d0,
     &                      byradian     = 1.d0/radian,
     &                      cpd          = 1.d6/avog
#ifdef SHINDELL_STRAT_CHEM
     &                     ,cfc_pppv     = 1722.d-12
     &                     ,n2o_pppv     = 290.d-9
     &                     ,cfc_rad95    = 794.d-12 
     &                     ,fact_cfc     = cfc_pppv/cfc_rad95
     &                     ,dtausub      = 1.d0
     &                     ,dsubdiv      = 1.d1
     &                     ,dlogp2       = 8.65964323d-1 !=10^(-.0625)
     &                     ,masfac=100.d0*6.022d23/28.97d0/9.8d0/10.d0
#endif
C Please note: since PCOalt is essentially the nominal 
C pressures for the 23-level GCM, I'm going to use it
C to define BrOx, ClOx,ClONOs,HCL,OxIC,CFCIC,N2OICX,CH4ICX too:
      REAL*8, PARAMETER, DIMENSION(LCOalt) :: PCOalt = (/
     & 0.9720D+03,0.9445D+03,0.9065D+03,
     & 0.8515D+03,0.7645D+03,0.6400D+03,0.4975D+03,0.3695D+03,
     & 0.2795D+03,0.2185D+03,0.1710D+03,0.1335D+03,0.1016D+03,
     & 0.7120D+02,0.4390D+02,0.2470D+02,0.1390D+02,0.7315D+01,
     & 0.3045D+01,0.9605D+00,0.3030D+00,0.8810D-01,0.1663D-01/)
      REAL*8, PARAMETER, DIMENSION(LcorrOx) :: PcorrOx =
     & (/278.7d0, 217.9d0, 170.5d0, 133.1d0, 101.3d0, 71.0d0, 43.8d0/)
      REAL*8, PARAMETER, DIMENSION(M__)  :: EMU = (/.06943184420297D0,
     &        .33000947820757D0,.66999052179243D0,.93056815579703D0/), 
     &                                    WTFASTJ=(/.17392742256873D0,
     &         .32607257743127D0,.32607257743127D0,.17392742256873D0/)
      REAL*8, PARAMETER, DIMENSION(LCOalt) ::  
     &     COaltIN = (/2d0,1.5625d0,1.375d0,1.25d0,1.125d0,1.0625d0,
     &     1d0,1d0,1d0,1d0,1d0,.5d0,.375d0,.2d0,.2d0,.2d0,.2d0,.2d0,
     &     .25d0,.4d0,2.5d0,12d0,60d0/)
#ifdef SHINDELL_STRAT_CHEM
     &     ,BrOxaltIN = (/1.d-2,1.d-2,1.d-2,1.d-2,1.d-2,1.d-2,1.d-2,
     &     1.d-2,1.d-2,1.d-2,1.d-2,0.12d0,0.12d0,0.12d0,0.12d0,0.06d0,
     &     0.06d0,0.06d0,0.06d0,0.06d0,0.06d0,0.06d0,0.06d0/)
     &     ,ClOxaltIN = (/1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,
     &     1.d0,1.d0,8.d0,8.d0,8.d0,8.d0,8.d1,8.d1,8.d1,8.d1,8.d0,8.d0,
     &     8.d0,8.d0/)
     &     ,ClONO2altIN = (/1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,
     &     1.d0,1.d0,1.d0,1.d0,1.d0,5.d1,5.d1,5.d1,5.d1,5.d1,5.d1,5.d1,
     &     5.d1,5.d1,5.d1/)
     &     ,HClaltIN = (/1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,
     &     1.d0,1.d0,2.5d1,4.0d1,9.0d1,1.7d2,1.9d2,2.5d2,2.5d2,2.5d2,
     &     2.5d2,2.5d2,2.5d2,2.5d2/)
#endif    
      REAL*8, PARAMETER, DIMENSION(LCH4alt) :: PCH4alt = 
     &                     (/569d0, 150d0, 100d0, 32d0, 3.2d0, 0.23d0/)
      REAL*8, PARAMETER, DIMENSION(LCH4alt) ::   
     &   CH4altINT =(/1.79d0, 1.75d0, 1.620d0,1.460d0,0.812d0,0.230d0/),
     &   CH4altINX =(/1.79d0, 1.75d0, 1.440d0,1.130d0,0.473d0,0.202d0/)    
      REAL*8, PARAMETER, DIMENSION(JCOlat)  :: COlat = (/40.,40.,40.,40.
     *     ,45.,50.,60.,70.,80.,90.,110.,125.,140.,165.,175.,180.,170.
     *     ,165.,150./)
     
!@dbparam ch4_init_sh,ch4_init_nh initial methane conc. (ppmv) 
!@+       defaults are for 1990
!@dbparam fix_CH4_chemistry (1=YES 0=NO) whether or not to used a fixed
!@+       value for methane in the chemistry code. USE -1 for initial
!@+       conditions from file CH4_IC (but L>LS1-1 only now!)
!@dbparam scale_ch4_IC_file multiplicative factor of CH4 IC if 
!@+       fix_CH4_chemistry=-1 (but only above LS1-1 !)
!@dbparam pfix_CH4_N fixed ratio of CH4/M in North. Hemis. (if used)
!@dbparam pfix_CH4_S fixed ratio of CH4/M in South. Hemis. (if used)
!@dbparam use_rad_ch4 =1 replaces CH4 surface sources with L=1
!+        overwriting with radiation code values.
!@dbparam use_rad_n2o =1 as ch4 case above
!@dbparam use_rad_cfc =1 as ch4 case above
!@dbparam which_trop 1=ls1-1 is tropopause, 0=LTROPO(I,J) is tropopause
!@dbparam PI_run used to turn on (1) and off (0) use of PI_ratio*
!@dbparam PIratio_N preindustrial ratio for NOx, HNO3, N2O5, HO2NO2
!@+       initial conditions and stratospheric overwriting.
!@dbparam PIratio_CO_T preindustrial ratio for tropospheric CO
!@dbparam PIratio_CO_S preindustrial ratio for stratospheric CO
!@dbparam PIratio_other PI ratio: PAN,Isoprene,AlkyNit,Alkenes,Paraffin
!@dbparam PIratio_indus preindustrial ratio for industrial sources
!@dbparam PIratio_bburn preindustrial ratio for biomass burning sources
!@dbparam PIratio_N2O preindustrial ratio for N2O ICs and L=1 overwrite
!@dbparam PIratio_CFC preindustrial ratio for CFC ICs and L=1 overwrite
!@dbparam rad_FL whether(>0) or not(=0) to have fastj photon flux vary 
!@+       with model time (JYEAR, JMON, JDAY) 

      INTEGER ::        fix_CH4_chemistry = 0
     &                 ,which_trop        = 0
     &                 ,PI_run            = 0
     &                 ,rad_FL            = 0
     &                 ,use_rad_ch4       = 0
     &                 ,use_rad_n2o       = 0
     &                 ,use_rad_cfc       = 0
     &                 ,checktracer_on    = 0
      REAL*8 ::             ch4_init_sh   = 1.750d0,
     &                      ch4_init_nh   = 1.855d0,
     &                      pfix_CH4_S    = 1.75d-6,
     &                      pfix_CH4_N    = 1.855d-6,
     &                      scale_ch4_IC_file= 1.d0, 
     &                      PIratio_N     = 0.667d0,
     &                      PIratio_CO_T  = 0.667d0,
     &                      PIratio_CO_S  = 0.500d0,
     &                      PIratio_other = 0.500d0,
     &                      PIratio_indus = 0.000d0,
     &                      PIratio_bburn = 0.100d0
#ifdef SHINDELL_STRAT_CHEM
     &                     ,PIratio_N2O   = 0.896d0
     &                     ,PIratio_CFC   = 0.000d0
#endif

      LOGICAL, PARAMETER :: luselb            = .false.,
     &                      correct_strat_Ox  = .true.

C**************  V  A  R  I  A  B  L  E  S *******************  
!@var nn reactant's number in mol list, first index reactant 1 or 2,
!@+      second - reaction number
!@var nnr reaction product's number in mol list, indicies like nn
!@var kss mollst number for product gases from photolysis
!@var ks mollst number for source gas in photolysis reaction
!@var nps reaction numbers by molecule, photolytic production
!@var nds reaction numbers by molecule, photolytic destruction
!@var npnr reaction numbers by molecule, photolytic production
!@var ndnr reaction numbers by molecule, photolytic destruction
!@var kps reaction numbers by molecule, chemical production
!@var kds reaction numbers by molecule, chemical destruction
!@var kpnr reaction numbers by molecule, chemical production
!@var kdnr reaction numbers by molecule, chemical destruction
!@var fam ___?
!@var nst reverse reaction number for dissociation reactions
!@var lbeg beginning of spectral range for photodissociation for ind.
!@+   gas (18 Sch-runge bands, then 200-730 nm in 5 nm steps)
!@var nir length of photodissociation spectra of gas (# of non-zero
!@+   absorption cross sections)
!@var lprn,jprn,iprn l, j, and i point for chemistry debugging
!@var ay name of gas being considered
!@var y concentration of gas, 1st index=gas number, 2nd=verticle level
!@var rr rate constant of chemical reaction, first index - reaction
!@+   number, 2nd is verticle level
!@var ss photodissociation coefficient, indicies; rxn #,L,I,J
!@var pe rate constant for bimolecular chemical reaction
!@var ea activation energy constant for bimolecular chemical reactions
!@var ro,r1,sn,sb rate parameters for trimolecular reactions
!@var sigg effective absorption cross section for radiation, first
!@+   index spectral interval number, second photolysis reaction number
!@var conc concentration of optically important gases (O2 & O3), first
!@+   vertivle level, second=gas number (1=O2,2=O3)
!@var qfu flux of solar radiation in the upper atmosphere in 18
!@+   schumann-runge bands 175 - 200 nm
!@var qf flux of solar radiation, first index for spectral interval
!@+   number, second - verticle level (photons/cm^2*c)
!@var wlt wavelength from 200 to 730 nm in 5 nm steps (107 intervals)
!@var sO3 absorption cross section of ozone (cm^2)
!@var sO2 absorption cross section of oxygen (cm^2)
!@var sech cross section of optically important gases, first index
!@+   gas number (1=O2,2=O3), second - spectral interval number 
!@var TXL temperature profile
!@var prnrts logical: print rate of each chemical reaction?
!@var prnchg logical: print chemical changes?
!@var prnls logical: print reaction lists by species?
!@var yNO3,pHOx,pNOx,pOx,yCH3O2,yC2O3,yROR,yXO2,yAldehyde,yXO2N,yRXPAR?
!@var yCl2,yCl2O2 3D arrays to remember some non-tracer species...
!@var NCFASTJ number of levels in the fastj atmosphere
!@var title_aer_pf titles read from aerosol phase function file
!@var TITLE0 blank title read in I think
!@var TITLEJ titles read from O2, O3, and other species X-sections
!@var jlabel Reference label identifying appropriate J-value to use
!@var jind mapping index for jvalues? 
!@var jndlev Levels at which we want J-values (centre of CTM levels)
!@var jaddlv Additional levels associated with each level
!@var jaddto Cumulative total of new levels to be added
!@var NW1, NW2 beginning, ending wavelength for wavelength "bins"
!@var MIEDX Type of aerosol scattering, currently 6 set up:
!@+   1=Rayly 2=iso 3=iso-equiv 4=bkgrd-sulf,5=volc-sulf,6=liq water
!@var NAA Number of categories for scattering phase functions
!@var npdep Number of pressure dependencies
!@var jpdep Index of cross sections requiring P dependence
!@var PFASTJ pressure sent to FASTJ
!@var PFASTJ2 pressure at level boundarie, sent to FASTJ2
!@var nss this is a copy of JPPJ that is read in from a file
!@var jfacta Quantum yield (or multiplication factor) for photolysis
!@var WBIN Boundaries of wavelength bins
!@var WL Centres of wavelength bins - 'effective wavelength'
!@var NWWW Number of wavelength bins, from NW1:NW2
!@var FL Solar flux incident on top of atmosphere (cm-2.s-1)
!@var   Rayleigh parameters (effective cross-section) (cm2)
!@var DUMMY placeholder for reading FL if rad_FL>0
!@var FLX temp array for varying FL if rad_FL>0   
!@var QBC Black Carbon abs. extinct. (specific cross-sect.m2/g)
!@var QO2      O2 cross-sections
!@var QO3      O3 cross-sections
!@var Q1D      O3 => O(1D) quantum yield
!@var TQQ      Temperature for supplied cross sections
!@var QQQ      Supplied cross sections in each wavelength bin (cm2)
!@var QAAFASTJ Aerosol scattering phase functions
!@var NK       Number of wavelengths at which functions are supplied
!@var WAAFASTJ Wavelengths for the NK supplied phase functions
!@var PAA Scaling for extinctions
!@var zpdep    Pressure dependencies by wavelength bin
!@var lpdep    Label for pressure dependence
!@var OREF     O3 reference profile
!@var TREF     temperature reference profile
!@var BREF     black carbon reference profile
!@var OREF2    fastj2 O3 reference profile
!@var TREF2    fastj2 temperature reference profile
!@var BREF2    fastj2 black carbon reference profile
!@var U0 cosine of the solar zenith angle
!@var ZFASTJ Altitude of each fastj pressure level (approx.) (cm)
!@var ZFASTJ2 Altitude of boundaries of model levels (cm) fastj2
!@var RFLECT Surface albedo (Lamertian) in fastj
!@var odtmp Optical depth (temporary array)
!@var odsum Column optical depth
!@var nlbatm Level of lower photolysis boundary - usually surface ('1')
!@var aer fastj aerosol profile?
!@var O3J Ozone profile on photolysis grid
!@var TJ Temperature profile on photolysis grid
!@var TJ2 Temperature profile on fastj2 photolysis grid
!@var DBC Mass of Black Carbon at each pressure level (g.cm-3)
!@var DBC2 fastj2 Mass of Black Carbon at each model level (g/cm-3)
!@var FFF Actinic flux at each level for each wavelength bin and level
!@var DMFASTJ Total number density at each pressure level (cm-3)
!@var DMFASTJ2 fastj2 Air column for each model level (molec/cm2)
!@var DO3 Ozone column number density at each pressure level (molec/cm2)
!@var DO32 fastj2 Ozone number density at each pressure level (")
!@var XQO2   Absorption cross-section of O2
!@var XQO3   Absorption cross-section of O3
!@var DTAUDZ   Local extinction at each point
!@var PIRAY    Contribution of Rayleigh scattering to extinction
!@var PIAER    Contribution of Aerosol scattering to extinction
!@var TTAU     Opt depth of air vert'y above each point(to top of atm)
!@var XLTAU    TTAU along the slant path
!@var FTAU     Attenuation of solar beam
!@var RZ      Distance from centre of Earth to each point (cm)
!@var RQ      Square of distance ratios
!@var TANHT   Tangent height for the current SZA
!@var XL      Slant path between points
!@var WTAU    Weighted slant path - each side of point
!@var dpomega   change in pomega per increment  (linear)
!@var POMEGAJ  Scattering phase function. the 2nd dimension on POMEGAJ
!@+   was 30 in the 9-layer code, which is (2*LM+2)+10. But it seems to
!@+   only need +1 (I.e. NCFASTJ+1)
!@var POMEGA,ZTAU,FZ,ZREFL,jndlv,FJFASTJ,EMU,ZFLUX,ZREFL,ZU0,WFASTJ ?
!@var PM0,PM,BFASTJ,AFASTJ,AAFASTJ,WTFASTJ,CC,HFASTJ,C1,SFASTJ,U1,V1 ?
!@var RR2 former RR from fastj ?
!@var nfam number of beginning molecule of each chemical family
!@var ZJ photodissociation coefficient? (level,reaction)
!@var RCLOUDFJ cloudiness (optical depth) parameter, radiation to fastj
!@var SALBFJ surface albedo parameter from radiation to fastj
!@var COlat carbon monoxide latitude distribution for init cond (ppbv)
!@+ [CO] ppbv based on 10deg lat-variation Badr & Probert 1994 fig 9.
!@var COaltIN adjustments of COlat by altitude (unitless,LCOalt levels)
!@+ Multiplier of free trop [CO] by layer Badr & Probert 94 fig10 & 11,
!@+ Lopez-Valverde et al 93 fig 3, and Warneck 88, ch1 fig14.
!@var OxICIN Ox initial conditions (unit=KG,LCOalt levels)
!@var OxICINL column version of OxICIN
!@var N2OICIN N2O initial conditions (unit=KG,LCOalt levels)
!@var N2OICINL column version of N2OICIN
!@var CH4ICIN CH4 initial conditions (unit=KG,LCOalt levels)
!@var CH4ICINL column version of CH4ICIN
!@var CFCICIN CFC initial conditions (unit=KG,LCOalt levels)
!@var CFCICINL column version of CFCICIN
!@var BrOxaltIN altitude dependence BrOx (unitless,LCOalt levels)
!@var ClOxaltIN altitude dependence ClOx (unitless,LCOalt levels)
!@var ClONO2altIN altitude dependence ClONO2 (unitless,LCOalt levels)
!@var HClaltIN altitude dependence HCl (unitless,LCOalt levels)
!@var CH4altINT tropical strat adjustments to CH4 (LCH4alt levels)
!@var CH4altINX xtra-tropical strat adjustments to CH4 LCH4alt levels)
!@var COalt adjustments of COlat by altitude (unitless,LM levels)
!@var OxIC Ox initial conditions (unit=KG,LM levels)
!@var OxICL column version of OxIC
!@var N2OICX N2O initial conditions (unit=KG,LM levels) X=not Jean's
!@var N2OICL column version of N2OICX
!@var CH4ICX CH4 initial conditions (unit=KG,LM levels) X=not Jean's
!@var CH4ICL column version of CH4ICX
!@var CFCIC CFC initial conditions (unit=KG,LM levels)
!@var CFCICL column version of CFCIC
!@var BrOxalt altitude dependence BrOx (unitless,LM levels)
!@var ClOxalt altitude dependence ClOx (unitless,LM levels)
!@var ClONO2alt altitude dependence ClONO2 (unitless,LM levels)
!@var HClalt altitude dependence HCl (unitless,LM levels)
!@var CH4altT tropical strat adjustments to CH4 (unitless, LM levels)
!@var CH4altX xtra-tropical strat adjustments to CH4 (LM levels)
!@var BYFJM = 1/JM
!@var MODPHOT if MODPHOT=0 do photolysis, else skip it
!@var TX temperature variable for master chem
!@var ta, pres local arrays to hold temperature,pressure
!@var TFASTJ temperature rpofile sent to FASTJ
!@var O3_FASTJ ozone sent to fastj
!@var FASTJLAT,FASTJLON latitude & LONGITUDE (degrees) for use in fastj
!@var SZA the solar zenith angle (degrees)
!@var JFASTJ photolysis rates
!@var sulfate N2O5 sulfate sink (formerly SRC(I,J,L,20) variable)   
!@var dms_offline DMS concentration for HOx sink reactions
!@var so2_offline SO2 concentration for HOx conversion reactions
!@var prod_sulfate  N2O5 change by sulfate reactions in mass units
!@var wprod_sulf N2O5 change by sulfate reactions in molecules/cm3/s
!@var DT2 variable chemical time step, set in masterchem
!@var nr total number of        reactions read in from gs_jpl00_trop_15
!@var nr3 #of trimolecular      reactions read in from gs_jpl00_trop_15
!@var nr2 #of mono+bi-molecular reactions read in from gs_jpl00_trop_15
!@var nmm #of monomolecular     reactions read in from gs_jpl00_trop_15
!@var nhet #of heterogenous     reactions read in from gs_jpl00_trop_15
!@var ratioNs,ratioN2,rNO2frac,rNOfrac,rNOdenom variables for nitrogen
!@+   conservation (strat)
!@var chemrate,photrate ?   
!@var MDOFM cumulative days at end of each month
!@var corrOxIN correction factor to tweak inital Ox in stratosphere
!@+   on LcorrOx levels
!@var corrOx correction factor to tweak inital Ox in stratosphere
!@+   on LM levels
!@var L75P first model level above nominal 75 hPa
!@var L75M first model level below nominal 75 hPa
!@var F75P interpolation coeff. of higher altitude value (units ln(P))
!@var F75M interpolation coeff. of lower altitude value (units ln(P))
!@var L569P first model level above nominal 569 hPa
!@var L569M first model level below nominal 569 hPa
!@var F569P interpolation coeff. of higher altitude value (units ln(P))
!@var F569M interpolation coeff. of lower altitude value (units ln(P))
!@var DU_O3 total column ozone in latitude band
!@var SF3 is H2O photolysis in Schumann-Runge Bands
!@var SF2 is NO photolysis in Schumann-Runge Bands
!@var SF3_fact used to alter SF3 in time (see comments in master)
!@var SF2_fact used to alter SF2 in time (see comments in master)
!@var bin4_1988 fastj2 bin#4 photon flux for year 1988
!@var bin4_1991 fastj2 bin#4 photon flux for year 1991
!@var bin5_1988 fastj2 bin#5 photon flux for year 1988
!@var AER2 fastj2 aerosol profile?
!@var odcol Optical depth at each model level
!@var AMF Air mass factor for slab between level and level above
!@var Jacet photolysis rate for acetone (not done through fastj)
!@var acetone 3D acetone mixing ratio (static for now)
      INTEGER :: nr,nr2,nr3,nmm,nhet,MODPHOT,L75P,L75M,L569P,L569M,
     &lprn,jprn,iprn,NW1,NW2,MIEDX,NAA,npdep,nss,NWWW,NK,nlbatm,NCFASTJ
#ifdef SHINDELL_STRAT_CHEM
      INTEGER, DIMENSION(n_fam)        :: nfam = 
     &     (/37+ntm_soa,40+ntm_soa,44+ntm_soa,50+ntm_soa,0/)
#else
      INTEGER, DIMENSION(n_fam)        :: nfam = 
     &     (/27+ntm_soa,30+ntm_soa,0,0/)
#endif
      INTEGER, DIMENSION(p_1,p_2)      :: nn, nnr, kss
      INTEGER, DIMENSION(p_2)          :: ks
      INTEGER, DIMENSION(p_3)          :: nps, nds, npnr, ndnr
      INTEGER, DIMENSION(p_4)          :: kps, kds, kpnr, kdnr
      INTEGER, DIMENSION(n_nst)        :: nst
      INTEGER, DIMENSION(n_bnd1)       :: lbeg, nir
      INTEGER, DIMENSION(LM)           :: jndlv,jndlev
      INTEGER, DIMENSION(JPPJ)         :: jind
      INTEGER, DIMENSION(NLFASTJ)      :: jaddlv
#ifdef SHINDELL_STRAT_CHEM
      INTEGER, DIMENSION(NLFASTJ)      :: jadsub
      INTEGER, DIMENSION(NLFASTJ+1)    :: jaddto
#else 
      INTEGER, DIMENSION(NLFASTJ)      :: jaddto
#endif
      INTEGER, DIMENSION(NJVAL)        :: jpdep  

C**************  Latitude-Dependant (allocatable) *******************
      REAL*8, ALLOCATABLE, DIMENSION(:)       :: DU_O3
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:)   :: corrOxIN,acetone
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:)   :: corrOx
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: ss
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:)   :: yNO3,pHOx,pNOx,pOx,
     & yCH3O2,yC2O3,yROR,yXO2,yAldehyde,yXO2N,yRXPAR,TX,sulfate,OxIC,
     & CH4ICX,dms_offline,so2_offline,yso2,ydms      
#ifdef SHINDELL_STRAT_CHEM
     & ,pClOx,pClx,pOClOx,pBrOx,yCl2,yCl2O2,N2OICX,CFCIC,SF3,SF2
#endif
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:)   :: OxICIN,CH4ICIN
#ifdef SHINDELL_STRAT_CHEM
     &                                         ,N2OICIN,CFCICIN
#endif

C**************  Not Latitude-Dependant ****************************      
      REAL*8 :: ZFLUX,ZREFL,ZU0,U0,RFLECT,odsum,XLTAU,TANHT,BYFJM,
     & FASTJLAT,FASTJLON,SZA,DT2,F75P,F75M,F569P,F569M
#ifdef SHINDELL_STRAT_CHEM
     & ,ratioNs,ratioN2,rNO2frac,rNOfrac,rNOdenom
     & ,bin4_1991,bin4_1988,bin5_1988
#endif
      REAL*8, DIMENSION(nc,LM)         :: y
      REAL*8, DIMENSION(n_rx,LM)       :: rr
      REAL*8, DIMENSION(n_bi)          :: pe, ea
      REAL*8, DIMENSION(n_tri)         :: ro, r1, sn, sb
      REAL*8, DIMENSION(n_bnd2,n_rx)   :: sigg
      REAL*8, DIMENSION(LM,n_oig)      :: conc
      REAL*8, DIMENSION(n_srb,p_5)     :: qfu
      REAL*8, DIMENSION(n_bnd3,LM)     :: qf
      REAL*8, DIMENSION(n_bnd3)        :: wlt, sO3, sO2
      REAL*8, DIMENSION(n_oig,n_bnd3)  :: sech
      REAL*8, DIMENSION(M__)           :: AFASTJ,C1,HFASTJ,V1
      REAL*8, DIMENSION(M__,M__)       :: BFASTJ,AAFASTJ,CC,SFASTJ,
     &                                    WFASTJ,U1
      REAL*8, DIMENSION(M__,2*M__)     :: PM
      REAL*8, DIMENSION(M__,M__,N__)   :: DD
      REAL*8, DIMENSION(M__,N__)       :: RR2    
      REAL*8, DIMENSION(2*M__)         :: PM0,dpomega
      REAL*8, DIMENSION(2*M__,N__)     :: POMEGA
      REAL*8, DIMENSION(N__)           :: ZTAU,FZ
      REAL*8, DIMENSION(2*M__,2*LM+2+1):: POMEGAJ
      REAL*8, DIMENSION(2*LM+2)        :: PFASTJ
      REAL*8, DIMENSION(NWFASTJ+1)     :: WBIN
      REAL*8, DIMENSION(NWFASTJ)       :: WL,FL,QRAYL,QBC,DUMMY,FLX
      REAL*8, DIMENSION(NWFASTJ,3)     :: QO3, QO2, Q1D, zpdep
      REAL*8, DIMENSION(3,NS)          :: TQQ
      REAL*8, DIMENSION(4,NP)          :: QAAFASTJ, WAAFASTJ
#ifdef SHINDELL_STRAT_CHEM
     &                                    ,SSA,RAA
#endif
      REAL*8, DIMENSION(8,4,NP)        :: PAA
      REAL*8, DIMENSION(31,18,12)      :: OREF
      REAL*8, DIMENSION(41,18,12)      :: TREF
      REAL*8, DIMENSION(41)            :: BREF
      REAL*8, DIMENSION(LM)            :: odtmp,ta,pres,TFASTJ,Jacet
      REAL*8, DIMENSION(NS)            :: VALJ
      REAL*8, DIMENSION(N__)           :: FJFASTJ
      REAL*8, DIMENSION(NWFASTJ,2,NS-3):: QQQ
      REAL*8, DIMENSION(NWFASTJ,jpnl)  :: FFF
      REAL*8, DIMENSION(JPPJ)          :: jfacta
      REAL*8, DIMENSION(JPNL,JPPJ)     :: zj, JFASTJ
      REAL*8, DIMENSION(p_2,LM)        :: chemrate, photrate 
      REAL*8, DIMENSION(2*LM)          :: O3_FASTJ    
      REAL*8, DIMENSION(ny,LM)         :: dest, prod
      REAL*8, DIMENSION(NLFASTJ,NLFASTJ):: WTAU
#ifdef SHINDELL_STRAT_CHEM
      REAL*8                            :: SF3_fact,SF2_fact
      REAL*8, DIMENSION(MXFASTJ,NBFASTJ):: AER2
      REAL*8, DIMENSION(NBFASTJ,NBFASTJ):: AMF
      REAL*8, DIMENSION(NBFASTJ)        :: TJ2,DO32,DBC2,ZFASTJ2,
     &                                     DMFASTJ2
      REAL*8, DIMENSION(LM+3)           :: PFASTJ2
      REAL*8, DIMENSION(51,18,12)       :: OREF2,TREF2
      REAL*8, DIMENSION(51)             :: BREF2
#endif
      REAL*8, DIMENSION(NLFASTJ)       :: aer,ZFASTJ,O3J,TJ,DBC,
     &  DMFASTJ,XQO3,XQO2,DTAUDZ,TTAU,FTAU,PIAER,RZ,RQ,DO3,PIRAY
      REAL*8, DIMENSION(LCOalt)        :: OxICINL,CH4ICINL
#ifdef SHINDELL_STRAT_CHEM
     &                                   ,N2OICINL,CFCICINL
#endif
      REAL*8, DIMENSION(LM)  :: COalt,CH4altT,CH4altX,OxICL,CH4ICL
#ifdef SHINDELL_STRAT_CHEM
     &                        ,BrOxalt,ClOxalt,ClONO2alt,HClalt,odcol                    
     &                        ,N2OICL,CFCICL
#endif

      LOGICAL                      :: fam,prnrts,prnchg,prnls      

      CHARACTER*20, DIMENSION(NP)  :: title_aer_pf !formerly TITLEA( )
      CHARACTER*78                 :: TITLE0
      CHARACTER*7, DIMENSION(3,NS) :: TITLEJ
      CHARACTER*7, DIMENSION(JPPJ) :: jlabel
      CHARACTER*7, DIMENSION(3)    :: lpdep
      CHARACTER*8, DIMENSION(nc)   :: ay
      
      COMMON/CHEM_LOC/chemrate,dest,FASTJLAT,FFF,O3_FASTJ,PFASTJ,
     & photrate,pres,prod,RFLECT,rr,SZA,ta,TANHT,TFASTJ,U0,VALJ,
     & WTAU,y,zj,jndlv,jndlev,jaddlv,jaddto,MIEDX,NCFASTJ!integers last
!$OMP THREADPRIVATE(/CHEM_LOC/)

      COMMON/FJAST_LOC/aer,ZFASTJ,O3J,TJ,DBC,DMFASTJ,XQO3,XQO2,DTAUDZ,
     & TTAU,FTAU,rr2,dd,PIAER,RZ,RQ,DO3,PIRAY,JFASTJ,odtmp,odsum,
     & XLTAU,dpomega,pomega,pomegaj,ztau,fz,zrefl,zu0,zflux,pm0,pm,
     & fjfastj,wfastj,BFASTJ,AFASTJ,AAFASTJ,CC,HFASTJ,C1,SFASTJ,U1,V1 
!$OMP THREADPRIVATE(/FJAST_LOC/)

#ifdef SHINDELL_STRAT_CHEM
      COMMON/SCHEM_LOC/ratioNs,rNO2frac,rNOfrac,rNOdenom,ratioN2
!$OMP THREADPRIVATE(/SCHEM_LOC/)

      COMMON/FJAST2_LOC/AER2,odcol,TJ2,DO32,DBC2,ZFASTJ2,
     &                  DMFASTJ2,PFASTJ2,AMF,jadsub
!$OMP THREADPRIVATE(/FJAST2_LOC/)
#endif
#endif
      END MODULE TRCHEM_Shindell_COM
      
      
      
      subroutine alloc_trchem_shindell_com(grid)
!@SUM  To alllocate arrays whose sizes now need to be determined
!@+    at run-time
!@auth G.Faluvegi
!@ver  1.0
      use domain_decomp, only : dist_grid, get
      use model_com, only     : im,lm
      use TRCHEM_Shindell_COM, only: DU_O3,corrOxIN,corrOx,ss,yNO3,
     & pHOx,pNOx,pOx,yCH3O2,yC2O3,yROR,yXO2,yAldehyde,yXO2N,yRXPAR,
     & TX,sulfate,OxIC,CH4ICX,dms_offline,so2_offline,yso2,ydms,
     & OxICIN,CH4ICIN,LcorrOx,JPPJ,LCOalt,acetone
#ifdef SHINDELL_STRAT_CHEM
     & ,pClOx,pClx,pOClOx,pBrOx,yCl2,yCl2O2,N2OICX,CFCIC,SF3,SF2,
     & N2OICIN,CFCICIN
#endif

      IMPLICIT NONE

      type (dist_grid), intent(in) :: grid
      integer :: ier, J_1H, J_0H
      logical :: init = .false.

      if(init)return
      init=.true.
    
      call get( grid , J_STRT_HALO=J_0H, J_STOP_HALO=J_1H )
 
      allocate(    corrOxIN(J_0H:J_1H,LcorrOx,12) )
      allocate(      corrOx(J_0H:J_1H,LM,12)      )
      allocate(          ss(JPPJ,LM,IM,J_0H:J_1H) )
      allocate(     acetone(IM,J_0H:J_1H,LM)      )
      allocate(        yNO3(IM,J_0H:J_1H,LM)      )
      allocate(        pHOx(IM,J_0H:J_1H,LM)      )
      allocate(        pNOx(IM,J_0H:J_1H,LM)      )
      allocate(         pOx(IM,J_0H:J_1H,LM)      )
      allocate(      yCH3O2(IM,J_0H:J_1H,LM)      )
      allocate(       yC2O3(IM,J_0H:J_1H,LM)      )
      allocate(        yROR(IM,J_0H:J_1H,LM)      )
      allocate(        yXO2(IM,J_0H:J_1H,LM)      )
      allocate(   yAldehyde(IM,J_0H:J_1H,LM)      )
      allocate(       yXO2N(IM,J_0H:J_1H,LM)      )
      allocate(      yRXPAR(IM,J_0H:J_1H,LM)      )
      allocate(          TX(IM,J_0H:J_1H,LM)      )
      allocate(     sulfate(IM,J_0H:J_1H,LM)      )
      allocate(        OxIC(IM,J_0H:J_1H,LM)      )
      allocate(      CH4ICX(IM,J_0H:J_1H,LM)      )
      allocate( dms_offline(IM,J_0H:J_1H,LM)      )
      allocate( so2_offline(IM,J_0H:J_1H,LM)      )
      allocate(        yso2(IM,J_0H:J_1H,LM)      )
      allocate(        ydms(IM,J_0H:J_1H,LM)      )
      allocate(      OxICIN(IM,J_0H:J_1H,LCOalt)  )
      allocate(     CH4ICIN(IM,J_0H:J_1H,LCOalt)  )
      allocate(       DU_O3(J_0H:J_1H)            )
#ifdef SHINDELL_STRAT_CHEM
      allocate(       pClOx(IM,J_0H:J_1H,LM)      )
      allocate(        pClx(IM,J_0H:J_1H,LM)      )
      allocate(      pOClOx(IM,J_0H:J_1H,LM)      ) 
      allocate(       pBrOx(IM,J_0H:J_1H,LM)      )
      allocate(        yCl2(IM,J_0H:J_1H,LM)      ) 
      allocate(      yCl2O2(IM,J_0H:J_1H,LM)      )
      allocate(      N2OICX(IM,J_0H:J_1H,LM)      )
      allocate(       CFCIC(IM,J_0H:J_1H,LM)      )
      allocate(         SF3(IM,J_0H:J_1H,LM)      )
      allocate(         SF2(IM,J_0H:J_1H,LM)      )
      allocate(     N2OICIN(IM,J_0H:J_1H,LCOalt)  )
      allocate(     CFCICIN(IM,J_0H:J_1H,LCOalt)  )
#endif
      
      return
      end subroutine alloc_trchem_shindell_com
      
