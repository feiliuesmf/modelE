      MODULE TRCHEM_Shindell_COM
!@sum  TRCHEM_Shindell_COM declares variables for tracer chemistry
!@+    and sources.
!@auth Drew Shindell (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on various chemistry modules of B436Tds3YM23 model)
c
      USE MODEL_COM, only :  im,jm,lm,ls1,istrat,psf,
     &                       ptop,sig,sige,dsig,bydsig,
     &                       dtsrc,Itime,ItimeI,T,JEQ
      USE CONSTANT, only   : pi, mair, mwat, radian
      USE DYNAMICS, only   : am, byam, PMID, PK
      USE GEOM, only       : BYDXYP,dxyp
C     USE RADPAR, only     : LX
      USE TRACER_COM, only : ntm, trm, TR_MM
c
      IMPLICIT NONE
      SAVE
c
C**************  P  A  R  A  M  E  T  E  R  S  *******************
C
!@param p_1 number of reactants per reaction
!@param p_2 total number of reactions
!@param p_3 number of rxns in assembled lists (check with print rxn list)
!@param p_4 total number of reactions
!@param p_5 number of levels from top down with SRB flux
!@param n_rx maximum number of chemical reactions
!@param n_bi maximum number of bimolecular reactions
!@param n_tri maximum number of trimolecular reactions
!@param n_bnd1 maximum number of spectral bands 1
!@param n_bnd2 maximum number of spectral bands 2
!@param n_bnd3 maximum number of spectral bands 3
!@param n_nst maximum number of reverse reactions
!@param n_fam maximum number of chemical families
!@param n_spc maximum number of chemical species
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
!@param nc "number of molecules considered" ?
!@param ny "number of calculated gases" ?
!@param numfam number of families
!@param n_igas uknown: number of prod/dest gases?
!@param n_phot how often to do photolysis (in increments of DTsrc)
!@param O3MULT =2.14D-2    !is 1.E10*2.69E13*48./6.02E26
!@param BYO3MULT = 1/O3MULT
!@param pfix_O2 fixed ratio of O2/M
!@param pfix_H2 fixed ratio of H2/M
!@param pfix_Aldehyde fixed ratio of Aldehyde/M for initial conditions
!@param fix_CH4_chemistry logical whether or not to used a fixed 
!@+     value for methane in the chemistry code
!@param checktracer_on logical to turn on the checktracer call
!@param pfix_CH4_S fixed ratio of CH4/M in South. Hemis. (if used)
!@param pfix_CH4_S fixed ratio of CH4/M in South. Hemis. (if used)
!@param MWabyMWw ratio of molecular weights of air/water
!@param O3_1_fact factor to alter surface O3 that is passed to FASTJ
!@+     this is fastj level 1, not model level 1.  Currently, it is 
!@+     decreased by a factor of (972/1000)mb
!@param RGAMMASULF unknown
!@param RKBYPIM=8.*RBOLTZ/pi/MASSN2O55=8.*1.38062D-23/3.14159/1.793D-25
!@param cboltz Boltzman's Constant = 1.3806d-19
!@param dlogp 10.d0**(-2./16.)
!@param byradian 1/radian = conversion from radians to degrees
!@param LCOalt number of levels in the COaltIN array 
!@param LCH4alt number of levels in the CH4altIN array
!@param LcorrOx number of levels in the corrOxIN array
!@param JCOlat number of latitudes in the COlat array 
!@param PCOalt pressures at LCOalt levels
!@param PCH4alt pressures at LCH4alt levels
!@param PcorrOx pressures at LcorrOx levels
      INTEGER, PARAMETER ::
     & LCOalt =   23,
     & JCOlat =   19,
     & LCH4alt=    5,
     & LcorrOX=    4,
     & p_1   =     2,
     & p_2   =   111,
     & p_3   =   200,
     & p_4   =    70,
     & p_5   =    14,
     & n_rx  =   101,
     & n_bi  =    96,
     & n_tri =    14,
     & n_bnd1=    31,
     & n_nst =     5,
     & n_fam =     4,
     & n_spc =    35,
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
c     & n_CH3OOH=    6,    ! TRACERS_DRV.f.
c     & n_HCHO=      7,    ! Note the UNDERSCORE !!!
c     & n_HO2NO2=    8,    !  T
c     & n_CO=        9,    !  R
c     & n_CH4=      10,    !  A
c     & n_PAN=      11,    !  C
c     & n_Isoprene= 12,    !  E
c     & n_AlkylNit= 13,    !  R
c     & n_Alkenes=  14,    !  S
c     & n_Paraffin= 15,    ! ---------------
     & nC2O3=     16,
     & nXO2=      17,    
     & nXO2N=     18,
     & nRXPAR=    19,
     & nROR=      20,
     & nAldehyde= 21,
     & nH2O=      22,
     & nCH3O2=    23,
     & nH2=       24,
     & nOH=       25,
     & nHO2=      26,
     & nO3=       27,
     & nO=        28,
     & nO1D=      29,
     & nNO=       30,
     & nNO2=      31,
     & nNO3=      32,
     & nHONO=     33,
     & nO2=       34,
     & nM=        35,     !you must always put nM last (highest number)
C ----------------------------------------------   
     & N__=     1800,     !jan00, was 450, then 900 in Nov99
     & M__=        4,
     & NLFASTJ=  350,     !300 is arbitrary for now
     & NS     =   51,
     & NWFASTJ=   15, 
     & JPPJ   =   16,
     & JPNL   =   12,
     & NJVAL  =   16,     !formerly read in from jv_spec00_15.dat
     & MFIT   =    8,
     & NFASTJ =    4,
     & MFASTJ =    1,
     & nc     =   35,     !formerly in param sub
     & ny     =   33,     !formerly in param sub  
     & numfam =    2,     !formerly in param sub  
     & n_igas =   43,
     & n_phot=     2      ! currently means 2 hours
C     
      REAL*8, PARAMETER ::  O3MULT       = 2.14d-2,
     &                      BYO3MULT     = 1./O3MULT,
     &                      pfix_O2      = 0.209476d0,
     &                      pfix_H2      = 560.d-9,
     &                      pfix_Aldehyde= 2.d-9,
     &                      pfix_CH4_S   = 1.75d-6,
     &                      pfix_CH4_N   = 1.855d-6,
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
     &                      byradian     = 1.d0/radian 

!@dbparam ch4_init_sh,ch4_init_nh initial methane conc. (ppmv) 
!@+       defaults are for 1990
      REAL*8 :: ch4_init_sh=1.750d0    , ch4_init_nh=1.855d0

      REAL*8, PARAMETER, DIMENSION(LCOalt) :: PCOalt = (/
     & 0.9720D+03,0.9445D+03,0.9065D+03,
     & 0.8515D+03,0.7645D+03,0.6400D+03,0.4975D+03,0.3695D+03,
     & 0.2795D+03,0.2185D+03,0.1710D+03,0.1335D+03,0.1016D+03,
     & 0.7120D+02,0.4390D+02,0.2470D+02,0.1390D+02,0.7315D+01,
     & 0.3045D+01,0.9605D+00,0.3030D+00,0.8810D-01,0.1663D-01/)
      REAL*8, PARAMETER, DIMENSION(LcorrOx) :: PcorrOx = 
     & (/133.5d0, 101.6d0, 71.2d0, 43.9d0/)
C  
C    These should really be defined in the run deck:
      LOGICAL, PARAMETER :: luselb            = .false.,
     &                      fix_CH4_chemistry = .false.,
     &                      checktracer_on    = .false.
c
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
!@+   1=Rayly 2=iso 30iso-equiv 4=bkgrd-sulf,5=volc-sulf,6=liq water
!@var NAA Number of categories for scattering phase functions
!@var npdep Number of pressure dependencies
!@var jpdep Index of cross sections requiring P dependence
!@var PFASTJ pressure sent to FASTJ
!@var nss this is a copy of JPPJ that is read in from a file
!@var jfacta Quantum yield (or multiplication factor) for photolysis
!@var WBIN Boundaries of wavelength bins
!@var WL Centres of wavelength bins - 'effective wavelength'
!@var NWWW Number of wavelength bins, from NW1:NW2
!@var FL Solar flux incident on top of atmosphere (cm-2.s-1)
!@var   Rayleigh parameters (effective cross-section) (cm2)
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
!@var U0 cosine of the solar zenith angle
!@var ZFASTJ Altitude of each fastj pressure level (approx.) (cm)
!@var RFLECT Surface albedo (Lamertian) in fastj
!@var odtmp Optical depth (temporary array)
!@var odsum Column optical depth
!@var nlbatm Level of lower photolysis boundary - usually surface ('1')
!@var aer obs fastj aerosol profile?
!@var O3J Ozone profile on photolysis grid
!@var TJ Temperature profile on photolysis grid
!@var DBC Mass of Black Carbon at each pressure level (g.cm-3)
!@var FFF Actinic flux at each level for each wavelength bin and level
!@var DMFASTJ Total number density at each pressure level (cm-3)
!@var DO3 Ozone number density at each pressure level (cm-3
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
!@var O3DLJI model ozone to radiation (calculated)
!@var O3DLJI_clim model ozone to radiation (climatological)
!@var COlat carbon monoxide latitude distribution for init cond (ppbv)
!@var COaltIN adjustments of COlat by altitude (unitless,LCOalt levels)
!@var CH4altINT tropical strat adjustments to CH4 (LCH4alt levels)
!@var CH4altINX xtra-tropical strat adjustments to CH4 LCH4alt levels)
!@var COalt adjustments of COlat by altitude (unitless,LM levels)
!@var CH4altT tropical strat adjustments to CH4 (unitless, LM levels)
!@var CH4altX xtra-tropical strat adjustments to CH4 (LM levels)
!@var BYFJM = 1/JM
!@var mass2vol local array to convert between mass and volume units.
!@var bymass2vol local array to convert between mass and volume units.
!@var MODPHOT if MODPHOT=0 do photolysis, else skip it
!@var TX temperature variable for master chem
!@var ta, pres local arrays to hold temperature,pressure
!@var TFASTJ temperature sent to FASTJ
!@var O3_FASTJ ozone sent to fastj
!@var FASTJLAT,FASTJLON latitude & LONGITUDE (degrees) for use in fastj
!@var SZA the solar zenith angle (degrees)
!@var JFASTJ photolysis rates
!@var RVELN2O5, prod_sulf unknown
!@var sulfate N2O5 sulfate sink (formerly SRC(I,J,L,20) variable)   
!@var prod_sulfate  N2O5 change by sulfate reactions in mass units
!@var wprod_sulf N2O5 change by sulfate reactions in molecules/cm3/s
!@var DT2 variable chemical time step, set in masterchem
!@var nr total number of        reactions read in from gs_jpl00_trop_15
!@var nr3 #of trimolecular      reactions read in from gs_jpl00_trop_15
!@var nr2 #of mono+bi-molecular reactions read in from gs_jpl00_trop_15
!@var nmm #of monomolecular     reactions read in from gs_jpl00_trop_15
!@var nhet #of heterogenous     reactions read in from gs_jpl00_trop_15
!@var pfactor to convert units on species chemical changes
!@var bypfactor to convert units on species chemical changes
!@var dNO3,gwprodHNO3,gprodHNO3,gwprodN2O5,changeAldehyde,
!@+   changeAlkenes,changeIsoprene,changeHCHO,changeAlkylNit,
!@+   changeHNO3,changeNOx,changeN2O5,wprodHCHO working variables to 
!@+   calculate nighttime chemistry changes
!@var rlossN,rprodN,ratioN variables for nitrogen conservation
!@var change change due to chemistry in mass/time
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
      INTEGER nr,nr2,nr3,nmm,nhet,MODPHOT,L75P,L75M,L569P,L569M, 
     & lprn,jprn,iprn,NW1,NW2,MIEDX,NAA,npdep,nss,
     & NWWW,NK,nlbatm,NCFASTJ                                 
      INTEGER, DIMENSION(n_fam)        :: nfam = (/27,30,0,0/) ! why 4?
      INTEGER, DIMENSION(p_1,p_2)      :: nn, nnr, kss
      INTEGER, DIMENSION(p_2)          :: ks
      INTEGER, DIMENSION(p_3)          :: nps, nds, npnr, ndnr
      INTEGER, DIMENSION(p_4)          :: kps, kds, kpnr, kdnr
      INTEGER, DIMENSION(n_nst)        :: nst
      INTEGER, DIMENSION(n_bnd1)       :: lbeg, nir
      INTEGER, DIMENSION(LM)           :: jndlv,jndlev
      INTEGER, DIMENSION(JPPJ)         :: jind
      INTEGER, DIMENSION(NLFASTJ)      :: jaddlv, jaddto
      INTEGER, DIMENSION(NJVAL)        :: jpdep  
      INTEGER, DIMENSION(12) , PARAMETER :: MDOFM =
     *     (/31,59,90,120,151,181,212,243,273,304,334,365/)
C      
      CHARACTER*8, DIMENSION(n_spc)    :: ay
C      
      REAL*8 ZFLUX,ZREFL,ZU0,U0,RFLECT,odsum,XLTAU,TANHT,
     & BYFJM,FASTJLAT,FASTJLON,SZA,RVELN2O5,
     & prod_sulf,DT2,wprod_sulf,dNO3,gwprodHNO3,gprodHNO3,gwprodN2O5,
     & changeAldehyde,changeAlkenes,changeIsoprene,changeHCHO,
     & wprodHCHO,changeAlkylNit,changeHNO3,changeNOx,changeN2O5,
     & wprodCO,rlossN,rprodN,ratioN,pfactor,bypfactor,F75P,F75M,
     & F569P,F569M
      REAL*8, DIMENSION(JM,LcorrOx,12) :: corrOxIN ! 12=month
      REAL*8, DIMENSION(JM,LM,12)      :: corrOx
      REAL*8, DIMENSION(n_spc,LM)      :: y
      REAL*8, DIMENSION(n_rx,LM)       :: rr
      REAL*8, DIMENSION(JPPJ,LM,IM,JM) :: ss
      REAL*8, DIMENSION(n_bi)          :: pe, ea
      REAL*8, DIMENSION(n_tri)         :: ro, r1, sn, sb
      REAL*8, DIMENSION(n_bnd2,n_rx)   :: sigg
      REAL*8, DIMENSION(LM,n_oig)      :: conc
      REAL*8, DIMENSION(n_srb,p_5)     :: qfu
      REAL*8, DIMENSION(n_bnd3,LM)     :: qf
      REAL*8, DIMENSION(n_bnd3)        :: wlt, sO3, sO2
      REAL*8, DIMENSION(n_oig,n_bnd3)  :: sech
      REAL*8, DIMENSION(IM,JM,LM)   :: yNO3,pHOx,pNOx,pOx,yCH3O2,yC2O3,
     &                                yROR,yXO2,yAldehyde,yXO2N,yRXPAR,
     &                                TX,sulfate
      REAL*8, DIMENSION(LM,IM,JM)      :: RCLOUDFJ
      REAL*8, DIMENSION(M__)         :: AFASTJ,C1,HFASTJ,V1
      REAL*8, DIMENSION(M__), PARAMETER ::
     *     EMU = (/.06943184420297D0, .33000947820757D0,
     *             .66999052179243D0,.93056815579703D0/), 
     *     WTFASTJ = (/.17392742256873D0,.32607257743127D0,
     *                 .32607257743127D0,.17392742256873D0/)
      REAL*8, DIMENSION(M__,M__)  :: BFASTJ,AAFASTJ,CC,SFASTJ,WFASTJ,U1
      REAL*8, DIMENSION(M__,2*M__)     :: PM
      REAL*8, DIMENSION(M__,M__,N__)   :: DD
      REAL*8, DIMENSION(M__,N__)       :: RR2    
      REAL*8, DIMENSION(2*M__)         :: PM0,dpomega
      REAL*8, DIMENSION(2*M__,N__)     :: POMEGA
      REAL*8, DIMENSION(N__)           :: ZTAU,FZ
      REAL*8, DIMENSION(2*M__,2*LM+2+1):: POMEGAJ
      REAL*8, DIMENSION(2*LM+2)        :: PFASTJ
      REAL*8, DIMENSION(NWFASTJ+1)     :: WBIN
      REAL*8, DIMENSION(NWFASTJ)       :: WL, FL, QRAYL, QBC
      REAL*8, DIMENSION(NWFASTJ,3)     :: QO3, QO2, Q1D, zpdep
      REAL*8, DIMENSION(3,NS)          :: TQQ
      REAL*8, DIMENSION(4,9)           :: QAAFASTJ, WAAFASTJ
      REAL*8, DIMENSION(8,4,9)         :: PAA
      REAL*8, DIMENSION(31,18,12)      :: OREF
      REAL*8, DIMENSION(41,18,12)      :: TREF
      REAL*8, DIMENSION(41)            :: BREF
      REAL*8, DIMENSION(NLFASTJ)       ::aer,ZFASTJ,O3J,TJ,DBC,DMFASTJ,
     &                                    XQO3,XQO2,DTAUDZ,TTAU,FTAU,
     &                                    PIAER,RZ,RQ,DO3,PIRAY
      REAL*8, DIMENSION(LM)            ::odtmp,ta,pres,TFASTJ  !,TXL
C Multiplier of free trop [CO] by layer Badr & Probert 94 fig10 & 11,
C Lopez-Valverde et al 93 fig 3, and Warneck 88, ch1 fig14 :
      REAL*8, DIMENSION(LCOalt), PARAMETER ::  COaltIN = 
     *     (/2d0,1.5625d0,1.375d0,1.25d0,1.125d0,1.0625d0,1d0,1d0,1d0
     *     ,1d0,1d0,.5d0,.375d0,.2d0,.2d0,.2d0,.2d0,.2d0,.25d0,.4d0,
     *     2.5d0,12d0,60d0/)
C**** additional levsls for CH4 to avoid extrapolation...
      REAL*8, PARAMETER, DIMENSION(LCH4alt) :: PCH4alt = 
     &     (/569d0, 150d0, 100d0, 32d0, 3.2d0, 0.23d0/)
      REAL*8, PARAMETER, DIMENSION(LCH4alt) ::   
     *   CH4altINT =(/1.79d0, 1.75d0, 1.620d0,1.460d0,0.812d0,0.230d0/),
     *   CH4altINX =(/1.79d0, 1.75d0, 1.440d0,1.130d0,0.473d0,0.202d0/)    

      REAL*8, DIMENSION(LM)            :: COalt,CH4altT,CH4altX
      REAL*8, DIMENSION(NS)            :: VALJ
      REAL*8, DIMENSION(N__)           :: FJFASTJ
      REAL*8, DIMENSION(NWFASTJ,2,NS-3):: QQQ
      REAL*8, DIMENSION(NWFASTJ,jpnl)  :: FFF
      REAL*8, DIMENSION(JPPJ)          :: jfacta
      REAL*8, DIMENSION(JPNL,JPPJ)      :: zj, JFASTJ
      REAL*8, DIMENSION(p_2,LM)         :: chemrate, photrate 
      REAL*8, DIMENSION(IM,JM)          :: SALBFJ
C     REAL*8, DIMENSION(LX,JM,IM)       :: O3DLJI, O3DLJI_clim
      REAL*8, DIMENSION(LM,JM,IM)       :: O3DLJI, O3DLJI_clim
      REAL*8, DIMENSION(2*(LS1-1))      :: O3_FASTJ
C [CO] ppbv based on 10deg lat-variation Badr & Probert 1994 fig 9:
      REAL*8, DIMENSION(JCOlat), PARAMETER  :: COlat = (/40.,40.,40.,40.
     *     ,45.,50.,60.,70.,80.,90.,110.,125.,140.,165.,175.,180.,170.
     *     ,165.,150./)
      REAL*8, DIMENSION(n_igas,LM)      :: dest, prod
      REAL*8, DIMENSION(NTM)            :: mass2vol,bymass2vol
      REAL*8, DIMENSION(IM,JM,LM,ntm)   :: change
      REAL*8, DIMENSION(NLFASTJ,NLFASTJ):: WTAU
C     
      LOGICAL                         fam,prnrts,prnchg,prnls      
C      
      CHARACTER*20, DIMENSION(9)   :: title_aer_pf !formerly TITLEA( )
      CHARACTER*78                    TITLE0
      CHARACTER*7, DIMENSION(3,NS) :: TITLEJ
      CHARACTER*7, DIMENSION(JPPJ) :: jlabel
      CHARACTER*7, DIMENSION(3)    :: lpdep
C
      COMMON/CHEM_LOC/bypfactor
     & ,changeAldehyde,changeAlkenes,changeAlkylNit
     & ,changeIsoprene,changeHCHO,changeHNO3,changeNOx,changeN2O5
     & ,chemrate,dest,dNO3,FASTJLAT
     & ,FFF,gwprodHNO3,gprodHNO3,gwprodN2O5,O3_FASTJ,PFASTJ,photrate
     & ,pfactor,pres,prod_sulf,prod,RFLECT,rlossN,rprodN,ratioN,rr
     & ,RVELN2O5,SZA
     & ,ta,TANHT,TFASTJ,U0,VALJ,wprod_sulf,wprodHCHO,wprodCO
     & ,WTAU,y,zj,jndlv
     & ,jndlev,jaddlv,jaddto,MIEDX,NCFASTJ  ! integers last (alignment)
!$OMP THREADPRIVATE(/CHEM_LOC/)

      COMMON/FJAST_LOC/aer,ZFASTJ,O3J,TJ,DBC,DMFASTJ,
     &     XQO3,XQO2,DTAUDZ,TTAU,FTAU,rr2,dd,
     &     PIAER,RZ,RQ,DO3,PIRAY , JFASTJ,odtmp,odsum,XLTAU,dpomega
     *     ,pomega,pomegaj,ztau,fz,zrefl,zu0,zflux,pm0,pm,fjfastj
     *     ,wfastj,BFASTJ,AFASTJ,AAFASTJ,CC,HFASTJ,C1,SFASTJ,U1,V1 
!$OMP THREADPRIVATE(/FJAST_LOC/)


      END MODULE TRCHEM_Shindell_COM
