      module FarquharBBpspar
      !@sum pfts  Plant functional type parameters for
      !Farqhuar-von Caemmerer (1982) photosynthesis and 
      !Ball-Berry (1985) stomatal conductance.
      use ent_const, only : N_PFT

      implicit none

      !======DECLARED TYPES====== !
      type pspartype
      integer :: pst            !Photosynth type.  1-C3, 2=C4
      real*8 :: PARabsorb       !Leaf PAR absorptance (fraction)

      !Photosynthesis/Conductance - Farquhar/Ball-Berry parameters
      real*8 :: Vcmax           !Maximum photosynthetic capacity (umol m-2 s-1)
!      real*8 :: Kc              !Michaelis-Menten constant for CO2 (Pa)
!      real*8 :: Ko              !Michaelis-Menten constant for O2 (Pa)
!      real*8 :: KcQ10           !Kc Q10 exponent
!      real*8 :: KoQ10           !Ko Q10 exponent
      !real*8 :: GammastarQ10    !CO2 compensation point Q10 (Pa)
      real*8 :: m               !Slope of Ball-Berry equation
      real*8 :: b               !Intercept of Ball-Berry equation (mol m-2 s-1)
!      real*8 :: Rdc             !Dark respiration c scaling factor (Harley&Tenhunen, 1991)
!      real*8 :: RdH             !Dark respiration deltaH (J mol-1) (Harley&Tenhunen, 1991)
      real*8 :: Nleaf           !g-N/m2[leaf]
      end type pspartype


      type psdrvtype
      real*8 :: ca              !Surface CO2 mole fraction (umol mol-1)
      real*8 :: ci              !Leaf internal CO2 mole fraction (umol mol-1)
      real*8 :: Tc              !Canopy (foliage) temperature (Celsius)
      real*8 :: Pa              !Atmospheric pressure (Pa)
      real*8 :: rh              !Relative humidity (fraction)
      end type psdrvtype


      !=======CONSTANTS========!
!      integer,parameter :: N_PFT = 13  !2-C3 herbaceous

!*********************************************************************
!* Ent PFTs
!* 1.  evergreen broadleaf early successional	
!* 2.  evergreen broadleaf late successional	
!* 3.  evergreen needleleaf early successional	
!* 4.  evergreen needleleaf late successional	
!* 5.  cold deciduous broadleaf early successional	
!* 6.  cold deciduous broadleaf late successional	
!* 7.  drought deciduous broadleaf	
!* 8.  decidous needleleaf	
!* 9.  cold adapted shrub
!* 10.  arid adapted shrub
!* 11.  C3 grass - perennial
!* 12.  C4 grass - perennial
!* 13.  C3 grass - annual
!* 14.  arctic C3 grass
!* 15.  crops - C4 herbaceous
!* 16.  crops - broadleaf woody
!*********************************************************************

#ifdef PSPAR_HACK
      type(pspartype),parameter :: pftpar(N_PFT) = !?PFT parameters for GISS veg types
     &!     pft PARabsorb Vcmax Kc Ko KcQ10 KoQ10 Gammastar  m b !Rdc RdH
     &     (/
     &     pspartype(1          !TUNDRA
     &     ,.89d0               !from leaf VIS albedo,CLM C3 arctic grass, Table 3.1 (Oleson, et al 2004)
     &     ,43.d0               !Vmax25, CLM
     &     ,9.d0                !m, CLM
     &     ,.002d0              !b, CLM
     &     ,7d0                 !Nleaf (gN/m2-leaf). Est. from Reich 1997 (big range).
     &),
     &     pspartype(1          !GRASSC3 - this is 2 in Ent, most grassland is C4.
     &     ,.86d0               !PARabsorb, Collatz et al. (1991)
!     &     ,30d0                !Vcmax, Wullschleger (1993) winter wheat, Triticum aestivum
!     &     ,100d0                !Vcmax, Sellers II (1996)
!     &     ,60d0                !Vcmax, von Caemmerer, CSIRO 2000, VARIOUS VALUES
!     &     ,45.d0                !Vcmax, best fit guess, see plots of 07/09/2007.
     &     ,60.d0               !Vcmax, best fit guess, see plots of 08/10/2007.
!     &     ,93.d0               !Vcmax, S. Verma and J. Berry,http://nigec.ucdavis.edu/publications/annual97/greatplains/project86.html
     &     ,11.d0               !m, X.Mo, et al. (2001)
     &     ,.008d0              !b, X.Mo, et al. (2001)
!     &     ,18.72d0,46390.d0    !Rdc,RdH, Bernacchi, et al. (2001) Nicotiana tabacum
     &     ,3.27d0),            !Ponca Ntot/LA average (actually seasonal curve); Rd not large deviance from direct daily Ntot 
     &     pspartype(2          !SHRUB
     &     ,.9d0               !1-leaf VIS albedo,CLM C3 arctic grass, Table 3.1 (Oleson, et al 2004)
     &     ,17d0               !Vmax25, CLM 
     &     ,9d0                !m, CLM
     &     ,.002d0              !b, CLM
     &     ,2d0                 !Nleaf (gN/m2-leaf). Est. from Reich 1997 (big range).
     &     ),
     &     pspartype(2          !SAVANNA
     &     ,.9d0               !1-leaf VIS albedo,CLM C3 arctic grass, Table 3.1 (Oleson, et al 2004)
     &     ,24d0               !Vmax25, CLM
     &     ,5d0                !m, CLM
     &     ,.002d0              !b, CLM
     &     ,4d0                 !Nleaf (gN/m2-leaf). Est. from Reich 1997 (big range).
     &     ),
     &     pspartype(2          !SHRUB
     &     ,.9d0               !1-leaf VIS albedo,CLM C3 arctic grass, Table 3.1 (Oleson, et al 2004)
     &     ,17d0               !Vmax25, CLM 
     &     ,9d0                !m, CLM
     &     ,.002d0              !b, CLM
     &     ,2d0                 !Nleaf (gN/m2-leaf). Est. from Reich 1997 (big range).
     &     ),
     &     pspartype(2          !SHRUB
     &     ,.9d0               !1-leaf VIS albedo,CLM C3 arctic grass, Table 3.1 (Oleson, et al 2004)
     &     ,17d0               !Vmax25, CLM 
     &     ,9d0                !m, CLM
     &     ,.002d0              !b, CLM
     &     ,2d0                 !Nleaf (gN/m2-leaf). Est. from Reich 1997 (big range).
     &     ),
     &     pspartype(2          !SHRUB
     &     ,.9d0               !1-leaf VIS albedo,CLM C3 arctic grass, Table 3.1 (Oleson, et al 2004)
     &     ,17d0               !Vmax25, CLM 
     &     ,9d0                !m, CLM
     &     ,.002d0              !b, CLM
     &     ,2d0                 !Nleaf (gN/m2-leaf). Est. from Reich 1997 (big range).
     &     ),
     &     pspartype(2          !SHRUB
     &     ,.9d0               !1-leaf VIS albedo,CLM C3 arctic grass, Table 3.1 (Oleson, et al 2004)
     &     ,17d0               !Vmax25, CLM 
     &     ,9d0                !m, CLM
     &     ,.002d0              !b, CLM
     &     ,2d0                 !Nleaf (gN/m2-leaf). Est. from Reich 1997 (big range).
     &     )
     &/)
#else !PSPAR_GISS
      type(pspartype),parameter :: pftpar(N_PFT) = !PFT parameters for GISS veg types
     &!     pft PARabsorb Vcmax Kc Ko KcQ10 KoQ10 Gammastar  m b !Rdc RdH
     &     (/
     &     pspartype(1          !1. TUNDRA
     &     ,.89d0               !from leaf VIS albedo,CLM C3 arctic grass, Table 3.1 (Oleson, et al 2004)
     &     ,43.d0               !Vmax25, CLM
     &     ,9.d0                !m, CLM
     &     ,.002d0              !b, CLM
     &     ,7d0                 !Nleaf (gN/m2-leaf). Est. from Reich 1997 (big range).
     &),
     &     pspartype(1          !2. GRASSC3 - this is 2 in Ent, most grassland is C4.
     &    ,.86d0                !PARabsorb, Collatz et al. (1991)
!     &     ,30d0                !Vcmax, Wullschleger (1993) winter wheat, Triticum aestivum
!     &     ,100d0                !Vcmax, Sellers II (1996)
!     &     ,60d0                !Vcmax, von Caemmerer, CSIRO 2000, VARIOUS VALUES
!     &     ,43.d0                !Vmax25, CLM , Table 8.2 (Oleson, et al 2004)
!     &     ,56.4d0                !Vcmax, Wang et al (2007) GCB, from inversion of flux data for oak-grass. Xu blue oak leaf chamber Vcmax reaches ~100.
     &     ,93.d0               !Vcmax, S. Verma and J. Berry,http://nigec.ucdavis.edu/publications/annual97/greatplains/project86.html
     &     ,11.d0               !m, X.Mo, et al. (2001)
     &     ,.008d0              !b, X.Mo, et al. (2001)
!     &     ,18.72d0,46390.d0    !Rdc,RdH, Bernacchi, et al. (2001) Nicotiana tabacum
     &     ,3.27d0              !Ponca Ntot/LA average (seasonal curve); Rd not large deviance from direct daily Ntot 
     &     ),
     &     pspartype(2          !3. SHRUB
     &     ,.9d0               !leaf VIS 1-albedo,CLM BES & BDS temperate 0.93, BDS 0.90, Table 3.1 (Oleson, et al 2004)
     &     ,17d0               !Vmax25, CLM BES & BDS temperate 17, BDS boreal 33, Table 8.2 (Oleson, et al 2004)
     &     ,9d0                !m, CLM
     &     ,.002d0              !b, CLM
     &     ,2d0                 !Nleaf (gN/m2-leaf). Est. from Reich 1997 (big range).
     &     ),
     &     pspartype(2          !4. SAVANNA
     &     ,.89d0               !leaf VIS 1-albedo,CLM BDT 0.90, grass 0.89, Table 3.1 (Oleson, et al 2004)
     &     ,32d0                !Vmax25, CLM C4 grass 24, BDT tropical 40,  Table 8.2 (Oleson, et al 2004)
!     &     ,47.d0              !Vmax25, CLM C3 grass 43, BDT temperate 51, Table 8.2 (Oleson, et al 2004)
!     &     ,56.4.d0             !Vcmax, Wang et al (2007) GCB, from inversion of flux data for oak-grass.
     &     ,9d0                !m, CLM
     &     ,.002d0              !b, CLM
     &     ,4d0                 !Nleaf (gN/m2-leaf). Est. from Reich 1997 (big range).
     &     ),
     &     pspartype(1          !5. DECIDFOREST
     &     ,0.90d0               !leaf VIS 1-albedo,CLM BDT temperate, Table 3.1 (Oleson, et al 2004)
     &     ,51.d0               !Vmax25, CLM BDT tropical 40, BDT temperate & boreal 51, Table 8.2 (Oleson, et al 2004)
     &     ,9.d0               !m, CLM BDT all, Table 8.2 (Oleson, et al 2004)
     &     ,.002d0              !b, CLM (Oleson, et al 2004, Section 8, p. 129)
     &     ,6.7d0                 !Nleaf (gN/m2-leaf). Est. from Reich 1997 (big range): (20 mg-N/g-leaf)/(300 cm2/g-leaf) =6.7 g-N/m2-leaf
     &     ),
     &     pspartype(1          !6. EVERGRNEEDLE
     &     ,0.93d0               !leaf VIS 1-albedo,CLM NET all, Table 3.1 (Oleson, et al 2004)
     &     ,47.d0               !Vmax25, CLM NET temperate 51, NET boreal 43, Table 8.2 (Oleson, et al 2004)
     &     ,6.d0               !m, CLM NET all, Table 8.2 (Oleson, et al 2004)
     &     ,.002d0              !b, CLM (Oleson, et al 2004, Section 8, p. 129)
     &     ,2.8d0                 !Nleaf (gN/m2-leaf). Est. from Reich 1997: (14 mg-N/g-leaf)/(50 cm2/g-leaf) = 2.8 g-N/m2-leaf
     &     ),
     &     pspartype(1          !7. TROPRAINF
     &     ,0.90d0               !leaf VIS 1-albedo,CLM BET & BDT tropical, Table 3.1 (Oleson, et al 2004)
     &     ,75.d0               !Vmax25, CLM BET tropical 75, BDT tropical 40, Table 8.2 (Oleson, et al 2004)
     &     ,9.d0               !m, CLM BET & BDT all, Table 8.2 (Oleson, et al 2004)
     &     ,.002d0              !b, CLM (Oleson, et al 2004, Section 8, p. 129)
     &     ,8.0d0                 !Nleaf (gN/m2-leaf). Est. from Reich 1997 (big range): (16 mg-N/g-leaf)/(200 cm2/g-leaf)=8 g-N/m2-leaf
     &     ),
     &     pspartype(2          !8. CROPS
     &     ,0.89d0               !leaf VIS 1-albedo,CLM Crop1 & Crop2, Table 3.1 (Oleson, et al 2004)
     &     ,50.d0               !Vmax25, CLM Crop1, Table 8.2 (Oleson, et al 2004)
     &     ,9.d0                !m, CLM Crop1, Table 8.2 (Oleson, et al 2004)
     &     ,.002d0              !b, CLM (Oleson, et al 2004, Section 8, p. 129)
     &     ,10.d0               !Nleaf (gN/m2-leaf). Est. from Reich 1997 (big range): (42 mg-N/g-leaf)/(400 cm2/g-leaf) =10.5 g-N/m2-leaf
     &     )
     &/)
#endif


      !NOTES:
        !--------Collatz, et al. (1991) Farquhar parameters-----
        !Vmax values:
        !Collatz, C3 grass. Vcmax = 200.
        !Harley, et al. (1992), cotton Vcmax=51-127 umol m-2 s-1
        !Ponca, Oklahoma (Fluxnet), winter wheat, Vcmax fit ~30.
!        pspar%Vmax = 200./(1 + exp((-220.e03+703.*(Tl+Kelvin))



!****************************************************************************
      end module FarquharBBpspar
