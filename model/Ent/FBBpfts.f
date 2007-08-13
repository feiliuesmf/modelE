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
      real*8 :: cf              !CO2 mole fraction at canopy surface (umol mol-1)
      real*8 :: ci              !Leaf internal CO2 mole fraction (umol mol-1)
      real*8 :: Tc              !Canopy (foliage) temperature (Celsius)
      real*8 :: Pa              !Atmospheric pressure (Pa)
      real*8 :: rh              !Relative humidity (fraction)
      end type psdrvtype


      !=======CONSTANTS========!
!      integer,parameter :: N_PFT = 8  !2-C3 herbaceous

      type(pspartype),parameter :: pftpar(N_PFT) = !PFT parameters for GISS veg types
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
     &     ,60.d0                !Vcmax, best fit guess, see plots of 08/10/2007.
     &     ,11.d0               !m, X.Mo, et al. (2001)
     &     ,.008d0              !b, X.Mo, et al. (2001)
!     &     ,18.72d0,46390.d0    !Rdc,RdH, Bernacchi, et al. (2001) Nicotiana tabacum
     &     ,3.27d0),            !Ponca Ntot/LA average (actually seasonal curve); Rd not large deviance from direct daily Ntot 
     &     pspartype(2          !SHRUB
     &     ,.9d0               !leaf VIS albedo,CLM C3 arctic grass, Table 3.1 (Oleson, et al 2004)
     &     ,17d0               !Vmax25, CLM 
     &     ,9d0                !m, CLM
     &     ,.002d0              !b, CLM
     &     ,2d0                 !Nleaf (gN/m2-leaf). Est. from Reich 1997 (big range).
     &     ),
     &     pspartype(2          !SAVANNA
     &     ,.9d0               !leaf VIS albedo,CLM C3 arctic grass, Table 3.1 (Oleson, et al 2004)
     &     ,24d0               !Vmax25, CLM
     &     ,5d0                !m, CLM
     &     ,.002d0              !b, CLM
     &     ,4d0                 !Nleaf (gN/m2-leaf). Est. from Reich 1997 (big range).
     &     ),
     &     pspartype(1          !DECIDFOREST
     &     ,0d0               !leaf VIS albedo,CLM C3 arctic grass, Table 3.1 (Oleson, et al 2004)
     &     ,0d0               !Vmax25, CLM
     &     ,0d0                !m, CLM
     &     ,0d0              !b, CLM
     &     ,0d0                 !Nleaf (gN/m2-leaf). Est. from Reich 1997 (big range).
     &     ),
     &     pspartype(1          !EVERGRNEEDLE
     &     ,0d0               !leaf VIS albedo,CLM C3 arctic grass, Table 3.1 (Oleson, et al 2004)
     &     ,0d0               !Vmax25, CLM
     &     ,0d0                !m, CLM
     &     ,0d0              !b, CLM
     &     ,0d0                 !Nleaf (gN/m2-leaf). Est. from Reich 1997 (big range).
     &     ),
     &     pspartype(1          !TROPRAINF
     &     ,0d0               !leaf VIS albedo,CLM C3 arctic grass, Table 3.1 (Oleson, et al 2004)
     &     ,0d0               !Vmax25, CLM
     &     ,0d0                !m, CLM
     &     ,0d0              !b, CLM
     &     ,0d0                 !Nleaf (gN/m2-leaf). Est. from Reich 1997 (big range).
     &     ),
     &     pspartype(2          !CROPS
     &     ,0d0               !leaf VIS albedo,CLM C3 arctic grass, Table 3.1 (Oleson, et al 2004)
     &     ,0d0               !Vmax25, CLM
     &     ,0d0                !m, CLM
     &     ,0d0              !b, CLM
     &     ,0d0                 !Nleaf (gN/m2-leaf). Est. from Reich 1997 (big range).
     &     )
     &/)

      !NOTES:
        !--------Collatz, et al. (1991) Farquhar parameters-----
        !Vmax values:
        !Collatz, C3 grass. Vcmax = 200.
        !Harley, et al. (1992), cotton Vcmax=51-127 umol m-2 s-1
        !Ponca, Oklahoma (Fluxnet), winter wheat, Vcmax fit ~30.
!        pspar%Vmax = 200./(1 + exp((-220.e03+703.*(Tl+Kelvin))



!****************************************************************************
      end module FarquharBBpspar
