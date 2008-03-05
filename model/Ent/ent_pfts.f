      module ent_pfts
!@sum GISS vegetation types for Ent model.

      use ent_const
      use ent_types
      implicit none


      !***************************************************
      !*            GISS VEGETATION TYPES                *
      !***************************************************
      !Types: These are GISS GCM types until we get full data sets for Ent.
      !1-tundra, 2-grassland, 3-shrubland, 4-savanna, 5-deciduous forest,
      !6-evergreen needleleaf forest, 7-tropical rainforest, 8-crops
      !**** Correspondence to CASA/LSM vegetation types for sla from 
      !**** Dickinson et al. (J. Clim., Nov. 1998).  See casa_bgfluxes.f and
      !**** 
      !#-GISS   = Dickinson = CASA/LSM
      !1-tundra = tundra and semidesert
      !2-grassland = avg short grass & tall grass = LSM cool C3, warm C4 grasses
      !3-shrubland = evergreen shrubs, deciduous shrubs
      !4-savanna = shrubs & broadleaf evergreen & grasses (avg 25,35,40)
      !5-deciduous forest = deciduous broadleaf & decid needleleaf
      !6-evergreen needleleaf = needleleaf evergreen
      !7-tropical rainforest = deciduous broadleaf = LSM tropical seasonal tree 
      !8-crops = crops

      integer, parameter :: COVEROFFSET = 1 !SAND in first position in GISS array
      integer, parameter :: TUNDRA = 1
      integer, parameter :: GRASSC3 = 2
      integer, parameter :: SHRUB = 3
      integer, parameter :: SAVANNA = 4
      integer, parameter :: DECIDFOREST = 5
      integer, parameter :: EVERGRNEEDLE = 6
      integer, parameter :: TROPRAINF = 7
      integer, parameter :: CROPS = 8

      !*-----------------------------------------
      !* Veg types correspondence between models:
*       LSM:   1  2    3  4  5  6  7  8  9 10 11 12 13 14^M
*      CASA:   4  5    1  2  6  7  9 11 10 10 12 12  6  8^M
*      GISS:   6  5,6  7  5  4  2  3 x  1  1  8  8   4  x
*     NOTE:  Need to distinguish: - needleleaf decid vs. evergreen
*                                 - bare soil
*                                 - savanna mix
*            Check for:  sla, lrage, woodage, albedo
      !*-----------------------------------------
      !For dependencies to work, these constants must go in ent_const.f
!      integer,parameter :: N_PFT = 8
!      integer,parameter :: N_SOILCOV = 2 !light sand, dark dirt (GISS)
!      integer,parameter :: N_OTHER = 1
!      integer,parameter :: N_COVERTYPES = N_PFT + N_SOILCOV + N_OTHER
       !pst - photosynthetic pathway, 1-C3, 2-C4
      !hwilt - wilting point (m)
      !sstar - soil moisture stress onset point (fraction of soil volumetric saturation)
      !swilt - wilting point (fraction of soil volumetric saturation)
      !nf - canopy nitrogen factor (dimensionless)
      !* Parameters from CASA:
      !sla - specific leaf area (m2 leaf area/kg leaf C)(CASA)
      !r - CASA respiration parameter (gC/gN)(CASA)
      !lrage - leaf and root litter turnover time (years) (CASA)
      !woodage - stem litter turnover time (years) (CASA)
      !lit_C2N - litter C:N (CASA)
      !lignin - lignin content
      !*for GISS pft, 
      !*phenotype, b1Cf, b2Cf, b1Cd, b2Cd, b1Ht, b2Ht are dummies!
      !phenotype - phenological type 
      !* Parameters for plant allomteries 
      !* (Albani et al. Global Change Biology 2006 & 
      !* estimated from KM67 (Santarem, Amazon) tree survey.)
      !b1Cf - para 1 for allometric relation between DBH & foliage C 
      !b2Cf - para 2 for allometric relation between DBH & foliage C
      !b1Cd - para 1 for allometric relation between DBH & structural(dead) C
      !b2Cd - para 2 for allometric relation between DBH & structural(dead) C
      !b1Ht - para 1 for allometric relation between DBH & height
      !b2Ht - para 2 for allometric relation between DBH & height

      type(pftype),parameter :: pfpar(N_PFT) =          !PFT parameters
     &!        pst,hwilt,sstar,swilt,nf,sla,r,lrage,woodage,lit_C2N,lignin,phenotype ! 
     &!        b1Cf, b2Cf, b1Cd, b2Cd, b1Ht, b2Ht
     &     (/
     &     pftype(1,   -153.d0,  .50d0, .30d0,  1.4d0, !tundra
     &     22.5d0, 0.6d0, 2.8d0, 5.5d0, 50.0d0,0.15d0,1,
     &     0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0),
     &     pftype(2,   -2030.d0,  .30d0, .10d0,  1.5d0, !grassC3 !hwilt=Vaira grassland final senescence soilmp layer2
     &     41.1d0, 0.6d0, 1.5d0, UNDEF, 50.0d0, 0.1d0,4, !Reich 1997 Figure 3b and Ponca, but 41.1 sla for Reich equation cited by CLM, leaf longevity 0.5 yr
     &     0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0),
!     &     pftype(2,   -100.d0,  .45d0, .27d0,  2.9d0, !grassC3 !NK Ponca nf and SLA
!     &     10.0d0, 0.6d0, 1.5d0, UNDEF, 50.0d0, 0.1d0,4),
     &     pftype(2,   -153.d0,  .65d0, .22d0,  1.3d0, !shrub
     &     32.5d0, 1.d0, 1.25d0, 5.5d0, 57.5d0, 0.15d0,3,
     &     0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0),
     &     pftype(2,   -2030.d0,  .65d0, .22d0,  1.3d0, !savanna
     &     32.5d0, 1.d0, 1.8d0, 25.d0, 50.0d0, 0.15d0,3,
     &     0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0),
     &     pftype(1,   -153.d0,  .55d0, .29d0,  1.5d0, !decidforest
     &     30.d0, 0.85d0, 1.5d0, 42.5d0, 50.0d0, 0.2d0,2,
     &     0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0),
     &     pftype(1,   -153.d0,  .60d0, .25d0,  0.9d0, !evergrneedle
     &     18.d0, 0.9d0, 5.d0,42.0d0, 80.0d0,0.25d0,1,   !SLA-CLM 10, Reich 18 for Nleaf 2.8 gN/m2 give C:N 20
     &     0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0),
     &     pftype(1,   -153.d0,  .55d0, .26d0,  1.1d0, !troprainf
     &     25.d0, 0.733d0, 1.8d0, 41.0d0, 40.0d0, 0.2d0,1,
     &     0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0),
     &     pftype(2,   -153.d0,  .45d0, .27d0,  1.3d0, !crops
     &     60.d0, 0.6d0, 1.1d0, 58.0d0, 52.5d0, 0.16d0,4,
     &     0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0)
!     &     pftype(1,   -100.d0,  .50d0, .30d0,  0.76d0)&
!     &     0.0d0,0.0d0,0.0d0,0.0d0)
     &     /)

      !**NOTE:  Above, sstar and swilt are guesses for all except grassland
      ! and savanna.


      !***************************************************
      !*         END - GISS VEGETATION TYPES             *
      !***************************************************
      end module ent_pfts
