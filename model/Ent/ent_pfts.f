      module ent_pfts
!@sum GISS vegetation types for Ent model.

      use ent_const
      use ent_types
      implicit none


      !***************************************************
      !*            GISS VEGETATION TYPES                *
      !***************************************************
      !Types: These are GISS GCM types until we get full data sets for Ent.
      !       In addition, boreal forest from BOREAS SOBS is added.
      !1-tundra, 2-grassland, 3-shrubland, 4-savanna, 5-deciduous forest,
      !6-evergreen needleleaf forest, 7-tropical rainforest, 8-crops
      !9-boreal forest
      !For dependencies to work, these constants must go in ent_const.f
!      integer,parameter :: N_PFT = 8
!      integer,parameter :: N_SOILCOV = 2 !light sand, dark dirt (GISS)
!      integer,parameter :: N_OTHER = 1
!      integer,parameter :: N_COVERTYPES = N_PFT + N_SOILCOV + N_OTHER
       !pst - photosynthetic pathway, 1-C3, 2-C4
      !hwilt - wilting point (m)
      !sstar - soil moisture stress onset point (fraction of soil volumetric saturation)
      !swilt - wilting point (fraction of soil volumetric saturation)
      !nf - canopy nitrogen factor (dimnesionless)
      !sla - specific leaf area (m2 leaf area/kg leaf C)(CASA)
      !stressCD - cold and drought stress factor (dimensionless) (CASA)
      !lrage - leaf and root litter age (years) (CASA)
      !woodage - stem litter age (years) (CASA)
      !lit_C2N - litter C:N (CASA)
      !lignin - lignin content

      type(pftype),parameter :: pfpar(N_PFT) =          !PFT parameters
     &!        pst,hwilt,sstar,swilt,nf,sla,stressCD, lrage,woodage,lit_C2N,lignin ! 
     &     (/
     &     pftype(1,   -100.d0,  .50d0, .30d0,  1.4d0,
     &     0.d0, 0.d0, 0.0d0,0.0d0,0.0d0,0.0d0),
     &     pftype(2,   -100.d0,  .45d0, .27d0,  1.5d0,
     &     0.d0, 0.d0, 0.0d0,0.0d0,0.0d0,0.0d0),
     &     pftype(2,   -100.d0,  .65d0, .22d0,  1.3d0,
     &     0.d0, 0.d0, 0.0d0,0.0d0,0.0d0,0.0d0),
     &     pftype(2,   -100.d0,  .65d0, .22d0,  1.3d0,
     &     0.d0, 0.d0, 0.0d0,0.0d0,0.0d0,0.0d0),
     &     pftype(1,   -100.d0,  .55d0, .29d0,  1.5d0,
     &     0.d0, 0.d0, 0.0d0,0.0d0,0.0d0,0.0d0),
     &     pftype(1,   -100.d0,  .60d0, .25d0,  0.9d0,
     &     0.d0, 0.d0, 0.0d0,0.0d0,0.0d0,0.0d0),
     &     pftype(1,   -100.d0,  .55d0, .26d0,  1.1d0,
     &     0.d0, 0.d0, 0.0d0,0.0d0,0.0d0,0.0d0),
     &     pftype(2,   -100.d0,  .45d0, .27d0,  1.3d0,
     &     0.d0, 0.d0, 0.0d0,0.0d0,0.0d0,0.0d0)
!     &     pftype(1,   -100.d0,  .50d0, .30d0,  0.76d0)&
!     &     0.0d0,0.0d0,0.0d0,0.0d0)
     &     /)

      !**NOTE:  Above, sstar and swilt are guesses for all except grassland
      ! and savanna.



      !***************************************************
      !*         END - GISS VEGETATION TYPES             *
      !***************************************************
      end module ent_pfts
