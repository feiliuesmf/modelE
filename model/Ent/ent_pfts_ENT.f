      module ent_pfts
!@sum Ent default supported 13 vegetation types

      !use ent_pftconst
      use ent_const
      use ent_types
      implicit none

      !***************************************************
      !*      ENT PLANT FUNCTIONAL TYPES                 *
      !***************************************************

      !* 1 - evergreen broadleaf early successional
      !* 2 - evergreen broadleaf late successional
      !* 3 - evergreen needleleaf early successional
      !* 4 - evergreen needleleaf late successional
      !* 5 - cold deciduous broadleaf early successional
      !* 6 - cold deciduous broadleaf late successional
      !* 7 - drought deciduous broadleaf
      !* 8 - cold adapted shrub
      !* 9 - arid adapted shrub
      !* 10- C3 grass
      !* 11- C4 grass
      !* 12- arctic C3 grass
      !* 13- C4 crops

!KIM-temp.
      integer, parameter :: GRASSC3 = 10
      integer, parameter :: CROPS = 13

      !pst - photosynthetic pathway, 1-C3, 2-C4
      !hwilt - wilting point (m)
      !sstar - soil moisture stress onset point (fraction of soil volumetric saturation)
      !swilt - wilting point (fraction of soil volumetric saturation)
      !nf - canopy nitrogen factor (dimensionless)
      !* Parameters from CASA:
      !sla - specific leaf area (m2 leaf area/kg leaf C)(CASA)
      !lrage - leaf and root litter age (years) (CASA)
      !woodage - stem litter age (years) (CASA)
      !lit_C2N - litter C:N (CASA)
      !lignin - lignin content

      !SOURCES:
      !soil  moisture points:  7 & 10 from Kiang (2002) dissertation.
      !                       Other values are guesses.
      !nf:  Kull&Kruijt ps cap/leaf N param. Guesses from Friend&Kiang (2005) 
      !sla:  from LSM parameterizations.
      !lrage, woodage,lit_C2N,lignin: from CASA parameterizations.
      
      !***************************************************
      !Temp values for Ent pfts (See ent_const.f for types)
      type(pftype),parameter :: pfpar(N_PFT) =         !PFT parameters
      !     pst,  hwilt, sstar, swilt,nf,sla,r,lrage,woodage,lit_C2N,lignin !
     &           (/      
     &     pftype(1, -100.d0, .60d0, .29d0, 1.3d0, 25.0d0, 0.85d0, 
     &     1.8d0, 41.0d0, 40.d0, 0.2d0),
     &     pftype(1, -100.d0, .60d0, .29d0, 1.2d0, 25.0d0, 0.85d0, 
     &     1.8d0,41.0d0, 40.d0, 0.2d0),
     &     pftype(1, -100.d0, .55d0, .25d0, 0.9d0, 10.0d0, 0.9d0, 
     &     5.0d0, 42.0d0, 80.d0, 0.25d0),
     &     pftype(1, -100.d0, .55d0, .25d0, 0.85d0,10.0d0, 0.9d0, 
     &     5.0d0, 42.0d0, 80.d0, 0.25d0),
     &     pftype(1, -100.d0, .50d0, .29d0, 1.5d0, 30.0d0, 0.9d0, 
     &     1.2d0, 58.0d0, 50.d0, 0.2d0),
     &     pftype(1, -100.d0, .50d0, .29d0, 1.4d0, 30.0d0, 0.9d0, 
     &     1.2d0, 58.0d0, 50.d0, 0.2d0),
     &     pftype(1, -100.d0, .45d0, .22d0, 1.4d0, 30.0d0, 0.85d0, 
     &     1.2d0,25.0d0, 60.d0, 0.2d0),
     &     pftype(1, -100.d0, .50d0, .22d0, 1.4d0, 25.0d0, 0.6d0, 
     &     2.8d0, 5.5d0, 50.d0, 0.15d0),
     &     pftype(1, -100.d0, .40d0, .22d0, 1.3d0, 25.0d0, 0.85d0, 
     &     1.0d0, 5.5d0, 65.d0, 0.2d0),
     &     pftype(1, -100.d0, .65d0, .27d0, 1.5d0, 40.0d0, 0.6d0, 
     &     1.5d0, 0.0d0, 50.d0, 0.1d0),
     &     pftype(2, -100.d0, .55d0, .22d0, 0.76d0,35.0d0, 0.6d0, 
     &     1.5d0, 0.0d0, 50.d0, 0.1d0),
     &     pftype(1, -100.d0, .60d0, .27d0, 1.4d0, 20.0d0, 0.6d0, 
     &     1.5d0, 0.0d0, 50.d0, 0.1d0),
     &     pftype(2, -100.d0, .65d0, .27d0, 0.76d0,60.0d0, 1.2d0, 
     &     1.2d0, 1.0d0, 40.d0, 0.1d0)
     &     /)

      ! other parameters needed for Ent to compile
      integer, parameter :: COVEROFFSET = 0

      end module ent_pfts
