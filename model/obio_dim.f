      MODULE obio_dim

!  Biological constituents and their units
!  P(1) = nitrate (uM or micro-moles/lt=mili-moles/m^3)
!  P(2) = ammonium (uM)
!  P(3) = silica (uM)
!  P(4) = iron (nM)
!  P(5) = diatoms (mg chl m-3)
!  P(6) = chlorophytes (mg chl m-3)
!  P(7) = cyanobacteria (mg chl m-3)
!  P(8) = coccolithophores (mg chl m-3)
!  P(9) = herbivores (mg chl m-3)
!  P(10)= inert tracer modelled after nitrate
!  Detrital components
!  det(1) = N/C detritus (uM)
!  det(2) = silica detritus (uM)
!  det(3) = iron detritus (nM)
!  Carbon components
!  C(1) = DOC (uM)
!  C(2) = DIC (uM)
!  pCO2 (uatm)

      implicit none

      SAVE

      integer, parameter :: nnut=4
     .                     ,nchl=4
     .                     ,nzoo=1
     .                     ,ntyp=nnut+nchl+nzoo
     .                     ,n_inert=1
     .                     ,ndet=3
     .                     ,ncar=2
     .                     ,FLUXDEF=0


      integer, parameter :: nlt=33,   !number of spectral channels
     .                      nh=200,   !number of depths for mean irradiance
     .                      nch=48,   !number of chl values for mean irrad
     .                      ncd=41    !number of cdom values for mean irrad

      integer, parameter :: npr=15,   !number of spectral values in par
     .                      nhn=12,   !number hourly oasim values per day
     .                      npar=npr  !same as npr


      integer, parameter :: nrg=13    !number of oceanographic basins

      END MODULE obio_dim
