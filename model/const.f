      MODULE CONSTANT
!@sum  CONSTANT definitions for physical constants and useful numbers
!@auth G. Schmidt
!@ver  1.0

C**** Note that there are often two definitions used here: the one that
C**** is current and the 'standard' definition to which we will move.
C**** Note also some non-double precision numbers, temporarily for
C**** consistency 

C**** Conventions: 'by' implies reciprocal, 'rt' implies square root

C**** Numerical constants

      real*8,parameter :: pi = 3.1415926535897932d0 !@param pi    pi
      real*8,parameter :: twopi = 2d0*pi           !@param twopi 2*pi
      real*8,parameter :: radian = pi/180d0        !@param radian pi/180
!@param rt2,byrt2   sqrt(2), 1/sqrt(2)
      real*8,parameter :: rt2 = 1.4142135623730950d0
      real*8,parameter :: byrt2 = 0.70710678118654752d0
!@param rt3,byrt3   sqrt(3), 1/sqrt(3)
      real*8,parameter :: rt3 = 1.7320508075688772d0
      real*8,parameter :: byrt3 = 0.57735026918962576d0
!@param rt12,byrt12   sqrt(12), 1/sqrt(12)
      real*8,parameter :: rt12 = 3.4641016151377546d0
      real*8,parameter :: byrt12 = 0.28867513459481288d0

      real*8,parameter :: by3 =0.33333333333333333d0  !@param by3  1/3
      real*8,parameter :: by9 =0.11111111111111111d0  !@param by9  1/9
      real*8,parameter :: by12=0.83333333333333333d-1 !@param by12 1/12

C**** Physical constants

!@param sday  sec per day (s)
      real*8,parameter :: sday = 86400.

!@param stbo Stefan-Boltzmann constant (W/m^2 K^4) 
c      real*8,parameter :: stbo =5.67051d-8 !current best estimate
      real*8,parameter :: stbo = 5.67032E-8

!@param lhe   latent heat of evap (J/kg)
c**** There should be a temperature dependence (Gill):
c**** lhe = 2.5008d6 - 2.3d3 T (in C) 
      real*8,parameter :: lhe = 2.5d6  
!@param lhm   latent heat of melt (J/kg)
c**** There should be a temperature dependence (Gill):
c**** lhe = 334590 + 2.05d3 T (in C) 
      real*8,parameter :: lhm = 3.34d5
!@param lhs  latent heat of sublimation (J/kg)
      real*8,parameter :: lhs = lhe+lhm

!@param rhow density of pure water (kg/m^3)
      real*8,parameter :: rhow = 1d3
!@param rhoi density of pure ice (kg/m^3)
      real*8,parameter :: rhoi = 916.6   !d0

!@param tf freezing point of water at 1 atm (K)
      real*8,parameter :: tf = 273.16   !d0
!@param bytf 1/tf (K^-1)
      real*8,parameter :: bytf = 1d0/tf

!@param shw heat capacity of water (at 20 C) (J/kg C)
      real*8,parameter :: shw  = 4185.
!@param byshw 1/shw 
      real*8,parameter :: byshw = 1d0/shw

!@param shi heat capacity of pure ice (at 0 C) (J/kg C)
      real*8,parameter :: shi  = 2060.
!@param byshi 1/shi 
      real*8,parameter :: byshi = 1d0/shi

!@param rgas gas constant (J/K kg)
c**** for values of CO2 much larger than present day (4x conc)
c**** the molecular weight of dry air M_A could change
c**** RGAS = R/M_A = 1000* 8.314510 J/mol K /28.9655 g/mol 
c**** (US Stand. Atm.) assuming M_w for O2 = 31.9988 and CO2 = 44.00995
c**** and current percentage 20.946% and 0.0350% 
c**** assuming CO2 displaces other gases equally M_A=28.9602 + n*0.00527
c**** where n is multiple of present day CO2 conc (350 ppm)
c**** For 4xCO2  M_A = 28.9813  => rgas = 286.89
c**** For 10xCO2 M_A = 29.0129  => rgas = 286.58
c      real*8,parameter :: rgas =287.05d0  
      real*8,parameter :: rgas = 287.

!@param rvap  gas constant for water vapour (J/K kg)
c**** defined as R/M_W = 1000* 8.314510 J/mol K /18.015 g/mol
c      real*8,parameter :: rvap =461.53d0
      real*8,parameter :: rvap = 461.5

!@param kapa ideal gas law exponent  (1)
c**** constant = (g-1)/g where g=1401 = c_p/c_v 
c      real*8,parameter :: kapa =.2862d0   ! best value
      real*8,parameter :: kapa = .286     !d0
!@param bykapa,bykapap1,bykapap2 various useful reciprocals of kapa
      real*8,parameter :: bykapa = 1./kapa
      real*8,parameter :: bykapap1 = 1./(kapa+1.)
      real*8,parameter :: bykapap2 = 1./(kapa+2.)
 
!@param sha specific heat of dry air (J/kg C)
      real*8,parameter :: sha = rgas/kapa
!@param bysha 1/sha 
      real*8,parameter :: bysha = 1./sha

!@param shv specific heat of vapour (J/kg C)
c**** shv is currently assumed to be zero to aid energy conservation in
c**** the atmosphere. Once the heat content associated with water
c**** vapour is included, this can be set to the standard value
c      real*8,parameter :: shv = sha/1401. ! = c_p/g
      real*8,parameter :: shv = 0.

C**** Astronomical constants

!@param omega earth's rotation rate (s^-1)
c      real*8,parameter :: omega = 7.2921151467d-5 ! NOVAS value
      real*8,parameter :: EDPERD = 1.
      real*8,parameter :: EDPERY = 365.
      real*8,parameter :: omega = TWOPI*(EDPERD+EDPERY)/
     *                            (EDPERD*EDPERY*SDAY)

!@param radius radius of the earth (m)
c**** radius of spherical earth, same volume = 6371000 m
      real*8,parameter :: radius = 6375000.
!@param grav gravitaional accelaration (m/s^2)
c**** SI reference gravity (at 45 deg) = 9.80665
      real*8,parameter :: grav = 9.81    !d0
!@param bygrav 1/grav 
      real*8,parameter :: bygrav = 1d0/grav

      END MODULE CONSTANT
