      MODULE PBLCOM

      USE E001M12_COM, only : im,jm
      USE SOCPBL, only : n

      IMPLICIT NONE

!@var uabl boundary layer profile for zonal wind
!@var vabl boundary layer profile for meridional wind
!@var tabl boundary layer profile for temperature
!@var qabl boundary layer profile for humidity
      real*8, dimension(n,im,jm,4) :: uabl,vabl,tabl,qabl
!@var eabl boundary layer profile for turbulent KE (calc. on sec. grid)
      real*8, dimension(n,im,jm,4) :: eabl

!@var cmgs drag coefficient (dimensionless surface momentum flux)
!@var chgs Stanton number   (dimensionless surface heat flux)
!@var cqgs Dalton number    (dimensionless surface moisture flux)
      real*8, dimension(im,jm,4) :: cmgs,chgs,cqgs

!@var ipbl flag for whether pbl properties were found at last timestep
      integer, dimension(im,jm,4) :: ipbl

C****
C**** ROUGHL    LOG10(ZGS/ROUGHNESS LENGTH), prescribed with ZGS=30 m.
C****
      double precision, dimension(im,jm) :: roughl

C**** BLDATA 1  COMPOSITE SURFACE WIND MAGNITUDE (M/S)
C****        2  COMPOSITE SURFACE AIR TEMPERATURE (K)
C****        3  COMPOSITE SURFACE AIR SPECIFIC HUMIDITY (1)
C****        4  LAYER TO WHICH DRY CONVECTION MIXES (1)
C****        5  MIXED LAYER DEPTH (Z1O NOT YET PART OF RESTART FILE)
C****        6  COMPOSITE SURFACE U WIND
C****        7  COMPOSITE SURFACE V WIND
C****        8  COMPOSITE SURFACE MOMENTUM TRANSFER (TAU)
C new     9-12  ustar for each ITYPE (sqrt of srfc momentum flux) (m/s)

      double precision, dimension(im,jm) ::
     &     wsavg,tsavg,qsavg,dclev,mld,usavg,vsavg,tauavg
      double precision, dimension(im,jm,4) :: ustar

      common /bleq/
     &     wsavg,tsavg,qsavg,dclev,mld,usavg,vsavg,tauavg,ustar
      DOUBLE PRECISION, DIMENSION(IM,JM,12) :: BLDATA
      EQUIVALENCE(BLDATA(1,1,1),WSAVG(1,1))

C**** model related constants (should really be taken from E001M12_COM)
      integer, parameter ::  iq1=im/4+1,iq2=im/2+1,iq3=3*im/4+1

C**** pressure gradient arrays
      double precision, dimension(im,jm) ::
     &     dpdxr,dpdyr,phi,dpdxr0,dpdyr0

      END MODULE PBLCOM
