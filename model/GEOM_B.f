      MODULE GEOM
!@sum  GEOM contains spherical geometric variables and arrays
!@auth Original development team
!@ver  1.0 (B grid version)
!@cont GEOM_B
c      USE CONSTANT, only : OMEGA,RADIUS,TWOPI
      USE E001M12_COM, only : IM,JM,FIM,RADIUS,TWOPI,OMEGA,AREAG
      IMPLICIT NONE
C**** The primary grid is the A grid (including both poles)
C**** The secondary grid is for the B grid velocities, located on the
C**** vertices of the A grid
C**** Polar boxes have different latitudinal size and are treated 
C**** as though they were 1/IM of their actual area

!@var  DLON grid spacing in longitude (deg)
c      REAL*8, PARAMETER :: DLON=TWOPI/IM
C**** For the wonderland model set DLON=DLON/3
c      REAL*8, PARAMETER :: DLON=TWOPI/(IM*3)
!@var  DLAT grid spacing in latitude (deg)
c      REAL*8, PARAMETER :: DLAT=.5*TWOPI/(JM-1)
!@var  FJEQ equatorial value of J 
      REAL*8, PARAMETER :: FJEQ=.5*(1+JM)

!@var  LAT latitude of mid point of primary grid box (radians)
      REAL*8, DIMENSION(JM) :: LAT
!@var  DXYP,BYDXYP area of grid box (+inverse) (m^2)
C**** Note that this is not the exact area, but is what is required for
C**** some B-grid conseravtion quantities
      REAL*8, DIMENSION(JM) :: DXYP,BYDXYP
!@var  AREAG global integral of area (m^2)
c      REAL*8 :: AREAG

!@var  DXYV area of grid box around a velocity point (m^2)
      REAL*8, DIMENSION(JM) :: DXYV

!@var  DXP,DYP distance between points on primary grid  
      REAL*8, DIMENSION(JM) :: DXP,DYP
!@var  DXP,DYP distance between velocity points (secondary grid)
      REAL*8, DIMENSION(JM) :: DXV,DYV
!@var  DXYN,DXYS half box areas to the North and South of primary grid point
      REAL*8, DIMENSION(JM) :: DXYS,DXYN
!@var  SINP sin of latitude at primary grid points 
      REAL*8, DIMENSION(JM) :: SINP
!@var  COSP, COSV cos of latitude at primary, secondary latitudes 
      REAL*8, DIMENSION(JM) :: COSP,COSV
!@var  RAPVS,RAPVN,RAVPS,RAVPN area scalings for primary and secondary grid
      REAL*8, DIMENSION(JM) :: RAPVS,RAPVN,RAVPS,RAVPN

!@var  RAJ scaling to get A grid velocities to B grid points
!      as function of latitude j
      REAL*8, DIMENSION(IM,JM) :: RAJ 
!@var  IDJJ J index of adjacent velocity points for A grid points
!      as function of latitude j
      INTEGER, DIMENSION(IM,JM) :: IDJJ
!@var  IDIJ I index of adjacent velocity points for A grid points
!      as function of longitude i and latitude j
      INTEGER, DIMENSION(IM,IM,JM) :: IDIJ
!@var  KMAXJ varying number of adjacent velocity points
      INTEGER, DIMENSION(JM) :: KMAXJ
!@var  IMAXJ varying number of used longitudes 
      INTEGER, DIMENSION(JM) :: IMAXJ

!@var  FCOR latitudinally varying coriolis parameter
      REAL*8, DIMENSION(JM) :: FCOR

      CONTAINS

      SUBROUTINE GEOM_B
!@sum  GEOM_B Calculate spherical geometry for B grid
!@auth Original development team (modifications by G. Schmidt)
!@ver  1.0 (B grid version)
      USE E001M12_COM
      IMPLICIT NONE
      REAL*8, PARAMETER :: EDPERD=1.,EDPERY = 365.

      INTEGER :: I,J,K,IM1  !@var I,J,K,IM1  loop variables
      INTEGER :: JVPO
      REAL*8  :: RAPO

      DLON=TWOPI/IM
      DLAT=.5*TWOPI/(JM-1)
      LAT(1)  = -.25*TWOPI
      LAT(JM) = -LAT(1)
      SINP(1)  = -1.
      SINP(JM) = 1.
      COSP(1)  = 0.
      COSP(JM) = 0.
      DXP(1)  = 0.
      DXP(JM) = 0.
      DO J=2,JM-1
         LAT(J)  = DLAT*(J-FJEQ)
         SINP(J) = SIN(LAT(J))
         COSP(J) = COS(LAT(J))
         DXP(J)  = RADIUS*DLON*COSP(J)
      END DO
      DO J=2,JM
         COSV(J) = .5*(COSP(J-1)+COSP(J))
         DXV(J)  = .5*(DXP(J-1)+DXP(J))
         DYV(J)  = RADIUS*(LAT(J)-LAT(J-1))
      END DO
      DYP(1)  = .5*DYV(2)
      DYP(JM) = .5*DYV(JM)
      DXYP(1) = .5*DXV(2)*DYP(1)
      BYDXYP(1) = 1./DXYP(1)
      DXYP(JM)= .5*DXV(JM)*DYP(JM)
      BYDXYP(JM) = 1./DXYP(JM)
      DXYS(1)  = 0.
      DXYS(JM) = DXYP(JM)
      DXYN(1)  = DXYP(1)
      DXYN(JM) = 0.
      AREAG = DXYP(1)+DXYP(JM)
      DO J=2,JM-1
         DYP(J)  = .5*(DYV(J)+DYV(J+1))
         DXYP(J) = .5*(DXV(J)+DXV(J+1))*DYP(J)
         BYDXYP(J) = 1./DXYP(J)
         DXYS(J) = .5*DXYP(J)
         DXYN(J) = .5*DXYP(J)
         AREAG = AREAG+DXYP(J)
      END DO
      AREAG = AREAG*FIM
      RAVPS(1)  = 0.
      RAVPN(JM) = 0.
      DO J=2,JM
         DXYV(J) = DXYN(J-1)+DXYS(J)
         RAPVS(J)   = .5*DXYS(J)/DXYV(J)
         RAPVN(J-1) = .5*DXYN(J-1)/DXYV(J)
         RAVPS(J)   = .5*DXYS(J)/DXYP(J)
         RAVPN(J-1) = .5*DXYN(J-1)/DXYP(J-1)
      END DO
C**** CALCULATE CORIOLIS PARAMETER
      OMEGA = TWOPI*(EDPERD+EDPERY)/(EDPERD*EDPERY*SDAY)
      FCOR(1)  = -RADIUS*OMEGA*.5*COSP(2)*DXV(2)
      FCOR(JM) = -FCOR(1)
      DO J=2,JM-1
         FCOR(J) = OMEGA*(DXV(J)*DXV(J)-DXV(J+1)*DXV(J+1))/DLON
      END DO

C**** Set indexes and scalings for the influence of A grid points on
C**** adjacent velocity points (assumes that U,V are sequential in a
C**** common block) 

C**** Conditions at the poles
      DO J=1,JM,JM-1
         IF(J.EQ.1) THEN
            JVPO=2
            RAPO=2.*RAPVN(1)
         ELSE
            JVPO=JM
            RAPO=2.*RAPVS(JM)
         END IF
         KMAXJ(J)=IM
         IMAXJ(J)=1
         RAJ(1:KMAXJ(J),J)=RAPO
         IDJJ(1:KMAXJ(J),J)=JVPO
         DO I=1,IM
            IDIJ(I,1:IM,J)=I
         END DO
      END DO
C**** Conditions at non-polar points
      DO J=2,JM-1
         KMAXJ(J)=4
         IMAXJ(J)=IM
         DO K=1,2
            RAJ(K,J)=RAPVS(J)
            IDJJ(K,J)=J
            RAJ(K+2,J)=RAPVN(J)
            IDJJ(K+2,J)=J+1
         END DO
         IM1=IM
         DO I=1,IM
            IDIJ(1,I,J)=IM1
            IDIJ(2,I,J)=I
            IDIJ(3,I,J)=IM1
            IDIJ(4,I,J)=I
            IM1=I
         END DO
      END DO

      RETURN
      END SUBROUTINE GEOM_B

      END MODULE GEOM

