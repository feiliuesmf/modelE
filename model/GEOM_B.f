      MODULE GEOM
!@sum  GEOM contains spherical geometric variables and arrays
!@auth Original development team
!@ver  1.0 (B grid version)
!@cont GEOM_B
      USE CONSTANT, only : OMEGA,RADIUS,TWOPI,SDAY,radian
      USE MODEL_COM, only : IM,JM,LM,FIM,BYIM
!     USE DOMAIN_DECOMP, only : grid, get
!     USE DOMAIN_DECOMP, only : halo_update, north, south, checksum
      IMPLICIT NONE
C**** The primary grid is the A grid (including both poles)
C**** The secondary grid is for the B grid velocities, located on the
C**** vertices of the A grid (Note: first velocity point is J=2)
C**** Polar boxes can have different latitudinal size and are treated
C**** as though they were 1/IM of their actual area
      SAVE
!@param  DLON grid spacing in longitude (deg)
      REAL*8, PARAMETER :: DLON = TWOPI*BYIM
C**** For the wonderland model set DLON=DLON/3
c      REAL*8, PARAMETER :: DLON=TWOPI/(IM*3)
!@param  DLAT,DLAT_DG grid spacing in latitude (rad,deg)
      REAL*8  :: DLAT,DLAT_DG
!@param  FJEQ equatorial value of J
      REAL*8, PARAMETER :: FJEQ=.5*(1+JM)
!@var  J1U index of southernmost latitude (currently 2, later 1)
      INTEGER, parameter :: J1U = 2
!@var  JRANGE_HEMI lowest,highest lat index for SH,NH for A,B grid
      INTEGER, parameter, dimension(2,2,2) :: JRANGE_HEMI = reshape(
     *  (/1,JM/2,  1+JM/2,JM,  J1U,J1U-1+JM/2, J1U-1+JM/2,JM+J1U-2/),
     *  (/2,2,2/))
!@var  LAT latitude of mid point of primary grid box (radians)
      REAL*8, ALLOCATABLE, DIMENSION(:) :: LAT
!@var  LAT_DG latitude of mid points of primary and sec. grid boxs (deg)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: LAT_DG
!@var  LON longitude of mid points of primary grid box (radians)
      REAL*8, ALLOCATABLE, DIMENSION(:) :: LON
!@var  LON_DG longitude of mid points of prim. and sec. grid boxes (deg)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: LON_DG
!@var  DXYP,BYDXYP area of grid box (+inverse) (m^2)
C**** Note that this is not the exact area, but is what is required for
C**** some B-grid conservation quantities
C**** Decomposition exception: DXYP not distributed because individual elements
C**** of the 1D array are needed by all processes.
      REAL*8  DXYP(JM)
      REAL*8, ALLOCATABLE, DIMENSION(:) :: BYDXYP
!@var AREAG global integral of area (m^2)
      REAL*8 :: AREAG
!@var WTJ area weighting used in JLMAP, JKMAP (for hemispheric means)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: WTJ

!@var DXYV,BYDXYV area of grid box around velocity point (recip.)(m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:) :: 
     &                  DXYV,BYDXYV

!@var  DXP,DYP,BYDXP,BYDYP distance between points on primary grid
!@+     (+inverse)
      REAL*8, ALLOCATABLE, DIMENSION(:) :: 
     &                  DXP,DYP,BYDXP,BYDYP
!@var  DXV,DYV distance between velocity points (secondary grid)
      REAL*8, ALLOCATABLE, DIMENSION(:) :: DXV,DYV
!@var  DXYN,DXYS half box areas to the North,South of primary grid point
      REAL*8, ALLOCATABLE, DIMENSION(:) :: DXYS,DXYN
!@var  SINP sin of latitude at primary grid points
      REAL*8, ALLOCATABLE, DIMENSION(:) :: SINP
!@var  COSP, COSV cos of latitude at primary, secondary latitudes
      REAL*8, ALLOCATABLE, DIMENSION(:) :: COSP,COSV
!@var  RAPVS,RAPVN,RAVPS,RAVPN area scalings for primary and sec. grid
      REAL*8, ALLOCATABLE, DIMENSION(:) :: 
     &                  RAPVS,RAPVN,RAVPS,RAVPN

!@var SINIV,COSIV,SINIP,COSIP longitud. sin,cos for wind,pressure grid
      REAL*8, ALLOCATABLE, DIMENSION(:) :: SINIV,COSIV,SINIP,COSIP
!@var  RAVJ scaling for A grid U/V to B grid points (func. of lat. j)
!@var  RAPJ scaling for B grid -> A grid conversion (1/4,1/im at poles)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: 
     &                  RAPJ,RAVJ
!@var  IDJJ J index of adjacent U/V points for A grid (func. of lat. j)
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: IDJJ
!@var  IDIJ I index of adjacent U/V points for A grid (func. of lat/lon)
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: 
     &                  IDIJ
!@var  KMAXJ varying number of adjacent velocity points
      INTEGER, ALLOCATABLE, DIMENSION(:) :: KMAXJ
!@var  IMAXJ varying number of used longitudes
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IMAXJ
!@var  FCOR latitudinally varying coriolis parameter
      REAL*8, ALLOCATABLE, DIMENSION(:) :: FCOR

      real*8 :: acor,acor2,polwt

      CONTAINS

      SUBROUTINE GEOM_B
!@sum  GEOM_B Calculate spherical geometry for B grid
!@auth Original development team (modifications by G. Schmidt)
!@ver  1.0 (B grid version)
      USE DOMAIN_DECOMP, only : grid, get
      USE DOMAIN_DECOMP, only : halo_update, north, south, checksum
      IMPLICIT NONE
      REAL*8, PARAMETER :: EDPERD=1.,EDPERY = 365.

      INTEGER :: I,J,K,IM1  !@var I,J,K,IM1  loop variables
      INTEGER :: JVPO,JMHALF
      REAL*8  :: RAVPO,LAT1,COSP1,DXP1
      INTEGER J_0,J_1,J_0S,J_1S,J_0STG,J_1STG

C**** latitudinal spacing depends on whether you have even spacing or
C**** a partial box at the pole
      DLAT_DG=180./JM                     ! even spacing (default)
      IF (JM.eq.46) DLAT_DG=180./(JM-1)   ! 1/2 box at pole for 4x5
cc    IF (JM.eq.24) DLAT_DG=180./(JM-1)   ! 1/2 box at pole, orig 8x10
      IF (JM.eq.24) DLAT_DG=180./(JM-1.5) ! 1/4 box at pole, 'real' 8x10
      DLAT=DLAT_DG*radian
      CALL GET(grid, J_STRT     =J_0,    J_STOP     =J_1,
     &               J_STRT_SKP =J_0S,   J_STOP_SKP =J_1S,
     &               J_STRT_STGR=J_0STG, j_STOP_STGR=J_1STG)

      IF (grid%HAVE_SOUTH_POLE) THEN
          LAT(1)  = -.25*TWOPI
          SINP(1)  = -1.
          COSP(1)  = 0.
          DXP(1)  = 0.
      END IF
      IF (grid%HAVE_NORTH_POLE) THEN
          LAT(JM) = .25*TWOPI
          SINP(JM) = 1.
          COSP(JM) = 0.
          DXP(JM) = 0.
      END IF
      DO J=J_0S,J_1S
        LAT(J)  = DLAT*(J-FJEQ)
        SINP(J) = SIN(LAT(J))
        COSP(J) = COS(LAT(J))
        DXP(J)  = RADIUS*DLON*COSP(J)
      END DO
      BYDXP(J_0S:J_1S) = 1.D0/DXP(J_0S:J_1S)
C****Update halos for arrays cosp, dxp, and lat
      CALL CHECKSUM(grid, COSP, __LINE__, __FILE__)
      CALL CHECKSUM(grid, DXP , __LINE__, __FILE__)
      CALL CHECKSUM(grid, LAT , __LINE__, __FILE__)
      CALL HALO_UPDATE( grid, COSP, from=SOUTH)
      CALL HALO_UPDATE( grid, DXP , from=SOUTH)
      CALL HALO_UPDATE( grid, LAT , from=SOUTH)

      LAT1    = DLAT*(1.-FJEQ)
      COSP1   = COS(LAT1)
      DXP1    = RADIUS*DLON*COSP1

      DO J=J_0STG,J_1STG
        COSV(J) = .5*(COSP(J-1)+COSP(J))
        DXV(J)  = .5*(DXP(J-1)+DXP(J))
        DYV(J)  = RADIUS*(LAT(J)-LAT(J-1))
C**** The following corrections have no effect for half polar boxes
C**** but are important for full and quarter polar box cases.
        IF (J.eq.2) THEN
          polwt = cosv(j)
          COSV(J) = .5*(COSP1+COSP(J))
          DXV(J)  = .5*(DXP1+DXP(J))
        END IF
        IF (J.eq.JM) THEN
          COSV(J) = .5*(COSP(J-1)+COSP1)
          DXV(J)  = .5*(DXP(J-1)+DXP1)
        END IF
C****
      END DO
      polwt = (cosv(3)-cosv(2))/(cosv(3)-polwt)
C****POLES
      IF (grid%HAVE_SOUTH_POLE) THEN
        DYP(1)  = RADIUS*(LAT(2)-LAT(1)-0.5*DLAT)
        DXYP(1) = .5*DXV(2)*DYP(1)
        BYDXYP(1) = 1./DXYP(1)
        DXYS(1)  = 0.
        DXYN(1)  = DXYP(1)
      ENDIF
      IF (grid%HAVE_NORTH_POLE) THEN
        DYP(JM) = RADIUS*(LAT(JM)-LAT(JM-1)-0.5*DLAT)
        DXYP(JM)= .5*DXV(JM)*DYP(JM)
        BYDXYP(JM) = 1./DXYP(JM)
        DXYS(JM) = DXYP(JM)
        DXYN(JM) = 0.
      END IF
      AREAG = DXYP(1)+DXYP(JM)

      DO J=J_0S,J_1S
        DYP(J)  = radius*dlat !.5*(DYV(J)+DYV(J+1))
        DXYP(J) = .5*(DXV(J)+DXV(J+1))*DYP(J)
        BYDXYP(J) = 1./DXYP(J)
        DXYS(J) = .5*DXYP(J)
        DXYN(J) = .5*DXYP(J)
        AREAG = AREAG+DXYP(J)
      END DO
      BYDYP(:) = 1.D0/DYP(:)
      AREAG = AREAG*FIM
C****EXCEPTION! 
C****      dxyp is not distributed -> make sure all procs. have right values for
C****                                 the whole array. 
C****  CALL ALL_GATHER(DXYP)

C****POLES
      IF (grid%HAVE_SOUTH_POLE) THEN
        RAVPS(1)  = 0.
        RAPVS(1)  = 0.
      END IF
      IF (grid%HAVE_NORTH_POLE) THEN
        RAVPN(JM) = 0.
        RAPVN(JM) = 0.
      END IF
      DO J=J_0STG,J_1STG
        DXYV(J) = DXYN(J-1)+DXYS(J)
        BYDXYV(J) = 1./DXYV(J)
        RAPVS(J)   = .5*DXYS(J)/DXYV(J)
        RAVPS(J)   = .5*DXYS(J)/DXYP(J)
      END DO
      acor = dxyv(2)/(.5*dxp(2)*dyv(2)) ! gridbox area correction factor
      acor2 = dxyv(2)/(dxv(2)*dyv(2))
      CALL CHECKSUM(grid,RAPVN, __LINE__, __FILE__)
      CALL HALO_UPDATE(grid,RAPVN, from=NORTH)
      DO J=J_0,J_1S
        RAPVN(J) = .5*DXYN(J)/DXYV(J+1)
        RAVPN(J) = .5*DXYN(J)/DXYP(J)
      END DO
C**** LONGITUDES (degrees); used in ILMAP
      LON_DG(1,1) = -180.+360./(2.*FLOAT(IM))
      LON_DG(1,2) = -180.+360./    FLOAT(IM)
      DO I=2,IM
        LON_DG(I,1) = LON_DG(I-1,1)+360./FLOAT(IM)
        LON_DG(I,2) = LON_DG(I-1,2)+360./FLOAT(IM)
      END DO
C**** LATITUDES (degrees); used extensively in the diagn. print routines
      IF (grid%HAVE_SOUTH_POLE) LAT_DG(1,1:2)=-90.
      IF (grid%HAVE_NORTH_POLE) LAT_DG(JM,1)=90.
      DO J=J_0S,J_1S
        LAT_DG(J,1)=DLAT_DG*(J-FJEQ)    ! primary (tracer) latitudes
      END DO
      DO J=J_0STG,J_1STG
        LAT_DG(J,2)=DLAT_DG*(J-JM/2-1)  ! secondary (velocity) latitudes
      END DO
C**** WTJ: area weighting for JKMAP, JLMAP hemispheres
      JMHALF= JM/2
      DO J=J_0,J_1
        WTJ(J,1,1)=1.
        WTJ(J,2,1)=2.*FIM*DXYP(J)/AREAG
      END DO
      DO J=J_0STG,J_1STG
        WTJ(J,1,2)=1.
        WTJ(J,2,2)=2.*FIM*DXYV(J)/AREAG
      END DO
C****EQUATOR
      IF (grid%HAVE_EQUATOR) THEN
        WTJ(JMHALF+1,1,2)=.5
        WTJ(JMHALF+1,2,2)=WTJ(JMHALF+1,2,2)/2.
      END IF
C****POLE
      IF (grid%HAVE_SOUTH_POLE) THEN
        WTJ(1,1,2)=0.
        WTJ(1,2,2)=0.
      END IF

C**** CALCULATE CORIOLIS PARAMETER (NOW RESOLUTION INDEPENDENT)
c      OMEGA = TWOPI*(EDPERD+EDPERY)/(EDPERD*EDPERY*SDAY)
      IF (grid%HAVE_SOUTH_POLE)
     &      FCOR(1)  = -OMEGA*DXV(2)*DXV(2)/DLON
      IF (grid%HAVE_NORTH_POLE)
     &      FCOR(JM) = OMEGA*DXV(JM)*DXV(JM)/DLON
C****Update halo for DXV array
      CALL CHECKSUM(grid,DXV, __LINE__, __FILE__)
      CALL HALO_UPDATE(grid,DXV, from=NORTH)
      DO J=J_0S,J_1S
        FCOR(J) = OMEGA*(DXV(J)*DXV(J)-DXV(J+1)*DXV(J+1))/DLON
      END DO

C**** Set indexes and scalings for the influence of A grid points on
C**** adjacent velocity points

C**** Calculate relative directions of polar box to nearby U,V points
      DO I=1,IM
        SINIV(I)=SIN((I-1)*DLON)
        COSIV(I)=COS((I-1)*TWOPI*BYIM) ! DLON)
        LON(I)=DLON*(I-.5)
        SINIP(I)=SIN(LON(I))
        COSIP(I)=COS(LON(I))
      END DO

C**** Conditions at the poles
      DO J=1,JM,JM-1
        IF (((J .EQ. 1) .AND. (grid%HAVE_SOUTH_POLE)) .OR.
     &      ((J .EQ.JM) .AND. (grid%HAVE_NORTH_POLE))) THEN
          IF(J.EQ.1) THEN
            JVPO=2
            RAVPO=2.*RAPVN(1)
          ELSE
            JVPO=JM
            RAVPO=2.*RAPVS(JM)
          END IF
          KMAXJ(J)=IM
          IMAXJ(J)=1
          RAVJ(1:KMAXJ(J),J)=RAVPO
          RAPJ(1:KMAXJ(J),J)=BYIM
          IDJJ(1:KMAXJ(J),J)=JVPO
          DO K=1,KMAXJ(J)
          IDIJ(K,1:IM,J)=K
        END DO
        END IF
      END DO
C**** Conditions at non-polar points
      DO J=J_0S,J_1S
        KMAXJ(J)=4
        IMAXJ(J)=IM
        DO K=1,2
          RAVJ(K,J)=RAPVS(J)
          RAPJ(K,J)=RAVPS(J)    ! = .25
          IDJJ(K,J)=J
          RAVJ(K+2,J)=RAPVN(J)
          RAPJ(K+2,J)=RAVPN(J)  ! = .25
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

      SUBROUTINE ALLOC_GEOM(grd_dum)
      USE MODEL_COM, only : IM
      USE DOMAIN_DECOMP, only: DYN_GRID, get
      USE GEOM, only: LAT,LAT_DG,
     &                LON,LON_DG,
     &                DXYP,BYDXYP,WTJ,
     &                DXYV,BYDXYV,
     &                DXP,DYP,BYDXP,BYDYP,
     &                DXV,DYV,DXYS,DXYN,
     &                SINP,COSP,COSV,
     &                RAPVS,RAPVN,RAVPS,RAVPN,
     &                SINIV,COSIV,SINIP,COSIP,
     &                RAPJ,RAVJ,
     &                IDJJ,IDIJ,KMAXJ,IMAXJ,
     &                FCOR

      TYPE (DYN_GRID),   INTENT(IN)    :: grd_dum

      CALL GET(grd_dum,J_STRT_HALO=J_0H,
     &                 J_STOP_HALO=J_1H)

      ALLOCATE( LAT       (J_0H:J_1H), 
     &          LAT_DG    (J_0H:J_1H, 2),
     &          LON    (IM),
     &          LON_DG (IM,           2),
     &                                    STAT = IER)
      ALLOCATE(
     &          BYDXYP    (J_0H:J_1H),
     &          WTJ       (J_0H:J_1H, 2, 2),
     &          DXYV      (J_0H:J_1H),
     &          BYDXYV    (J_0H:J_1H),
     &          DXP       (J_0H:J_1H),
     &          DYP       (J_0H:J_1H),
     &          BYDXP     (J_0H:J_1H),
     &          BYDYP     (J_0H:J_1H),
     &          DXV       (J_0H:J_1H),
     &          DYV       (J_0H:J_1H),
     &          DXYS      (J_0H:J_1H),
     &          DXYN      (J_0H:J_1H),
     &          SINP      (J_0H:J_1H),
     &          COSP      (J_0H:J_1H),
     &          COSV      (J_0H:J_1H),
     &          RAPVS     (J_0H:J_1H),
     &          RAPVN     (J_0H:J_1H),
     &          RAVPS     (J_0H:J_1H),
     &          RAVPN     (J_0H:J_1H),
     &                                    STAT = IER)
      ALLOCATE(
     &          SINIV  (IM),
     &          COSIV  (IM),
     &          SINIP  (IM),
     &          COSIP  (IM),
     &          RAPJ   (IM,J_0H:J_1H),
     &          RAVJ   (IM,J_0H:J_1H),
     &          IDJJ   (IM,J_0H:J_1H),
     &          IDIJ(IM,IM,J_0H:J_1H),
     &          KMAXJ     (J_0H:J_1H),
     &          IMAXJ     (J_0H:J_1H),
     &          FCOR      (J_0H:J_1H),
     &                                    STAT = IER)

      RETURN
      END SUBROUTINE ALLOC_GEOM
