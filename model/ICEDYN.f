#include "rundeck_opts.h"

      MODULE ICEDYN
!@sum  ICEDYN holds variables for dynamic sea ice
!@auth Gary Russell/Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : im,jm,dts=>dtsrc
      USE CONSTANT, only : rhoi,omega,radian
      IMPLICIT NONE
      SAVE

C**** variables used in advection
!@var RSIX,RSIY first order moments for seaice concentration
!@var USI,VSI east-west, and north-south sea ice velocities (m/s)
      REAL*8, DIMENSION(IM,JM) :: RSIX,RSIY,USI,VSI
!@var RSISAVE saved value of sea ice concentration before DYNSI
      REAL*8, DIMENSION(IM,JM) :: RSISAVE

C**** local grid variables for ice rheology scheme
!@var nx1 number of grid points in the longitudinal direction
!@+  (calculated points from 2 through nx1-1. End points are boundaries)
!@var ny1 number of grid points in the latitudinal direction
!@+  (calculated points from 2 through ny1-1. End points are boundaries)
      integer, parameter :: nx1=im+2, ny1=jm

C**** ice dynamics input/output/work variables
      REAL*8, DIMENSION(NX1,NY1) :: PRESS,HEFFM,UVM,DWATN,COR,SINEN
     *     ,BYDXDY
      REAl*8, DIMENSION(NX1) :: DXT,DXU,BYDX2,BYDXR
      REAl*8, DIMENSION(NY1) :: DYT,DYU,BYDY2,BYDYR,CST,CSU,TNGT,TNG
     *     ,SINE
!@var RATIC,RICAT area ratios for atmospheric tracer grid to ice grid
      REAL*8, DIMENSION(NY1) :: RATIC, RICAT

!@var OIPHI ice-ocean turning angle (25 degrees)
!@var ECCEN value of eccentricity for ellipse
      REAL*8, PARAMETER :: ECCEN=2.0, OIPHI=25d0*radian
!@var SINWAT,COSWAT sin and cos of ice-ocean turning angle
      REAL*8 SINWAT,COSWAT

!@var USIDT,VSIDT sea ice mass fluxes, saved for advection
      REAL*8, DIMENSION(IM,JM) :: USIDT,VSIDT

!@var PSTAR maximum sea ice pressure
      REAL*8, PARAMETER :: PSTAR=2.75d4

!@var BYDTS reciprocal of timestep in ice dynamics code
      REAl*8 :: BYDTS

      END MODULE ICEDYN

      SUBROUTINE io_icedyn(kunit,iaction,ioerr)
!@sum  io_icedyn reads and writes dynamic ice arrays to file
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : ioread,iowrite,irsfic,irerun,lhead
      USE ICEDYN
      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var HEADER Character string label for individual records
      CHARACTER*80 :: HEADER, MODULE_HEADER = "ICEDYN01"

      MODULE_HEADER(lhead+1:80) = 'R8 dim(im,jm):RSIX,RSIY,USI,VSI'

      SELECT CASE (IACTION)
      CASE (:IOWRITE)            ! output to standard restart file
        WRITE (kunit,err=10) MODULE_HEADER,RSIX,RSIY,USI,VSI
      CASE (IOREAD:)            ! input from restart file
        SELECT CASE (IACTION)
        CASE (IRSFIC)           ! initial conditions
        CASE (ioread,irerun)    ! restarts
          READ (kunit,err=10) HEADER,RSIX,RSIY,USI,VSI
          IF (HEADER(1:LHEAD).NE.MODULE_HEADER(1:LHEAD)) THEN
            PRINT*,"Discrepancy in module version ",HEADER,MODULE_HEADER
            GO TO 10
          END IF
        END SELECT
      END SELECT

      RETURN
 10   IOERR=1
      RETURN
C****
      END SUBROUTINE io_icedyn

      SUBROUTINE DYNSI
!@sum  DYNSI calculate ice velocites
!@auth Jiping Liu/Gavin Schmidt (based on code from J. Zhang)
!@ver  1.0
      USE CONSTANT, only : radian,rhoi,radius,grav
      USE MODEL_COM, only : p,ptop,itime
C**** Dynamic sea ice should be on the ocean grid
      USE OCEAN, only : im,jm,focean,lmu,lmv,uo,vo,dxyno,dxyso,dxyvo
     *     ,dxypo,bydxypo,ratoc,rocat,dxpo,dyvo,opress,ogeoz,imaxj
      USE ICEDYN, only : rsix,rsiy,usi,vsi,nx1,ny1,press,heffm,uvm
     *     ,dwatn,cor,sinen,dxt,dxu,dyt,dyu,cst,csu,tngt,tng,sine,usidt
     *     ,vsidt,bydx2,bydxr,bydxdy,bydy2,bydyr,dts,sinwat,coswat,oiphi
     *     ,ratic,ricat,bydts,rsisave,omega
      USE FLUXES, only : dmua,dmva,dmui,dmvi,UI2rho
      USE SEAICE, only : ace1i
      USE SEAICE_COM, only : rsi,msi,snowi
      USE ODIAG, only : oij,ij_usi,ij_vsi,ij_dmui,ij_dmvi,ij_pice !,ij_rsi
      IMPLICIT NONE
      SAVE
C****
      REAL*8, PARAMETER :: BYRHOI=1D0/RHOI
C**** working arrays for V-P rheology scheme (on B grid)
      REAL*8, DIMENSION(NX1,NY1,3) :: UICE,VICE
      REAL*8, DIMENSION(NX1,NY1,3) :: HEFF,AREA
      REAL*8, DIMENSION(NX1,NY1) ::  ETA,ZETA,DRAGS,DRAGA,GAIRX,GAIRY
     *     ,GWATX,GWATY,PGFUB,PGFVB,FORCEX,FORCEY,AMASS,UICEC,VICEC,UIB
     *     ,VIB,DMU,DMV,USAVE,VSAVE
C**** intermediate calculation for pressure gradient term
      REAL*8, DIMENSION(IM,JM) :: PGFU,PGFV

C**** should be defined based on ocean grid
      REAL*8 :: dlat,dlon,phit,phiu,hemi,rms
      INTEGER, SAVE :: IFIRST = 1
      INTEGER I,J,n,k,kki,ip1,im1,sumk,l
      REAL*8 USINP,DMUINP,RAT,duA,dvA

C**** set up initial parameters
      IF(IFIRST.gt.0) THEN
      IFIRST = 0
c****
      dlat=nint(180./(jm-1))*radian
      dlon=nint(360./im)*radian
      bydts = 1./dts
c****
      do j = 1,ny1
        dyt(j) = dlat*radius
        dyu(j) = dlat*radius
      enddo
      do i=1,nx1
        dxt(i) = dlon*radius
        dxu(i) = dlon*radius
      enddo
      dxt(1) = dxt(nx1-1)
      dxt(nx1) = dxt(2)
      dxu(1) = dxu(nx1-1)
      dxu(nx1) = dxu(2)

      do i=1,nx1
        bydx2(i)=0.5/(dxu(i)*dxu(i))
        bydxr(i)=0.5/(dxu(i)*radius)
      end do
      do j=1,ny1
        bydy2(j)=0.5/(dyu(j)*dyu(j))
        bydyr(j)=0.5/(dyu(j)*radius)
        do i=1,nx1
          bydxdy(i,j) = 0.5/(dxu(i)*dyu(j))
        end do
      end do

      do j = 1,ny1
       phit = (-90.+(j-1)*4.)*radian
       phiu = (-88.+(j-1)*4.)*radian
       cst(j) = cos(phit)
       csu(j) = cos(phiu)
       sine(j) = sin(phiu)
       tng(j) = sine(j)/csu(j)
       TNGT(J)=SIN(PHIT)/CST(J)
       DO I=1,NX1
         SINEN(I,J)=SIN(PHIU)
       enddo
      enddo
      TNGT(NY1)=TNGT(NY1-1)
      TNG(NY1)=TNG(NY1-1)
      CSU(NY1)=CSU(NY1-1)
C**** set area ratios for converting fluxes
      RATIC = RATOC  ! currently the same as for ocean
      RICAT = ROCAT

C**** Set land masks for tracer and velocity points
C**** This is defined based on the ocean grid
      do j=1,ny1
        do i=2,nx1-1
          heffm(i,j)=nint(focean(i-1,j))
        enddo
        heffm(1,j)=heffm(nx1-1,j)
        heffm(nx1,j)=heffm(2,j)  
      enddo
C**** define velocity points (including exterior corners)
      do j=1,ny1-1
        do i=1,nx1-1
c          sumk=heffm(i,j)+heffm(i+1,j)+heffm(i,j+1)+heffm(i+1,j+1)
c          if (sumk.ge.3) uvm(i,j)=1  ! includes exterior corners
          uvm(i,j) = nint(min(heffm(i,j), heffm(i+1,j), heffm(i,j+1),
     *         heffm(i+1,j+1)))
        end do
      end do
C**** reset tracer points to surround velocity points (except for single
      do j=2,ny1-1
        do i=2,nx1-1
          k = nint(max (uvm(i,j), uvm(i-1,j), uvm(i,j-1), uvm(i-1,j-1)))
c         sumk = nint(uvm(i,j)+uvm(i+1,j)+uvm(i,j+1)+uvm(i+1,j+1))
c set to k except if an island 
c         if (.not. (sumk.eq.4.and.focean(i-1,j).eq.0) ) then
            heffm(i,j) = k
c         end if
        enddo
      enddo
C**** final sweep to reinstate islands
c      do j=2,ny1-1
c        do i=2,nx1-1
c          sumk = nint(uvm(i,j)+uvm(i+1,j)+uvm(i,j+1)+uvm(i+1,j+1))
c          if (sumk.eq.4.and.heffm(i,j).eq.0.) then
c            uvm(i,j)=0 ; uvm(i+1,j)=0 ; uvm(i,j+1)=0 ; uvm(i+1,j+1)=0
c          end if
c        enddo
c      enddo
c set lateral boundary conditions
      do j=1,ny1
        heffm(1,j)   = heffm(nx1-1,j)
        heffm(nx1,j) = heffm(2,j)
      enddo
      do j=1,ny1-1
        do i=1,nx1-1
          uvm(i,j) = nint(min(heffm(i,j), heffm(i+1,j), heffm(i,j+1),
     *         heffm(i+1,j+1)))
        end do
      end do
c set cyclic conditions on eastern and western boundary
      do j=1,ny1
        uvm(1,j) = uvm(nx1-1,j)
        uvm(nx1,j) = uvm(2,j)
      enddo

c**** initialize heff and area
      do j=1,NY1
        do i=2,NX1-1
          HEFF(I,J,3)=0.
          HEFF(I,J,2)=0.
          AREA(I,J,3)=1.
          AREA(I,J,2)=1.
        enddo
        AREA(1,j,3)=1.
        AREA(1,j,2)=1.
        AREA(NX1,j,3)=1.
        AREA(NX1,j,2)=1.
      enddo
C**** sin/cos ice-ocean turning angle
      SINWAT=SIN(OIPHI)
      COSWAT=COS(OIPHI)

      END IF

C**** Start main loop
C**** Replicate polar boxes
      RSI(2:IM,JM)=RSI(1,JM)
      MSI(2:IM,JM)=MSI(1,JM)
      UO(2:IM,JM,1)=UO(1,JM,1)
      VO(1:IM,JM,1)=0.
      DMUA(2:IM,JM,2) = DMUA(1,JM,2)
      DMVA(2:IM,JM,2) = DMVA(1,JM,2)

      do j=1,NY1
        do i=2,NX1-1
          HEFF(I,J,1)=RSI(I-1,J)*RATIC(J)*(ACE1I+MSI(I-1,J))*BYRHOI
          AREA(I,J,1)=RSI(I-1,J)
c if heffm(i,j).eq.1 .and. focean(i-1,j).eq.0 interpolate.....
c          if (heffm(i,j).eq.1 .and. focean(i-1,j).eq.0) then
c     rsifix=(rsi(i-2,j)*focean(i-2,j)+rsi(i,j)*focean(i,j)+rsi(i-1,j+1)
c    *  *focean(i-1,j+1)+rsi(i-1,j-1)*focean(i-1,j-1))/(focean(i-2,j)
c    *  +focean(i,j)+focean(i-1,j+1)+focean(i-1,j-1))
c     msifix=(msi(i-2,j)*focean(i-2,j)+msi(i,j)*focean(i,j)+msi(i-1,j+1)
c    *  *focean(i-1,j+1)+msi(i-1,j-1)*focean(i-1,j-1))/(focean(i-2,j)
c    *  +focean(i,j)+focean(i-1,j+1)+focean(i-1,j-1))
c      HEFF(I,J,1)=rsifix**RATIC(J)*(ACE1I+msifix)*BYRHOI
c      AREA(I,J,1)=rsifix
c          end if
          HEFF(I,J,1)=HEFF(I,J,1)*HEFFM(I,J)
        enddo
      enddo
C**** fill in overlap regions
      DO J=1,NY1
       HEFF(1,J,1)=HEFF(NX1-1,J,1)
       AREA(1,J,1)=AREA(NX1-1,J,1)
       HEFF(NX1,J,1)=HEFF(2,J,1)
       AREA(NX1,J,1)=AREA(2,J,1)
      END DO

C**** save current value of sea ice concentration for ADVSI
      RSISAVE(:,:)=RSI(:,:)
C**** Calculate pressure anomaly at ocean surface (and scale for areas)
C**** OPRESS is on ocean grid
      DO J=1,JM
        DO I=1,IMAXJ(J)
          OPRESS(I,J) = RATOC(J)*(100.*(P(I,J)+PTOP-1013.25d0)+
     *         RSI(I,J)*(SNOWI(I,J)+ACE1I+MSI(I,J))*GRAV)
        END DO
      END DO
      OPRESS(2:IM,1)  = OPRESS(1,1)
      OPRESS(2:IM,JM) = OPRESS(1,JM)

C**** calculate sea surface tilt (incl. atmospheric pressure term)
C**** on ocean velocity grid
C**** PGF is an accelaration
      PGFU(1:IM,JM)=0
      DO J=2,JM-1
        I=IM   
        DO IP1=1,IM
          IF(LMU(I,J).gt.0. .and. RSI(I,J)+RSI(IP1,J).gt.0.) THEN
            PGFU(I,J)=-((OPRESS(IP1,J)-OPRESS(I,J))*BYRHOI
     *           +OGEOZ(IP1,J)-OGEOZ(I,J))/DXPO(J)
          ELSE
            PGFU(I,J)=0. 
          END IF
          I=IP1 
        END DO
      END DO
      DO J=1,JM-1 
        DO I=1,IM 
          IF(LMV(I,J).gt.0. .and. RSI(I,J)+RSI(I,J+1).gt.0.) THEN
            PGFV(I,J)=-((OPRESS(I,J+1)-OPRESS(I,J))*BYRHOI
     *           +OGEOZ(I,J+1)-OGEOZ(I,J))/DYVO(J)
          ELSE
            PGFV(I,J)=0. 
          END IF
        END DO
      END DO

c**** interpolate air, current and ice velocity from C grid to B grid
C**** This should be more generally from ocean grid to ice grid
      do j=1,jm-1
        im1=im
        do i=1,im
          UIB  (i,j)=0.5*(USI (im1,j)  +USI (im1,j+1))
          GWATX(i,j)=0.5*(UO  (im1,j,1)+UO  (im1,j+1,1))
          PGFUB(i,j)=0.5*(PGFU(im1,j)  +PGFU(im1,j+1))
          VIB  (i,j)=0.5*(VSI (im1,j)  +VSI (i,j))
          GWATY(i,j)=0.5*(VO  (im1,j,1)+VO  (i,j,1))
          PGFVB(i,j)=0.5*(PGFV(im1,j)  +PGFV(i,j))
          im1=i
        enddo
      enddo
c**** set north pole
      do i=1,im
        UIB  (i,jm)=USI(1,jm)
        GWATX(i,jm)=UO(1,jm,1)
        PGFUB(i,jm)=PGFU(1,jm)
        VIB  (i,jm)=0.
        GWATY(i,jm)=0.
        PGFVB(i,jm)=0.
      enddo
      DO J=1,NY1
        UIB(nx1-1,J)=UIB(1,J)
        VIB(nx1-1,J)=VIB(1,J)
        GWATX(nx1-1,J)=GWATX(1,J)
        GWATY(nx1-1,J)=GWATY(1,J)
        PGFUB(nx1-1,J)=PGFUB(1,J)
        PGFVB(nx1-1,J)=PGFVB(1,J)
        UIB(nx1,J)=UIB(2,J)
        VIB(nx1,J)=VIB(2,J)
        GWATX(nx1,J)=GWATX(2,J)
        GWATY(nx1,J)=GWATY(2,J)
        PGFUB(nx1,J)=PGFUB(2,J)
        PGFVB(nx1,J)=PGFVB(2,J)
      ENDDO

C**** DMUA is defined over the whole box (not just over ptype)
C**** Convert to stress over ice fraction only, & convert to ocean area
      DO J=1,JM
        DO I=1,IM
          IF (RSI(I,J).gt.0) THEN
            RAT=RATIC(J)/RSI(I,J)
            DMUA(I,J,2) = DMUA(I,J,2)*RAT
            DMVA(I,J,2) = DMVA(I,J,2)*RAT
          ELSE
            DMUA(I,J,2) = 0.
            DMVA(I,J,2) = 0.
          END IF
        END DO
      END DO

c**** interpolate air stress from A grid in atmos, to B grid in ocean
      do j=1,jm-1
        im1=im
        do i=1,im
          GAIRX(i,j)=0.25*(dmua(i,j,2)+dmua(im1,j,2)+dmua(im1,j+1,2)
     #         +dmua(i,j+1,2))*bydts
          GAIRY(i,j)=0.25*(dmva(i,j,2)+dmva(im1,j,2)+dmva(im1,j+1,2)
     #         +dmva(i,j+1,2))*bydts
          im1=i
        enddo
      enddo
      GAIRX(1:im,jm)=dmua(1,jm,2)*bydts
      GAIRY(1:im,jm)=dmva(1,jm,2)*bydts

      do j=1,ny1
       GAIRX(nx1-1,j)=GAIRX(1,j)
       GAIRY(nx1-1,j)=GAIRY(1,j)
       GAIRX(nx1,j)=GAIRX(2,j)
       GAIRY(nx1,j)=GAIRY(2,j)
      enddo

C****
C**** Set up mass per unit area and coriolis term
C****
      DO J=1,NY1-1
      DO I=1,NX1-1
        AMASS(I,J)=RHOI*0.25*(HEFF(I,J,1)
     *       +HEFF(I+1,J,1)+HEFF(I,J+1,1)+HEFF(I+1,J+1,1))
        COR(I,J)=AMASS(I,J)*2.0*OMEGA*SINEN(I,J)
      END DO
      END DO

c**** read in sea ice velocity
      DO J=1,NY1
      DO I=1,NX1
       UICE(I,J,1)=UIB(I,J)
       VICE(I,J,1)=VIB(I,J)
       UICE(I,J,2)=0.
       VICE(I,J,2)=0.
       UICE(I,J,3)=0.
       VICE(I,J,3)=0.
      END DO
      END DO

C KKI LOOP IS FOR PSEUDO-TIMESTEPPING
      KKI=0.
 10   KKI=KKI+1

C FIRST DO PREDICTOR
      DO J=1,NY1
      DO I=1,NX1
       UICE(I,J,3)=UICE(I,J,1)
       VICE(I,J,3)=VICE(I,J,1)
       UICEC(I,J)=UICE(I,J,1)
       VICEC(I,J)=VICE(I,J,1)
      END DO
      END DO

      CALL FORM(UICE,VICE,ETA,ZETA,DRAGS,DRAGA,GAIRX,GAIRY,GWATX,GWATY
     *     ,PGFUB,PGFVB,FORCEX,FORCEY,HEFF,AMASS,AREA)
      CALL RELAX(UICE,VICE,ETA,ZETA,DRAGS,DRAGA,AMASS,FORCEX,FORCEY
     1,UICEC,VICEC)

      DO J=1,NY1
       UICE(1,J,1)=UICE(NX1-1,J,1)
       VICE(1,J,1)=VICE(NX1-1,J,1)
       UICE(NX1,J,1)=UICE(2,J,1)
       VICE(NX1,J,1)=VICE(2,J,1)
      END DO

C NOW DO REGULAR TIME STEP
C NOW DO MODIFIED EULER STEP
c
      DO J=1,NY1
      DO I=1,NX1
       UICE(I,J,1)=0.5*(UICE(I,J,1)+UICE(I,J,2))
       VICE(I,J,1)=0.5*(VICE(I,J,1)+VICE(I,J,2))
      END DO
      END DO

      CALL FORM(UICE,VICE,ETA,ZETA,DRAGS,DRAGA,GAIRX,GAIRY,GWATX,GWATY
     *     ,PGFUB,PGFVB,FORCEX,FORCEY,HEFF,AMASS,AREA)
c
C NOW SET U(1)=U(2) AND SAME FOR V
      DO J=1,NY1
      DO I=1,NX1
       UICE(I,J,3)=UICE(I,J,1)
       VICE(I,J,3)=VICE(I,J,1)
       UICEC(I,J)=UICE(I,J,1)
       VICEC(I,J)=VICE(I,J,1)
       UICE(I,J,1)=UICE(I,J,2)
       VICE(I,J,1)=VICE(I,J,2)
      END DO
      END DO

      CALL RELAX(UICE,VICE,ETA,ZETA,DRAGS,DRAGA,AMASS,FORCEX,FORCEY
     1,UICEC,VICEC)

      DO J=1,NY1
       UICE(1,J,1)=UICE(NX1-1,J,1)
       VICE(1,J,1)=VICE(NX1-1,J,1)
       UICE(NX1,J,1)=UICE(2,J,1)
       VICE(NX1,J,1)=VICE(2,J,1)
      END DO

      if (kki.gt.1) then ! test convergence
        rms=0.
        do i=1,nx1
          do j=1,ny1
            rms=rms+(USAVE(i,j)-UICE(i,j,1))**2+(VSAVE(i,J)-VICE(i,j,1))
     *           **2
          end do
        end do
      end if

      if (kki.eq.20) then
        write(6,*) "Too many iterations in DYNSI. kki:",kki,rms
      elseif (kki.eq.1 .or. rms.gt.0.01d0) then
        USAVE=UICE(:,:,1)
        VSAVE=VICE(:,:,1)
        goto 10 
      end if

C**** Calculate stress on ice velocity grid
      DO J=1,NY1
        hemi=-1.
        if (J.gt.NY1/2) hemi=1.
        DO I=1,NX1
          DMU(i,j)=DTS*dwatn(i,j)*(COSWAT*(UICE(i,j,1)-GWATX(i,j))-
     *         HEMI*SINWAT*(VICE(i,j,1)-GWATY(i,j)))
          DMV(i,j)=DTS*dwatn(i,j)*(HEMI*SINWAT*(UICE(i,j,1)-GWATX(i,j))
     *         +COSWAT*(VICE(i,j,1)-GWATY(i,j)))
        END DO
      END DO

C**** interpolate ice velocity and stress from B grid to C grid in ocean
      do j=2,jm
        i=im
        do ip1=1,im
          usi(i,j)=0.5*(uice(i+1,j-1,1)+uice(i+1,j,1))
          IF (abs(USI(I,J)).lt.1d-10) USI(I,J)=0
          DMUI(I,J) = 0.5*(dmu(i+1,j-1)+dmu(i+1,j))
C**** Rescale DMUI to be net momentum into ocean
          DMUI(I,J) = 0.5*DMUI(I,J)*(RSI(I,J)+RSI(IP1,J))
          i=ip1
        enddo
      enddo
      do j=1,jm-1
        do i=1,im
          vsi(i,j)=0.5*(vice(i,j,1)+vice(i+1,j,1))
          IF (abs(VSI(I,J)).lt.1d-10) VSI(I,J)=0
          DMVI(I,J) = 0.5*(dmv(i,j)+dmv(i+1,j))
C**** Rescale DMVI to be net momentum into ocean
          IF (J.lt.JM-1) DMVI(I,J) = 0.5*DMVI(I,J)*(RSI(I,J)*DXYNO(J)
     *         +RSI(I,J+1)*DXYSO(J+1))/DXYVO(J)
          IF (J.eq.JM-1) DMVI(I,JM-1) = 0.5*DMVI(I,JM-1)*(RSI(I,JM-1)
     *         *DXYNO(JM-1)+RSI(1,JM)*DXYSO(JM))/DXYVO(JM-1)
        enddo
      enddo
      usi(1:im,1)=0.
      vsi(1:im,jm)=0.
      dmui(1:im,1)=0.
      dmvi(1:im,jm)=0.

C**** Calculate ustar*2*rho for ice-ocean fluxes on atmosphere grid
C**** UI2rho = | tau |
      do j=1,jm
        do i=1,imaxj(j)
          UI2rho(i,j)=0
          if (FOCEAN(I,J)*RSI(i,j).gt.0) THEN
C**** calculate 4 point average of B grid values of stresses (CHECK!)
            duA = 0.5*(DXYNO(J)*(dmu(i+1,j)+dmu(i,j))+DXYSO(j)*(dmu(i+1
     *           ,j-1)+dmu(i,j-1)))*BYDXYPO(J)   
            dvA = 0.5*(DXYNO(J)*(dmv(i+1,j)+dmv(i,j))+DXYSO(j)*(dmv(i+1
     *           ,j-1)+dmv(i,j-1)))*BYDXYPO(J)   
            UI2rho(i,j)= sqrt (duA**2 + dvA**2) * bydts
          end if
        end do
      end do

C**** set north pole in C grid
      USINP=0.
      DMUINP=0.
      do i=1,im
        USINP = USINP + USI(i,jm)
        DMUINP = DMUINP + DMUI(i,jm)
      enddo
      USINP=USINP/IM
      DMUINP=DMUINP/IM
      USI(1:im,jm)=USINP
      DMUI(1:im,jm)=DMUINP
      VSI(1:im,jm)=0.
      DMVI(1:im,jm)=0.

C**** calculate mass fluxes for the ice advection
      DO J=1,JM-1
        I=IM
        DO IP1=1,IM
          USIDT(I,J)=0.
          IF (LMU(I,J).gt.0. .and. RSI(I,J)+RSI(IP1,J).gt.1d-4) THEN
            USIDT(I,J)=USI(I,J)*DTS
            OIJ(I,J,IJ_USI) =OIJ(I,J,IJ_USI) +(RSI(I,J)+RSI(IP1,J))
     *           *USI(i,j)
            OIJ(I,J,IJ_DMUI)=OIJ(I,J,IJ_DMUI)+DMUI(i,j)
          END IF
          VSIDT(I,J)=0.
          IF (LMV(I,J).gt.0. .and. RSI(I,J)+RSI(I,J+1).gt.1d-4) THEN
            VSIDT(I,J)=VSI(I,J)*DTS
            OIJ(I,J,IJ_VSI) =OIJ(I,J,IJ_VSI) +(RSI(I,J)+RSI(I,J+1))
     *           *VSI(i,j)
            OIJ(I,J,IJ_DMVI)=OIJ(I,J,IJ_DMVI)+DMVI(i,j)
          END IF
          OIJ(I,J,IJ_PICE)=OIJ(I,J,IJ_PICE)+ RSI(I,J)*press(i+1,j)
c         OIJ(I,J,IJ_RSI) =OIJ(I,J,IJ_RSI) + RSI(I,J)
          I=IP1
        END DO
      END DO
      VSIDT(1:IM,JM)=0.
      USIDT(1:IM,JM)=USI(1,JM)*DTS
      OIJ(1,JM,IJ_USI) =OIJ(1,JM,IJ_USI) +RSI(1,JM)*USI(1,JM)
      OIJ(1,JM,IJ_DMUI)=OIJ(1,JM,IJ_DMUI)+DMUI(1,JM)
c     OIJ(1,JM,IJ_RSI)=OIJ(1,JM,IJ_RSI)+ RSI(1,JM)
C****
      END SUBROUTINE DYNSI

C*************************************************************
C PLEASE KEEP THIS NOTE OF MODEL-DEVELOPMENT HISTORY
C Matrix solve uses Thomas algorithm, 10/1991, Jinlun Zhang
C Spherical coordinate system, 10/27/93, Jinlun Zhang
C Latest finite differencing scheme for treatment of NP,
C 9/9/1996,Jinlun Zhang
C Alternating direction implicit (ADI) method is used, 10/1998,
C Jinlun Zhang
C For details about ADI dynamics model, see Zhang and Rothrock,
C "Modeling
C Arctic Sea ice with an efficient plastic solution",
C submitted to JGR, 1999
C Adapted for GISS coupled model Jiping Liu/Gavin Schmidt 2000
C*************************************************************

      SUBROUTINE FORM(UICE,VICE,ETA,ZETA,DRAGS,DRAGA,GAIRX,GAIRY,
     *     GWATX,GWATY,PGFUB,PGFVB,FORCEX,FORCEY,HEFF,AMASS,AREA)
!@sum  FORM calculates ice dynamics input parameters for relaxation
!@auth Jiping Liu/Gavin Schmidt (based on code from J. Zhang)
!@ver  1.0
      USE ICEDYN
      IMPLICIT NONE
      REAL*8, DIMENSION(NX1,NY1,3), INTENT(INOUT) :: UICE,VICE,HEFF,AREA
      REAL*8, DIMENSION(NX1,NY1), INTENT(INOUT) :: ETA,ZETA,DRAGS,DRAGA
     *     ,GAIRX,GAIRY,GWATX,GWATY,FORCEX,FORCEY,AMASS,PGFUB,PGFVB
      REAL*8, DIMENSION(NX1,NY1) :: ZMAX,ZMIN
      INTEGER I,J
      REAL*8 AAA
C****
C**** Set up non linear water drag
C****
      DO J=1,NY1-1
      DO I=1,NX1-1
        DWATN(I,J)=5.5*SQRT((UICE(I,J,1)-GWATX(I,J))**2
     1       +(VICE(I,J,1)-GWATY(I,J))**2)
      END DO
      END DO
C NOW SET UP SYMMETRIC DRAG
      DO J=1,NY1-1
      DO I=1,NX1-1
        DRAGS(I,J)=DWATN(I,J)*COSWAT
      END DO
      END DO
C NOW SET UP ANTI SYMMETRIC DRAG PLUS CORIOLIS
      DO J=1,NY1
      DO I=1,NX1
        IF(J.GT.NY1/2) THEN
          DRAGA(I,J)=DWATN(I,J)*SINWAT+COR(I,J)
        ELSE
          DRAGA(I,J)=DWATN(I,J)*(-SINWAT)+COR(I,J)
        END IF
      END DO
      END DO
C NOW SET UP FORCING FIELD
      DO J=1,NY1
      DO I=1,NX1

C FIRST DO WIND
       FORCEX(I,J)=GAIRX(i,j)
       FORCEY(I,J)=GAIRY(i,j)

C NOW ADD IN CURRENT FORCE
       IF(J.GT.NY1/2) THEN
         FORCEX(I,J)=FORCEX(I,J)+DWATN(I,J)*(COSWAT*GWATX(I,J)
     1        -SINWAT*GWATY(I,J))
         FORCEY(I,J)=FORCEY(I,J)+DWATN(I,J)*(SINWAT*GWATX(I,J)
     1        +COSWAT*GWATY(I,J))
       ELSE
         FORCEX(I,J)=FORCEX(I,J)+DWATN(I,J)*(COSWAT*GWATX(I,J)
     1        +SINWAT*GWATY(I,J))
         FORCEY(I,J)=FORCEY(I,J)+DWATN(I,J)*(-SINWAT*GWATX(I,J)
     1        +COSWAT*GWATY(I,J))
       END IF

C NOW ADD IN TILT 
C**** This assumes explicit knowledge of sea surface tilt
       FORCEX(I,J)=FORCEX(I,J)+AMASS(I,J)*PGFUB(I,J)
       FORCEY(I,J)=FORCEY(I,J)+AMASS(I,J)*PGFVB(I,J)
C**** Otherwise estimate tilt using geostrophy
c       FORCEX(I,J)=FORCEX(I,J)-COR(I,J)*GWATY(I,J)
c       FORCEY(I,J)=FORCEY(I,J)+COR(I,J)*GWATX(I,J)

      END DO
      END DO

C NOW SET UP ICE PRESSURE AND VISCOSITIES
      DO J=1,NY1
      DO I=1,NX1
        PRESS(I,J)=PSTAR*HEFF(I,J,1)*EXP(-20.0*(1.0-AREA(I,J,1)))
        ZMAX(I,J)=(5d12/2d4)*PRESS(I,J)
c       ZMIN(I,J)=0.0D+00
        ZMIN(I,J)=4d8
      END DO
      END DO

      CALL PLAST(UICE,VICE,PRESS,ETA,ZETA,ZMAX,ZMIN)

      AAA=0.0
      DO I=2,NX1-1
        AAA=AAA+PRESS(I,NY1-1)
      END DO
      AAA=AAA/FLOAT(NX1-2)
      DO I=1,NX1
        PRESS(I,NY1)=AAA
      END DO

 8481 CONTINUE

      DO J=1,NY1
        PRESS(1,J)=PRESS(NX1-1,J)
        PRESS(NX1,J)=PRESS(2,J)
      END DO

C NOW SET VISCOSITIES AND PRESSURE EQUAL TO ZERO AT OUTFLOW PTS

      DO J=1,NY1
      DO I=1,NX1
        PRESS(I,J)=PRESS(I,J)*HEFFM(I,J)
        ETA(I,J)=ETA(I,J)*HEFFM(I,J)
        ZETA(I,J)=ZETA(I,J)*HEFFM(I,J)
      END DO
      END DO

C NOW CALCULATE PRESSURE FORCE AND ADD TO EXTERNAL FORCE
      DO J=1,NY1-1
      DO I=1,NX1-1
      FORCEX(I,J)=FORCEX(I,J)-(0.25/(DXU(I)*CSU(J)))
     1  *(PRESS(I+1,J)+PRESS(I+1,J+1)-PRESS(I,J)-PRESS(I,J+1))
      FORCEY(I,J)=FORCEY(I,J)-0.25/DYU(J)
     1  *(PRESS(I,J+1)+PRESS(I+1,J+1)-PRESS(I,J)-PRESS(I+1,J))
C NOW PUT IN MINIMAL MASS FOR TIME STEPPING CALCULATIONS
      END DO
      END DO

      DO J=1,NY1
        FORCEX(1,J)=FORCEX(NX1-1,J)
        FORCEY(1,J)=FORCEY(NX1-1,J)
        FORCEX(NX1,J)=FORCEX(2,J)
        FORCEY(NX1,J)=FORCEY(2,J)
      END DO

      RETURN
      END SUBROUTINE FORM

      SUBROUTINE PLAST(UICE,VICE,PRESS,ETA,ZETA,ZMAX,ZMIN)
!@sum  PLAST Calculates strain rates and viscosity for dynamic ice
!@auth Jiping Liu/Gavin Schmidt (based on code from J. Zhang)
!@ver  1.0
      USE CONSTANT, only : radius
      USE ICEDYN, only : eccen,dxt,cst,dyt,tngt,nx1,ny1
      IMPLICIT NONE
      REAL*8, DIMENSION(NX1,NY1,3), INTENT(IN) :: UICE,VICE
      REAL*8, DIMENSION(NX1,NY1), INTENT(IN) :: PRESS,ZMAX,ZMIN
      REAL*8, DIMENSION(NX1,NY1), INTENT(OUT) :: ETA,ZETA
      REAL*8, DIMENSION(NX1,NY1) :: E11,E22,E12
c      REAL*8 :: SS11
      REAL*8, PARAMETER :: ECM2 = 1.0/(ECCEN**2),GMIN=1d-20
      REAL*8 DELT,DELT1,AAA
      INTEGER I,J

C EVALUATE STRAIN RATES
      DO J=2,NY1-1
        DO I=2,NX1-1
          E11(I,J)=0.5/(DXT(I)*CST(J))*(UICE(I,J,1)+UICE(I,J-1,1)
     *         -UICE(I-1,J,1)-UICE(I-1,J-1,1))-0.25*(VICE(I,J,1)+
     *         VICE(I-1,J,1)+VICE(I-1,J-1,1)+VICE(I,J-1,1))*TNGT(J)
     *         /RADIUS
          E22(I,J)=0.5/DYT(J)*(VICE(I,J,1)+VICE(I-1,J,1)
     *         -VICE(I,J-1,1)-VICE(I-1,J-1,1))
          E12(I,J)=0.5*(0.5/DYT(J)*(UICE(I,J,1)+UICE(I-1,J,1)-
     *         UICE(I,J-1,1)-UICE(I-1,J-1,1))+0.5/(DXT(I)*CST(J))*
     *         (VICE(I,J,1)+VICE(I,J-1,1)-VICE(I-1,J,1)-VICE(I-1,J-1,1))
     *         +0.25*(UICE(I,J,1)+UICE(I-1,J,1)+UICE(I-1,J-1,1)+UICE(I,J
     *         -1,1))*TNGT(J)/RADIUS)
C NOW EVALUATE VISCOSITIES
          DELT=(E11(I,J)**2+E22(I,J)**2)*(1.0+ECM2)+4.0*ECM2*E12(I,J)**2
     *         +2.0*E11(I,J)*E22(I,J)*(1.0-ECM2)
          DELT1=SQRT(DELT)
          DELT1=MAX(GMIN,DELT1)
          ZETA(I,J)=0.5*PRESS(I,J)/DELT1
        END DO
      END DO

C NOW PUT MIN AND MAX VISCOSITIES IN
      DO J=1,NY1
        DO I=1,NX1
          ZETA(I,J)=MIN(ZMAX(I,J),ZETA(I,J))
          ZETA(I,J)=MAX(ZMIN(I,J),ZETA(I,J))
        END DO
      END DO

      AAA=0.0
      DO I=2,NX1-1
        AAA=AAA+ZETA(I,NY1-1)
      END DO
      AAA=AAA/FLOAT(NX1-2)
      DO I=1,NX1
        ZETA(I,NY1)=AAA
      END DO

      DO J=1,NY1
        ZETA(1,J)=ZETA(NX1-1,J)
        ZETA(NX1,J)=ZETA(2,J)
      END DO

      DO J=1,NY1
        DO I=1,NX1
          ETA(I,J)=ECM2*ZETA(I,J)
c         E11(I,J)=E11(I,J)*HEFFM(I,J)
c         E22(I,J)=E22(I,J)*HEFFM(I,J)
c         E12(I,J)=E12(I,J)*HEFFM(I,J)
c         SS11=(ZETA(I,J)-ETA(I,J))*(E11(I,J)+E22(I,J))-PRESS(I,J)*0.5
        END DO
      END DO
      RETURN
      END SUBROUTINE PLAST

      SUBROUTINE RELAX(UICE,VICE,ETA,ZETA,DRAGS,DRAGA
     *     ,AMASS,FORCEX,FORCEY,UICEC,VICEC)
!@sum  RELAX calculates ice dynamics relaxation method
!@auth Jiping Liu/Gavin Schmidt (based on code from J. Zhang)
!@ver  1.0
      USE CONSTANT, only : radius
      USE ICEDYN
      IMPLICIT NONE
      REAL*8, DIMENSION(NX1,NY1,3), INTENT(INOUT) :: UICE,VICE
      REAL*8, DIMENSION(NX1,NY1), INTENT(INOUT) :: ETA,ZETA,DRAGS,DRAGA
     *     ,FORCEX,FORCEY,AMASS,UICEC,VICEC
      REAL*8, DIMENSION(NX1,NY1) :: AU,BU,CU,AV,BV,CV,FXY,FXY1
      REAL*8, DIMENSION(NX1) :: CUU,URT
      REAL*8, DIMENSION(NY1) :: CVV,VRT
      INTEGER, PARAMETER :: NYPOLE=NY1-1,NPOL=1,NXLCYC=NX1-1,LCYC=1
      REAL*8, PARAMETER :: BYRAD2 = 1./(RADIUS*RADIUS)
      INTEGER I,J,J1,J2,IMD,JMD
      REAL*8 DELXY,DELXR,DELX2,DELY2,DELYR,ETAMEAN,ZETAMEAN,AA1,AA2
     *     ,AA3,AA4,AA5,AA6,AA9

      DO J=1,NY1
        DO I=1,NX1
          FORCEX(I,J)=FORCEX(I,J)*UVM(I,J)
          FORCEY(I,J)=FORCEY(I,J)*UVM(I,J)
        END DO
      END DO
C MUST UPDATE HEFF BEFORE CALLING RELAC
C FIRST SET U(2)=U(1)
      DO J=1,NY1
        DO I=1,NX1
C NOW MAKE SURE BDRY PTS ARE EQUAL TO ZERO
          UICE(I,J,2)=UICE(I,J,1)
          VICE(I,J,2)=VICE(I,J,1)
          UICE(I,J,1)=UICE(I,J,3)*UVM(I,J)
          VICE(I,J,1)=VICE(I,J,3)*UVM(I,J)
        END DO
      END DO

      DO I=1,NX1/2
      UICE(I,NY1,1)=-UICEC(I+(NX1-2)/2,NY1-1)
      VICE(I,NY1,1)=-VICEC(I+(NX1-2)/2,NY1-1)
      UICE(I,NY1,3)=-UICEC(I+(NX1-2)/2,NY1-1)
      VICE(I,NY1,3)=-VICEC(I+(NX1-2)/2,NY1-1)

      UICEC(I,NY1)=-UICEC(I+(NX1-2)/2,NY1-1)
      VICEC(I,NY1)=-VICEC(I+(NX1-2)/2,NY1-1)
      END DO

      DO I=NX1/2+1,NX1-1
      UICE(I,NY1,1)=-UICEC(I-(NX1-2)/2,NY1-1)
      VICE(I,NY1,1)=-VICEC(I-(NX1-2)/2,NY1-1)
      UICE(I,NY1,3)=-UICEC(I-(NX1-2)/2,NY1-1)
      VICE(I,NY1,3)=-VICEC(I-(NX1-2)/2,NY1-1)

      UICEC(I,NY1)=-UICEC(I-(NX1-2)/2,NY1-1)
      VICEC(I,NY1)=-VICEC(I-(NX1-2)/2,NY1-1)
      END DO

      DO J=1,NY1
      UICE(1,J,1)=UICEC(NX1-1,J)
      VICE(1,J,1)=VICEC(NX1-1,J)
      UICE(NX1,J,1)=UICEC(2,J)
      VICE(NX1,J,1)=VICEC(2,J)
      UICE(1,J,3)=UICEC(NX1-1,J)
      VICE(1,J,3)=VICEC(NX1-1,J)
      UICE(NX1,J,3)=UICEC(2,J)
      VICE(NX1,J,3)=VICEC(2,J)

      UICEC(1,J)=UICEC(NX1-1,J)
      VICEC(1,J)=VICEC(NX1-1,J)
      UICEC(NX1,J)=UICEC(2,J)
      VICEC(NX1,J)=VICEC(2,J)
      END DO

C FIRST DO UICE
C THE FIRST HALF

      DO J=2,NYPOLE
      DO I=2,NXLCYC
      DELXY=BYDXDY(I,J)    ! 0.5/(DXU(I)*DYU(J))
      DELXR=BYDXR(I)       ! 0.5/(DXU(I)*RADIUS)
      DELX2=BYDX2(I)       ! 0.5/(DXU(I)*DXU(I))
      ETAMEAN=0.25*(ETA(I,J+1)+ETA(I+1,J+1)+ETA(I,J)+ETA(I+1,J))
      ZETAMEAN=0.25*(ZETA(I,J+1)+ZETA(I+1,J+1)+ZETA(I,J)+ZETA(I+1,J))

      FXY(I,J)=DRAGA(I,J)*VICEC(I,J)+FORCEX(I,J)
     3+0.5*(ZETA(I+1,J+1)*(VICEC(I+1,J+1)+VICEC(I,J+1)
     3-VICEC(I+1,J)-VICEC(I,J))+ZETA(I+1,J)*(VICEC(I+1,J)
     3+VICEC(I,J)-VICEC(I+1,J-1)-VICEC(I,J-1))+ZETA(I,J+1)
     3*(VICEC(I,J)+VICEC(I-1,J)-VICEC(I,J+1)-VICEC(I-1,J+1))
     3+ZETA(I,J)*(VICEC(I,J-1)+VICEC(I-1,J-1)-VICEC(I,J)
     3-VICEC(I-1,J)))*DELXY/CSU(J)
     3
     4-0.5*(ETA(I+1,J+1)*(VICEC(I+1,J+1)+VICEC(I,J+1)
     4-VICEC(I+1,J)-VICEC(I,J))+ETA(I+1,J)*(VICEC(I+1,J)
     4+VICEC(I,J)-VICEC(I+1,J-1)-VICEC(I,J-1))+ETA(I,J+1)
     4*(VICEC(I,J)+VICEC(I-1,J)-VICEC(I,J+1)-VICEC(I-1,J+1))
     4+ETA(I,J)*(VICEC(I,J-1)+VICEC(I-1,J-1)-VICEC(I,J)
     4-VICEC(I-1,J)))*DELXY/CSU(J)
     4
     5+0.5*(VICEC(I+1,J)-VICEC(I-1,J))*(ETA(I,J+1)+ETA(I+1,J+1)
     5-ETA(I,J)-ETA(I+1,J))*DELXY/CSU(J)+0.5*ETAMEAN*((VICEC(I+1,J+1)
     5-VICEC(I-1,J+1))/CSU(J+1)-(VICEC(I+1,J-1)-VICEC(I-1,J-1))
     5/CSU(J-1))*DELXY
     5
     6-((ZETA(I+1,J+1)+ZETA(I+1,J)-ZETA(I,J)-ZETA(I,J+1))
     6+(ETA(I+1,J+1)+ETA(I+1,J)-ETA(I,J)-ETA(I,J+1)))
     6*TNG(J)*VICEC(I,J)*DELXR/CSU(J)
     6-(ETAMEAN+ZETAMEAN)*TNG(J)*(VICEC(I+1,J)-VICEC(I-1,J))
     6*DELXR/CSU(J)
     6
     7-ETAMEAN*2.0*TNG(J)*(VICEC(I+1,J)-VICEC(I-1,J))*DELXR/CSU(J)

      END DO
      END DO

      DO J=2,NYPOLE
      DO I=2,NXLCYC
      DELX2=BYDX2(I)       ! 0.5/(DXU(I)*DXU(I))
      DELY2=BYDY2(J)       ! 0.5/(DYU(J)*DYU(J))
      DELYR=BYDYR(J)       ! 0.5/(DYU(J)*RADIUS)
      ETAMEAN=0.25*(ETA(I,J+1)+ETA(I+1,J+1)+ETA(I,J)+ETA(I+1,J))
      AA1=((ETA(I+1,J)+ZETA(I+1,J))/CSU(J)+(ETA(I+1,J+1)+ZETA(I+1,J+1))
     &     /CSU(J))/CSU(J)
      AA2=((ETA(I,J)+ZETA(I,J))/CSU(J)+(ETA(I,J+1)+ZETA(I,J+1))
     &     /CSU(J))/CSU(J)
      AA3=ETA(I,J+1)+ETA(I+1,J+1)
      AA4=ETA(I,J)+ETA(I+1,J)
      AA5=-(ETA(I,J+1)+ETA(I+1,J+1)-ETA(I,J)-ETA(I+1,J))*TNG(J)
      AA6=2.0*ETAMEAN*TNG(J)*TNG(J)
      AU(I,J)=-AA2*DELX2*UVM(I,J)
      BU(I,J)=((AA1+AA2)*DELX2+AA6*BYRAD2
     &+AMASS(I,J)*BYDTS*2.0+DRAGS(I,J))*UVM(I,J)+(1.0-UVM(I,J))
      CU(I,J)=-AA1*DELX2*UVM(I,J)
      END DO
      END DO
      DO J=2,NYPOLE
      AU(2,J)=0.0
      CU(NXLCYC,J)=0.0
      CU(2,J)=CU(2,J)/BU(2,J)
      END DO

      DO 1200 J=2,NYPOLE
      DO I=2,NXLCYC

      DELX2=BYDX2(I)       ! 0.5/(DXU(I)*DXU(I))
      DELY2=BYDY2(J)       ! 0.5/(DYU(J)*DYU(J))
      DELYR=BYDYR(J)       ! 0.5/(DYU(J)*RADIUS)
      ETAMEAN=0.25*(ETA(I,J+1)+ETA(I+1,J+1)+ETA(I,J)+ETA(I+1,J))

      AA1=((ETA(I+1,J)+ZETA(I+1,J))/CSU(J)+(ETA(I+1,J+1)+ZETA(I+1,J+1))
     &     /CSU(J))/CSU(J)
      AA2=((ETA(I,J)+ZETA(I,J))/CSU(J)+(ETA(I,J+1)+ZETA(I,J+1))
     &     /CSU(J))/CSU(J)
      AA3=ETA(I,J+1)+ETA(I+1,J+1)
      AA4=ETA(I,J)+ETA(I+1,J)
      AA5=-(ETA(I,J+1)+ETA(I+1,J+1)-ETA(I,J)-ETA(I+1,J))*TNG(J)
      AA6=2.0*ETAMEAN*TNG(J)*TNG(J)

      IF(I.EQ.2) THEN
      AA9=AA2*DELX2*UICEC(I-1,J)*UVM(I,J)*FLOAT(LCYC-0)
      ELSE IF(I.EQ.NXLCYC) THEN
      AA9=AA1*DELX2*UICEC(I+1,J)*UVM(I,J)*FLOAT(LCYC-0)
      ELSE
      AA9=0.0
      END IF

      URT(I)=AA9+FXY(I,J)-AA5*DELYR*UICE(I,J,2)
     1-(AA3+AA4)*DELY2*UICE(I,J,2)
     1+(ETA(I,J+1)+ETA(I+1,J+1))*UICE(I,J+1,2)*DELY2
     2+(ETA(I,J)+ETA(I+1,J))*UICE(I,J-1,2)*DELY2
     3+ETAMEAN*DELYR*(UICE(I,J+1,2)*TNG(J+1)-UICE(I,J-1,2)*TNG(J-1))
     4-ETAMEAN*DELYR*2.0*TNG(J)*(UICE(I,J+1,2)-UICE(I,J-1,2))
      URT(I)=(URT(I)+AMASS(I,J)*BYDTS*UICE(I,J,2)*2.0)*UVM(I,J)
      END DO

      DO I=2,NXLCYC
      CUU(I)=CU(I,J)
      END DO
      URT(2)=URT(2)/BU(2,J)
      DO I=3,NXLCYC
      IMD=I-1
      CUU(I)=CUU(I)/(BU(I,J)-AU(I,J)*CUU(IMD))
      URT(I)=(URT(I)-AU(I,J)*URT(IMD))/(BU(I,J)-AU(I,J)*CUU(IMD))
      END DO
      DO I=1,NXLCYC-2
      J1=NXLCYC-I
      J2=J1+1
      URT(J1)=URT(J1)-CUU(J1)*URT(J2)
      END DO
      DO I=2,NXLCYC
      UICE(I,J,1)=URT(I)
      END DO
 1200 CONTINUE

      DO I=2,NXLCYC
      DO J=2,NYPOLE
      UICE(I,J,3)=UICE(I,J,1)
      END DO
      END DO

C NOW THE SECOND HALF
      DO I=2,NXLCYC
      DO J=2,NYPOLE
      DELX2=BYDX2(I)       ! 0.5/(DXU(I)*DXU(I))
      DELY2=BYDY2(J)       ! 0.5/(DYU(J)*DYU(J))
      DELYR=BYDYR(J)       ! 0.5/(DYU(J)*RADIUS)
      ETAMEAN=0.25*(ETA(I,J+1)+ETA(I+1,J+1)+ETA(I,J)+ETA(I+1,J))
      ZETAMEAN=0.25*(ZETA(I,J+1)+ZETA(I+1,J+1)+ZETA(I,J)+ZETA(I+1,J))

      AA1=ETA(I,J+1)+ETA(I+1,J+1)
      AA2=ETA(I,J)+ETA(I+1,J)
      AA5=-(ETA(I,J+1)+ETA(I+1,J+1)-ETA(I,J)-ETA(I+1,J))*TNG(J)
      AA6=2.0*ETAMEAN*TNG(J)*TNG(J)

      AV(I,J)=(-AA2*DELY2+ETAMEAN*DELYR*(TNG(J-1)-2.0*TNG(J)))*UVM(I,J)
      BV(I,J)=((AA1+AA2)*DELY2+AA5*DELYR+AA6*BYRAD2
     &+AMASS(I,J)*BYDTS*2.0+DRAGS(I,J))*UVM(I,J)+(1.0-UVM(I,J))
      CV(I,J)=(-AA1*DELY2-ETAMEAN*DELYR*(TNG(J+1)-2.0*TNG(J)))*UVM(I,J)
      END DO
      END DO
      DO I=2,NXLCYC
      AV(I,2)=0.0
      CV(I,NYPOLE)=0.0
      CV(I,2)=CV(I,2)/BV(I,2)
      END DO

      DO I=2,NXLCYC
      DO J=2,NYPOLE
      DELX2=BYDX2(I)       ! 0.5/(DXU(I)*DXU(I))
      DELY2=BYDY2(J)       ! 0.5/(DYU(J)*DYU(J))
      DELYR=BYDYR(J)       ! 0.5/(DYU(J)*RADIUS)
      ETAMEAN=0.25*(ETA(I,J+1)+ETA(I+1,J+1)+ETA(I,J)+ETA(I+1,J))

      AA1=((ETA(I+1,J)+ZETA(I+1,J))/CSU(J)+(ETA(I+1,J+1)+ZETA(I+1,J+1))
     &     /CSU(J))/CSU(J)
      AA2=((ETA(I,J)+ZETA(I,J))/CSU(J)+(ETA(I,J+1)+ZETA(I,J+1))
     &     /CSU(J))/CSU(J)

      IF(J.EQ.NYPOLE) THEN
      AA9=( (ETA(I,J+1)+ETA(I+1,J+1))*DELY2*UICEC(I,J+1)
     &     +ETAMEAN*DELYR*(TNG(J+1)-2.0*TNG(J))*UICEC(I,J+1) )*UVM(I,J)
     &*FLOAT(NPOL-0)
      ELSE
      AA9=0.0
      END IF

      FXY1(I,J)=AA9+AMASS(I,J)*BYDTS*UICE(I,J,1)*2.0
     5-(AA1+AA2)*DELX2*UICE(I,J,1)
     6+((ETA(I+1,J)+ZETA(I+1,J)+ETA(I+1,J+1)+ZETA(I+1,J+1))
     6*UICE(I+1,J,1)
     6+(ETA(I,J)+ZETA(I,J)+ETA(I,J+1)+ZETA(I,J+1))*UICE(I-1,J,1))
     6*DELX2/CSU(J)/CSU(J)

      END DO
      END DO

      DO 1300 I=2,NXLCYC
      DO J=2,NYPOLE
      VRT(J)=FXY(I,J)+FXY1(I,J)
      VRT(J)=VRT(J)*UVM(I,J)
      END DO

      DO J=2,NYPOLE
      CVV(J)=CV(I,J)
      END DO
      VRT(2)=VRT(2)/BV(I,2)
      DO J=3,NYPOLE
      JMD=J-1
      CVV(J)=CVV(J)/(BV(I,J)-AV(I,J)*CVV(JMD))
      VRT(J)=(VRT(J)-AV(I,J)*VRT(JMD))/(BV(I,J)-AV(I,J)*CVV(JMD))
      END DO
      DO J=1,NYPOLE-2
      J1=NYPOLE-J
      J2=J1+1
      VRT(J1)=VRT(J1)-CVV(J1)*VRT(J2)
      END DO
      DO J=2,NYPOLE
      UICE(I,J,1)=VRT(J)
      END DO
 1300 CONTINUE

C NOW DO VICE
C THE FIRST HALF

      DO I=2,NXLCYC
      DO J=2,NYPOLE
      DELXY=BYDXDY(I,J)    ! 0.5/(DXU(I)*DYU(J))
      DELXR=BYDXR(I)       ! 0.5/(DXU(I)*RADIUS)
      DELY2=BYDY2(J)       ! 0.5/(DYU(J)*DYU(J))
      DELYR=BYDYR(J)       ! 0.5/(DYU(J)*RADIUS)
      ETAMEAN=0.25*(ETA(I,J+1)+ETA(I+1,J+1)+ETA(I,J)+ETA(I+1,J))
      ZETAMEAN=0.25*(ZETA(I,J+1)+ZETA(I+1,J+1)+ZETA(I,J)+ZETA(I+1,J))

      FXY(I,J)=-DRAGA(I,J)*UICEC(I,J)+FORCEY(I,J)
     3+(0.5*(UICEC(I+1,J)-UICEC(I-1,J))*(ZETA(I,J+1)+ZETA(I+1,J+1)
     3-ZETA(I,J)-ZETA(I+1,J))*DELXY/CSU(J)+0.5*ZETAMEAN*
     3((UICEC(I+1,J+1)
     3-UICEC(I-1,J+1))/CSU(J+1)-(UICEC(I+1,J-1)-UICEC(I-1,J-1))
     3/CSU(J-1))*DELXY)
     3
     4-(0.5*(UICEC(I+1,J)-UICEC(I-1,J))*(ETA(I,J+1)+ETA(I+1,J+1)
     4-ETA(I,J)-ETA(I+1,J))*DELXY/CSU(J)+0.5*ETAMEAN*((UICEC(I+1,J+1)
     4-UICEC(I-1,J+1))/CSU(J+1)-(UICEC(I+1,J-1)-UICEC(I-1,J-1))
     4/CSU(J-1))*DELXY)
     4
     5+0.5*(ETA(I+1,J+1)*(UICEC(I+1,J+1)+UICEC(I,J+1)
     5-UICEC(I+1,J)-UICEC(I,J))+ETA(I+1,J)*(UICEC(I+1,J)
     5+UICEC(I,J)-UICEC(I+1,J-1)-UICEC(I,J-1))+ETA(I,J+1)
     5*(UICEC(I,J)+UICEC(I-1,J)-UICEC(I,J+1)-UICEC(I-1,J+1))
     5+ETA(I,J)*(UICEC(I,J-1)+UICEC(I-1,J-1)-UICEC(I,J)
     5-UICEC(I-1,J)))*DELXY/CSU(J)
     5
     6+(ETA(I+1,J+1)+ETA(I+1,J)-ETA(I,J)-ETA(I,J+1))
     6*TNG(J)*UICEC(I,J)*DELXR/CSU(J)
     6+ETAMEAN*TNG(J)*(UICEC(I+1,J)-UICEC(I-1,J))*DELXR/CSU(J)
     6
     7+ETAMEAN*2.0*TNG(J)*(UICEC(I+1,J)-UICEC(I-1,J))*DELXR/CSU(J)

      END DO
      END DO

      DO I=2,NXLCYC
      DO J=2,NYPOLE
      DELX2=BYDX2(I)       ! 0.5/(DXU(I)*DXU(I))
      DELY2=BYDY2(J)       ! 0.5/(DYU(J)*DYU(J))
      DELYR=BYDYR(J)       ! 0.5/(DYU(J)*RADIUS)
      ETAMEAN=0.25*(ETA(I,J+1)+ETA(I+1,J+1)+ETA(I,J)+ETA(I+1,J))
      ZETAMEAN=0.25*(ZETA(I,J+1)+ZETA(I+1,J+1)+ZETA(I,J)+ZETA(I+1,J))
      AA1=ETA(I,J+1)+ZETA(I,J+1)+ETA(I+1,J+1)+ZETA(I+1,J+1)
      AA2=ETA(I,J)+ZETA(I,J)+ETA(I+1,J)+ZETA(I+1,J)
      AA3=(ETA(I+1,J)/CSU(J)+ETA(I+1,J+1)/CSU(J))/CSU(J)
      AA4=(ETA(I,J)/CSU(J)+ETA(I,J+1)/CSU(J))/CSU(J)
      AA5=((ZETA(I,J+1)-ETA(I,J+1))+(ZETA(I+1,J+1)-ETA(I+1,J+1))
     &-(ZETA(I,J)-ETA(I,J))-(ZETA(I+1,J)-ETA(I+1,J)))*TNG(J)
      AA6=2.0*ETAMEAN*TNG(J)*TNG(J)

      AV(I,J)=(-AA2*DELY2-(ZETAMEAN-ETAMEAN)*TNG(J-1)*DELYR
     &-ETAMEAN*2.0*TNG(J)*DELYR)*UVM(I,J)
      BV(I,J)=((AA1+AA2)*DELY2+AA5*DELYR+AA6*BYRAD2
     &+AMASS(I,J)*BYDTS*2.0+DRAGS(I,J))*UVM(I,J)+(1.0-UVM(I,J))
      CV(I,J)=(-AA1*DELY2+(ZETAMEAN-ETAMEAN)*TNG(J+1)*DELYR
     &+ETAMEAN*2.0*TNG(J)*DELYR)*UVM(I,J)
      END DO
      END DO
      DO I=2,NXLCYC
      AV(I,2)=0.0
      CV(I,NYPOLE)=0.0
      CV(I,2)=CV(I,2)/BV(I,2)
      END DO

      DO 1301 I=2,NXLCYC
      DO J=2,NYPOLE
      DELX2=BYDX2(I)       ! 0.5/(DXU(I)*DXU(I))
      DELY2=BYDY2(J)       ! 0.5/(DYU(J)*DYU(J))
      DELYR=BYDYR(J)       ! 0.5/(DYU(J)*RADIUS)

      ETAMEAN=0.25*(ETA(I,J+1)+ETA(I+1,J+1)+ETA(I,J)+ETA(I+1,J))
      ZETAMEAN=0.25*(ZETA(I,J+1)+ZETA(I+1,J+1)+ZETA(I,J)+ZETA(I+1,J))

      AA1=ETA(I,J+1)+ZETA(I,J+1)+ETA(I+1,J+1)+ZETA(I+1,J+1)
      AA2=ETA(I,J)+ZETA(I,J)+ETA(I+1,J)+ZETA(I+1,J)
      AA3=(ETA(I+1,J)/CSU(J)+ETA(I+1,J+1)/CSU(J))/CSU(J)
      AA4=(ETA(I,J)/CSU(J)+ETA(I,J+1)/CSU(J))/CSU(J)
      AA5=((ZETA(I,J+1)-ETA(I,J+1))+(ZETA(I+1,J+1)-ETA(I+1,J+1))
     &-(ZETA(I,J)-ETA(I,J))-(ZETA(I+1,J)-ETA(I+1,J)))*TNG(J)
      AA6=2.0*ETAMEAN*TNG(J)*TNG(J)

      IF(J.EQ.NYPOLE) THEN
      AA9=(AA1*DELY2-(ZETAMEAN-ETAMEAN)*TNG(J+1)*DELYR
     &-ETAMEAN*2.0*TNG(J)*DELYR)*VICEC(I,J+1)*UVM(I,J)*FLOAT(NPOL-0)
      ELSE
      AA9=0.0
      END IF

      VRT(J)=AA9+FXY(I,J)-(AA3+AA4)*DELX2*VICE(I,J,2)
     6+((ETA(I+1,J)/CSU(J)+ETA(I+1,J+1)/CSU(J))*VICE(I+1,J,2)*DELX2
     7+(ETA(I,J)/CSU(J)+ETA(I,J+1)/CSU(J))*VICE(I-1,J,2)*DELX2)/CSU(J)
      VRT(J)=(VRT(J)+AMASS(I,J)*BYDTS*VICE(I,J,2)*2.0)*UVM(I,J)
      END DO

      DO J=2,NYPOLE
      CVV(J)=CV(I,J)
      END DO
      VRT(2)=VRT(2)/BV(I,2)
      DO J=3,NYPOLE
      JMD=J-1
      CVV(J)=CVV(J)/(BV(I,J)-AV(I,J)*CVV(JMD))
      VRT(J)=(VRT(J)-AV(I,J)*VRT(JMD))/(BV(I,J)-AV(I,J)*CVV(JMD))
      END DO
      DO J=1,NYPOLE-2
      J1=NYPOLE-J
      J2=J1+1
      VRT(J1)=VRT(J1)-CVV(J1)*VRT(J2)
      END DO
      DO J=2,NYPOLE
      VICE(I,J,1)=VRT(J)
      END DO
 1301 CONTINUE

      DO I=2,NXLCYC
      DO J=2,NYPOLE
      VICE(I,J,3)=VICE(I,J,1)
      END DO
      END DO

C NOW THE SECOND HALF

      DO J=2,NYPOLE
      DO I=2,NXLCYC
      DELX2=BYDX2(I)       ! 0.5/(DXU(I)*DXU(I))
      DELY2=BYDY2(J)       ! 0.5/(DYU(J)*DYU(J))
      DELYR=BYDYR(J)       ! 0.5/(DYU(J)*RADIUS)
      ETAMEAN=0.25*(ETA(I,J+1)+ETA(I+1,J+1)+ETA(I,J)+ETA(I+1,J))
      ZETAMEAN=0.25*(ZETA(I,J+1)+ZETA(I+1,J+1)+ZETA(I,J)+ZETA(I+1,J))

      AA1=ETA(I,J+1)+ZETA(I,J+1)+ETA(I+1,J+1)+ZETA(I+1,J+1)
      AA2=ETA(I,J)+ZETA(I,J)+ETA(I+1,J)+ZETA(I+1,J)
      AA3=(ETA(I+1,J)/CSU(J)+ETA(I+1,J+1)/CSU(J))/CSU(J)
      AA4=(ETA(I,J)/CSU(J)+ETA(I,J+1)/CSU(J))/CSU(J)
      AA6=2.0*ETAMEAN*TNG(J)*TNG(J)

      AU(I,J)=-AA4*DELX2*UVM(I,J)
      BU(I,J)=((AA3+AA4)*DELX2+AA6*BYRAD2
     &+AMASS(I,J)*BYDTS*2.0+DRAGS(I,J))*UVM(I,J)+(1.0-UVM(I,J))
      CU(I,J)=-AA3*DELX2*UVM(I,J)
      END DO
      END DO
      DO J=2,NYPOLE
      AU(2,J)=0.0
      CU(NXLCYC,J)=0.0
      CU(2,J)=CU(2,J)/BU(2,J)
      END DO

      DO J=2,NYPOLE
      DO I=2,NXLCYC

      DELX2=BYDX2(I)       ! 0.5/(DXU(I)*DXU(I))
      DELY2=BYDY2(J)       ! 0.5/(DYU(J)*DYU(J))
      DELYR=BYDYR(J)       ! 0.5/(DYU(J)*RADIUS)
      ETAMEAN=0.25*(ETA(I,J+1)+ETA(I+1,J+1)+ETA(I,J)+ETA(I+1,J))
      ZETAMEAN=0.25*(ZETA(I,J+1)+ZETA(I+1,J+1)+ZETA(I,J)+ZETA(I+1,J))

      AA1=ETA(I,J+1)+ZETA(I,J+1)+ETA(I+1,J+1)+ZETA(I+1,J+1)
      AA2=ETA(I,J)+ZETA(I,J)+ETA(I+1,J)+ZETA(I+1,J)
      AA3=(ETA(I+1,J)/CSU(J)+ETA(I+1,J+1)/CSU(J))/CSU(J)
      AA4=(ETA(I,J)/CSU(J)+ETA(I,J+1)/CSU(J))/CSU(J)
      AA5=((ZETA(I,J+1)-ETA(I,J+1))+(ZETA(I+1,J+1)-ETA(I+1,J+1))
     &-(ZETA(I,J)-ETA(I,J))-(ZETA(I+1,J)-ETA(I+1,J)))*TNG(J)
      AA6=2.0*ETAMEAN*TNG(J)*TNG(J)

      IF(I.EQ.2) THEN
      AA9=AA4*DELX2*VICEC(I-1,J)*UVM(I,J)*FLOAT(LCYC-0)
      ELSE IF(I.EQ.NXLCYC) THEN
      AA9=AA3*DELX2*VICEC(I+1,J)*UVM(I,J)*FLOAT(LCYC-0)
      ELSE
      AA9=0.0
      END IF
      FXY1(I,J)=AA9+AMASS(I,J)*BYDTS*VICE(I,J,1)*2.0
     1-AA5*DELYR*VICE(I,J,1)
     1-(AA1+AA2)*DELY2*VICE(I,J,1)
     1+AA1*DELY2*VICE(I,J+1,1)-((ZETAMEAN-ETAMEAN)*TNG(J+1)*DELYR
     1+ETAMEAN*2.0*TNG(J)*DELYR)*VICE(I,J+1,1)
     2+AA2*DELY2*VICE(I,J-1,1)+((ZETAMEAN-ETAMEAN)*TNG(J-1)*DELYR
     2+ETAMEAN*2.0*TNG(J)*DELYR)*VICE(I,J-1,1)
      END DO
      END DO

      DO 1201 J=2,NYPOLE
      DO I=2,NXLCYC
      URT(I)=FXY(I,J)+FXY1(I,J)
      URT(I)=URT(I)*UVM(I,J)
      END DO

      DO I=2,NXLCYC
      CUU(I)=CU(I,J)
      END DO
      URT(2)=URT(2)/BU(2,J)
      DO I=3,NXLCYC
      IMD=I-1
      CUU(I)=CUU(I)/(BU(I,J)-AU(I,J)*CUU(IMD))
      URT(I)=(URT(I)-AU(I,J)*URT(IMD))/(BU(I,J)-AU(I,J)*CUU(IMD))
      END DO
      DO I=1,NXLCYC-2
      J1=NXLCYC-I
      J2=J1+1
      URT(J1)=URT(J1)-CUU(J1)*URT(J2)
      END DO
      DO I=2,NXLCYC
      VICE(I,J,1)=URT(I)
      END DO
 1201 CONTINUE

      DO J=2,NYPOLE
      DO I=2,NXLCYC
      UICE(I,J,1)=UICE(I,J,1)*UVM(I,J)
      VICE(I,J,1)=VICE(I,J,1)*UVM(I,J)
      END DO
      END DO

      RETURN
      END SUBROUTINE RELAX

      SUBROUTINE ADVSI
!@sum  ADVSI advects sea ice
!@auth Gary Russell
!@ver  1.0
      USE CONSTANT, only : byshi,lhm
      USE MODEL_COM, only : im,jm,focean,itocean,itoice,ftype
      USE OCEAN, only : dxyp=>dxypo,dyp=>dypo,dxp=>dxpo,dxv=>dxvo,
     *     bydxyp=>bydxypo
      USE ICEDYN, only : usidt,vsidt,rsix,rsiy,rsisave
      USE SEAICE, only : ace1i,xsi_glob=>xsi
      USE SEAICE_COM, only : rsi,msi,snowi,hsi,ssi,lmi
#ifdef TRACERS_WATER
     *     ,trsi,ntm
#endif
      USE ODIAG, only : oij,ij_musi,ij_mvsi
      USE FLUXES, only : gtemp
#ifdef TRACERS_WATER
     *     ,gtracer
#endif
      IMPLICIT NONE
      REAL*8, DIMENSION(IM) :: FAW,FASI,FXSI,FYSI
!@var NTRICE max. number of tracers to be advected (mass/heat/salt+)
#ifndef TRACERS_WATER
      INTEGER, PARAMETER :: NTRICE=2+2*LMI
#else
      INTEGER, PARAMETER :: NTRICE=2+(2+NTM)*LMI
      INTEGER ITR
#endif
      REAL*8 FMSI(NTRICE,IM),SFMSI(NTRICE),AMSI(NTRICE)
      INTEGER I,J,L,IM1,IP1,K
      REAL*8 SFASI,DMHSI,ASI,YSI,XSI,FRSI
!@var MHS mass/heat/salt content of sea ice
      REAL*8, DIMENSION(NTRICE,IM,JM) :: MHS
C****
C**** FLUXCB  USIDT  U compon of time integrated sea ice velocity (m)
C****         VSIDT  V compon of time integrated sea ice velocity (m)
C****
C**** WORK01  FAW    flux of surface water area (m^2) = USIDT*DYP
C****         FASI   flux of sea ice area (m^2) = USIDT*DYP*RSIedge
C****         FMSI   flux of sea ice mass (kg) or heat (J) or salt (kg)

C**** Regularise ice concentration gradients to prevent advection errors
      DO J=2,JM-1
      DO I=1,IM
        IF (RSI(I,J).gt.1d-4) THEN
          IF (RSISAVE(I,J).gt.RSI(I,J)) THEN ! reduce gradients
            FRSI=(RSISAVE(I,J)-RSI(I,J))/RSISAVE(I,J)
            RSIX(I,J)=RSIX(I,J)*(1.-FRSI)
            RSIY(I,J)=RSIY(I,J)*(1.-FRSI)
          END IF
          IF(RSI(I,J)-RSIX(I,J).lt.0.)  RSIX(I,J) =    RSI(I,J)
          IF(RSI(I,J)+RSIX(I,J).lt.0.)  RSIX(I,J) =   -RSI(I,J)
          IF(RSI(I,J)-RSIX(I,J).gt.1d0) RSIX(I,J) =    RSI(I,J)-1d0
          IF(RSI(I,J)+RSIX(I,J).gt.1d0) RSIX(I,J) =1d0-RSI(I,J)
          IF(RSI(I,J)-RSIY(I,J).lt.0.)  RSIY(I,J) =    RSI(I,J)
          IF(RSI(I,J)+RSIY(I,J).lt.0.)  RSIY(I,J) =   -RSI(I,J)
          IF(RSI(I,J)-RSIY(I,J).gt.1d0) RSIY(I,J) =    RSI(I,J)-1d0
          IF(RSI(I,J)+RSIY(I,J).gt.1d0) RSIY(I,J) =1d0-RSI(I,J)
        ELSE
          RSIX(I,J) = 0.  ; RSIY(I,J) = 0.
        END IF
      END DO
      END DO

C**** set up local MHS array to contain all advected quantities
C**** MHS(1:2) = MASS, MHS(3:6) = HEAT, MHS(7:10)=SALT
      MHS(1,:,:) = ACE1I + SNOWI
      MHS(2,:,:) = MSI
      DO L=1,LMI
        MHS(L+2,:,:) = HSI(L,:,:)
        MHS(L+2+LMI,:,:) = SSI(L,:,:)
#ifdef TRACERS_WATER
C**** add tracers to advected arrays
        DO ITR=1,NTM
          MHS(L+2+(1+ITR)*LMI,:,:)=TRSI(ITR,L,:,:)
        END DO
#endif
      END DO
C****
C**** North-South Advection of Sea Ice
C****
      SFASI  = 0.
      DO 5 K=1,NTRICE
    5 SFMSI(K) = 0.
      DO 340 I=1,IM
C****
C**** Calculate south-north sea ice fluxes at grid box edges
C****
      DO 120 J=2,JM-2
      IF(VSIDT(I,J).eq.0.)  GO TO 120
      FAW(J) = VSIDT(I,J)*DXV(J)
      IF(VSIDT(I,J).gt.0.)  GO TO 110
C**** Sea ice velocity is southward at grid box edge
      FASI(J) =FAW(J)*(RSI(I,J+1)-(1d0+FAW(J)*BYDXYP(J+1))*RSIY(I,J+1))
      FXSI(J) = FAW(J)*RSIX(I,J+1)
      FYSI(J) = FAW(J)*
     *         (FAW(J)*BYDXYP(J+1)*FAW(J)*RSIY(I,J+1) - 3d0*FASI(J))
      DO 105 K=1,NTRICE
  105 FMSI(K,J) = FASI(J)*MHS(K,I,J+1)
         OIJ(I,J,IJ_MVSI) = OIJ(I,J,IJ_MVSI) + (FMSI(1,J)+FMSI(2,J))
      GO TO 120
C**** Sea ice velocity is northward at grid box edge
  110 FASI(J) = FAW(J)*(RSI(I,J)+(1d0-FAW(J)*BYDXYP(J))*RSIY(I,J))
      FXSI(J) = FAW(J)*RSIX(I,J)
      FYSI(J) = FAW(J)*(FAW(J)*BYDXYP(J)*FAW(J)*RSIY(I,J)-3d0*FASI(J))
      DO 115 K=1,NTRICE
  115 FMSI(K,J) = FASI(J)*MHS(K,I,J)
         OIJ(I,J,IJ_MVSI) = OIJ(I,J,IJ_MVSI) + (FMSI(1,J)+FMSI(2,J))
  120 CONTINUE
C****
C**** Calculate south-north sea ice fluxes near North Pole
C****
      IF(VSIDT(I,JM-1).eq.0.)  GO TO 200
      FAW(JM-1) = VSIDT(I,JM-1)*DXV(JM-1)
      IF(VSIDT(I,JM-1).gt.0.)  GO TO 130
C**** Sea ice velocity is southward from North Pole box
      FASI(JM-1) = FAW(JM-1)*RSI(1,JM)
      FXSI(JM-1) = 0.
      FYSI(JM-1) = -FAW(JM-1)*FASI(JM-1)
      DO 125 K=1,NTRICE
  125 FMSI(K,JM-1) = FASI(JM-1)*MHS(K,1,JM)
         OIJ(I,JM-1,IJ_MVSI) = OIJ(I,JM-1,IJ_MVSI) +
     *     (FMSI(1,JM-1)+FMSI(2,JM-1))
      GO TO 140
C**** Sea ice velocity is northward into North Pole box
  130 FASI(JM-1) = FAW(JM-1)*
     *  (RSI(I,JM-1)+(1d0-FAW(JM-1)*BYDXYP(JM-1))*RSIY(I,JM-1))
      FXSI(JM-1) = FAW(JM-1)*RSIX(I,JM-1)
      FYSI(JM-1) = FAW(JM-1)*
     *  (FAW(JM-1)*BYDXYP(JM-1)*FAW(JM-1)*RSIY(I,JM-1)-3d0*FASI(JM-1))
      DO 135 K=1,NTRICE
  135 FMSI(K,JM-1) = FASI(JM-1)*MHS(K,I,JM-1)
         OIJ(I,JM-1,IJ_MVSI) = OIJ(I,JM-1,IJ_MVSI) +
     *     (FMSI(1,JM-1)+FMSI(2,JM-1))
C**** Accumulate sea ice leaving and entering North Pole box
  140 SFASI = SFASI + FASI(JM-1)
      DO 145 K=1,NTRICE
  145 SFMSI(K) = SFMSI(K) + FMSI(K,JM-1)
C****
C**** Update sea ice variables due to south-north fluxes
C****
  200 DO 330 J=2,JM-1
      IF(VSIDT(I,J-1)) 240,210,280
C**** VSIDT(J-1)=0.
  210 IF(VSIDT(I,J))  220,330,230
C**** VSIDT(J-1)=0, VSIDT(J)<0.
  220 ASI = RSI(I,J)*DXYP(J) -  FASI(J)
      DO 225 K=1,NTRICE
  225 AMSI(K) = RSI(I,J)*DXYP(J)*MHS(K,I,J) - FMSI(K,J)
      IF(ASI.gt.DXYP(J))  GO TO 320
      YSI = (RSIY(I,J)*DXYP(J)*DXYP(J) - FYSI(J)
     *    + 3d0*(FAW(J)*ASI-DXYP(J)*FASI(J))) / (DXYP(J)-FAW(J))
      RSI(I,J)  = ASI*BYDXYP(J)
      RSIY(I,J) = YSI*BYDXYP(J)
      RSIX(I,J) = RSIX(I,J) - FXSI(J)*BYDXYP(J)
      DO 226 K=1,NTRICE
  226 MHS(K,I,J) = AMSI(K)/ASI
      GO TO 310
C**** VSIDT(J-1)=0, VSIDT(J)>0.
  230 RSI(I,J)  =  RSI(I,J) -  FASI(J)*BYDXYP(J)
      RSIX(I,J) = RSIX(I,J)*(1d0-FAW(J)*BYDXYP(J))
      RSIY(I,J) = RSIY(I,J)*(1d0-FAW(J)*BYDXYP(J))**2
      GO TO 310
C**** VSIDT(J-1)<0.
  240 IF(VSIDT(I,J))  260,250,270
C**** VSIDT(J-1)<0, VSIDT(J)=0.
  250 RSI(I,J)  =  RSI(I,J) +  FASI(J-1)*BYDXYP(J)
      RSIX(I,J) = RSIX(I,J)*(1d0+FAW(J-1)*BYDXYP(J))
      RSIY(I,J) = RSIY(I,J)*(1d0+FAW(J-1)*BYDXYP(J))**2
      GO TO 310
C**** VSIDT(J-1)<0, VSIDT(J)<0  or  VSIDT(J-1)>0, VSIDT(J)?0.
  260 ASI = RSI(I,J)*DXYP(J) + ( FASI(J-1)- FASI(J))
      DO 265 K=1,NTRICE
  265 AMSI(K) = RSI(I,J)*DXYP(J)*MHS(K,I,J) + (FMSI(K,J-1)-FMSI(K,J))
      IF(ASI.gt.DXYP(J))  GO TO 320
      YSI = (RSIY(I,J)*DXYP(J)*DXYP(J) + (FYSI(J-1)-FYSI(J))
     *    + 3d0*((FAW(J-1)+FAW(J))*ASI-DXYP(J)*(FASI(J-1)+FASI(J))))
     *    / (DXYP(J) + (FAW(J-1)-FAW(J)))
      RSI(I,J)  = ASI*BYDXYP(J)
      RSIY(I,J) = YSI*BYDXYP(J)
      RSIX(I,J) = RSIX(I,J) + (FXSI(J-1)-FXSI(J))*BYDXYP(J)
      DO 266 K=1,NTRICE
  266 MHS(K,I,J) = AMSI(K)/ASI
      GO TO 310
C**** VSIDT(J-1)<0, VSIDT(J)>0.
  270 RSI(I,J)  =  RSI(I,J) + (FASI(J-1)-FASI(J))*BYDXYP(J)
      RSIX(I,J) = RSIX(I,J)*(1d0+(FAW(J-1)-FAW(J))*BYDXYP(J))
      RSIY(I,J) = RSIY(I,J)*(1d0+(FAW(J-1)-FAW(J))*BYDXYP(J))**2
      GO TO 310
C**** VSIDT(J-1)>0.
  280 IF(VSIDT(I,J).ne.0.)  GO TO 260
C**** VSIDT(J-1)>0, VSIDT(J)=0.
      ASI = RSI(I,J)*DXYP(J) + FASI(J-1)
      DO 285 K=1,NTRICE
  285 AMSI(K) = RSI(I,J)*DXYP(J)*MHS(K,I,J) + FMSI(K,J-1)
      IF(ASI.gt.DXYP(J))  GO TO 320
      YSI = (RSIY(I,J)*DXYP(J)*DXYP(J) + FYSI(J-1)
     *    + 3d0*(FAW(J-1)*ASI-DXYP(J)*FASI(J-1))) / (DXYP(J)+FAW(J-1))
      RSI(I,J)  = ASI*BYDXYP(J)
      RSIY(I,J) = YSI*BYDXYP(J)
      RSIX(I,J) = RSIX(I,J) + FXSI(J-1)*BYDXYP(J)
      DO 286 K=1,NTRICE
  286 MHS(K,I,J) = AMSI(K)/ASI
C**** Limit RSIX and RSIY so that sea ice is positive at the edges
  310 RSI(I,J) = MAX(0d0,RSI(I,J))
      IF(RSI(I,J)-RSIX(I,J).lt.0.)  RSIX(I,J) =    RSI(I,J)
      IF(RSI(I,J)+RSIX(I,J).lt.0.)  RSIX(I,J) =   -RSI(I,J)
      IF(RSI(I,J)-RSIX(I,J).gt.1d0) RSIX(I,J) =    RSI(I,J)-1d0
      IF(RSI(I,J)+RSIX(I,J).gt.1d0) RSIX(I,J) =1d0-RSI(I,J)
      IF(RSI(I,J)-RSIY(I,J).lt.0.)  RSIY(I,J) =    RSI(I,J)
      IF(RSI(I,J)+RSIY(I,J).lt.0.)  RSIY(I,J) =   -RSI(I,J)
      IF(RSI(I,J)-RSIY(I,J).gt.1d0) RSIY(I,J) =    RSI(I,J)-1d0
      IF(RSI(I,J)+RSIY(I,J).gt.1d0) RSIY(I,J) =1d0-RSI(I,J)
      GO TO 330
C**** Sea ice crunches into itself and completely covers grid box
  320 RSI(I,J)   = 1d0
      RSIX(I,J)  = 0.
      RSIY(I,J)  = 0.
      MHS(1,I,J) = AMSI(1)/ASI
      MHS(2,I,J) =(AMSI(1)+AMSI(2))*BYDXYP(J) - MHS(1,I,J)
      DO K=1,(NTRICE-2)/LMI
        MHS(3+LMI*(K-1),I,J) = AMSI(3+LMI*(K-1)) / ASI
        MHS(4+LMI*(K-1),I,J) = AMSI(4+LMI*(K-1)) / ASI
        DMHSI = (AMSI(3+LMI*(K-1))+AMSI(4+LMI*(K-1))+AMSI(5+LMI*(K-1))
     *       +AMSI(6+LMI*(K-1)))*(BYDXYP(J) -1d0 / ASI )
        MHS(5+LMI*(K-1),I,J) = AMSI(5+LMI*(K-1)) / ASI +
     *       XSI_GLOB(3)*DMHSI
        MHS(6+LMI*(K-1),I,J) = AMSI(6+LMI*(K-1)) / ASI +
     *       XSI_GLOB(4)*DMHSI
      END DO
C**** End of loop over J
  330 CONTINUE
C**** End of loop over I
  340 CONTINUE
C****
C**** Advection of Sea Ice leaving and entering North Pole box
C****
      ASI = RSI(1,JM)*DXYP(JM) + SFASI/IM
      DO 345 K=1,NTRICE
  345 AMSI(K) = RSI(1,JM)*DXYP(JM)*MHS(K,1,JM) + SFMSI(K)/IM
      IF(ASI.gt.DXYP(JM))  GO TO 350
      RSI(1,JM)   = ASI*BYDXYP(JM)
      DO 346 K=1,NTRICE
  346 MHS(K,1,JM) = AMSI(K)/ASI
      GO TO 400
C**** Sea ice crunches into itself at North Pole box
  350 RSI(1,JM)   = 1d0
      MHS(1,1,JM) = AMSI(1)/ASI
      MHS(2,1,JM) =(AMSI(1)+AMSI(2))*BYDXYP(JM) - MHS(1,1,JM)
      DO K=1,(NTRICE-2)/LMI
        MHS(3+LMI*(K-1),1,JM) = AMSI(3+LMI*(K-1)) / ASI
        MHS(4+LMI*(K-1),1,JM) = AMSI(4+LMI*(K-1)) / ASI
        DMHSI = (AMSI(3+LMI*(K-1))+AMSI(4+LMI*(K-1))+AMSI(5+LMI*(K-1))
     *       +AMSI(6+LMI*(K-1)))*(BYDXYP(JM) -1d0/ ASI)
        MHS(5+LMI*(K-1),1,JM) = AMSI(5+LMI*(K-1)) / ASI +
     *       XSI_GLOB(3)*DMHSI
        MHS(6+LMI*(K-1),1,JM) = AMSI(6+LMI*(K-1)) / ASI +
     *       XSI_GLOB(4)*DMHSI
      END DO
C****
C**** East-West Advection of Sea Ice
C****
  400 DO 640 J=2,JM-1
C****
C**** Calculate west-east sea ice fluxes at grid box edges
C****
      I=IM
      DO 420 IP1=1,IM
      IF(USIDT(I,J).eq.0.)  GO TO 420
      FAW(I) = USIDT(I,J)*DYP(J)
      IF(USIDT(I,J).gt.0.)  GO TO 410
C**** Sea ice velocity is westward at grid box edge
      FASI(I) =FAW(I)*(RSI(IP1,J)-(1d0+FAW(I)*BYDXYP(J))*RSIX(IP1,J))
      FXSI(I) =FAW(I)*(FAW(I)*BYDXYP(J)*FAW(I)*RSIX(IP1,J)-3d0*FASI(I))
      FYSI(I) = FAW(I)*RSIY(IP1,J)
      DO 405 K=1,NTRICE
  405 FMSI(K,I) = FASI(I)*MHS(K,IP1,J)
         OIJ(I,J,IJ_MUSI) = OIJ(I,J,IJ_MUSI) + (FMSI(1,I)+FMSI(2,I))
      GO TO 420
C**** Sea ice velocity is eastward at grid box edge
  410 FASI(I) = FAW(I)*(RSI(I,J)+(1d0-FAW(I)*BYDXYP(J))*RSIX(I,J))
      FXSI(I) = FAW(I)*(FAW(I)*BYDXYP(J)*FAW(I)*RSIX(I,J)-3d0*FASI(I))
      FYSI(I) = FAW(I)*RSIY(I,J)
      DO 415 K=1,NTRICE
  415 FMSI(K,I) = FASI(I)*MHS(K,I,J)
         OIJ(I,J,IJ_MUSI) = OIJ(I,J,IJ_MUSI) + (FMSI(1,I)+FMSI(2,I))
  420 I=IP1
C****
C**** Update sea ice variables due to west-east fluxes
C****
      IM1=IM
      DO 630 I=1,IM
      IF(USIDT(IM1,J)) 540,510,580
C**** USIDT(IM1)=0.
  510 IF(USIDT(I,J))  520,630,530
C**** USIDT(IM1)=0, USIDT(I)<0.
  520 ASI = RSI(I,J)*DXYP(J) -  FASI(I)
      DO 525 K=1,NTRICE
  525 AMSI(K) = RSI(I,J)*DXYP(J)*MHS(K,I,J) - FMSI(K,I)
      IF(ASI.gt.DXYP(J))  GO TO 620
      XSI = (RSIX(I,J)*DXYP(J)*DXYP(J) - FXSI(I)
     *    + 3d0*(FAW(I)*ASI-DXYP(J)*FASI(I))) / (DXYP(J)-FAW(I))
      RSI(I,J)  = ASI*BYDXYP(J)
      RSIX(I,J) = XSI*BYDXYP(J)
      RSIY(I,J) = RSIY(I,J) - FYSI(I)*BYDXYP(J)
      DO 526 K=1,NTRICE
  526 MHS(K,I,J) = AMSI(K)/ASI
      GO TO 610
C**** USIDT(IM1)=0, USIDT(I)>0.
  530 RSI(I,J)  =  RSI(I,J) -  FASI(I)*BYDXYP(J)
      RSIX(I,J) = RSIX(I,J)*(1d0-FAW(I)*BYDXYP(J))**2
      RSIY(I,J) = RSIY(I,J)*(1d0-FAW(I)*BYDXYP(J))
      GO TO 610
C**** USIDT(IM1)<0.
  540 IF(USIDT(I,J))  560,550,570
C**** USIDT(IM1)<0, USIDT(I)=0.
  550 RSI(I,J)  =  RSI(I,J) +  FASI(IM1)*BYDXYP(J)
      RSIX(I,J) = RSIX(I,J)*(1d0+FAW(IM1)*BYDXYP(J))**2
      RSIY(I,J) = RSIY(I,J)*(1d0+FAW(IM1)*BYDXYP(J))
      GO TO 610
C**** USIDT(IM1)<0, USIDT(I)<0  or  USIDT(IM1)>0, USIDT(I)?0.
  560 ASI = RSI(I,J)*DXYP(J) + (FASI(IM1)- FASI(I))
      DO 565 K=1,NTRICE
  565 AMSI(K) = RSI(I,J)*DXYP(J)*MHS(K,I,J) + (FMSI(K,IM1)-FMSI(K,I))
      IF(ASI.gt.DXYP(J))  GO TO 620
      XSI = (RSIX(I,J)*DXYP(J)*DXYP(J) + (FXSI(IM1)-FXSI(I))
     *    + 3d0*((FAW(IM1)+FAW(I))*ASI-DXYP(J)*(FASI(IM1)+FASI(I))))
     *    / (DXYP(J) + (FAW(IM1)-FAW(I)))
      RSI(I,J)  = ASI*BYDXYP(J)
      RSIX(I,J) = XSI*BYDXYP(J)
      RSIY(I,J) = RSIY(I,J) + (FYSI(IM1)-FYSI(I))*BYDXYP(J)
      DO 566 K=1,NTRICE
  566 MHS(K,I,J) = AMSI(K)/ASI
      GO TO 610
C**** USIDT(IM1)<0, USIDT(I)>0.
  570 RSI(I,J)  =  RSI(I,J) + (FASI(IM1)-FASI(I))*BYDXYP(J)
      RSIX(I,J) = RSIX(I,J)*(1d0+(FAW(IM1)-FAW(I))*BYDXYP(J))**2
      RSIY(I,J) = RSIY(I,J)*(1d0+(FAW(IM1)-FAW(I))*BYDXYP(J))
      GO TO 610
C**** USIDT(IM1)>0.
  580 IF(USIDT(I,J).ne.0.)  GO TO 560
C**** USIDT(IM1)>0, USIDT(I)=0.
      ASI = RSI(I,J)*DXYP(J) + FASI(IM1)
      DO 585 K=1,NTRICE
  585 AMSI(K) = RSI(I,J)*DXYP(J)*MHS(K,I,J) + FMSI(K,IM1)
      IF(ASI.gt.DXYP(J))  GO TO 620
      XSI = (RSIX(I,J)*DXYP(J)*DXYP(J) + FXSI(IM1)
     *    + 3d0*(FAW(IM1)*ASI-DXYP(J)*FASI(IM1))) / (DXYP(J)+FAW(IM1))
      RSI(I,J)  = ASI*BYDXYP(J)
      RSIX(I,J) = XSI*BYDXYP(J)
      RSIY(I,J) = RSIY(I,J) + FYSI(IM1)*BYDXYP(J)
      DO 586 K=1,NTRICE
  586 MHS(K,I,J) = AMSI(K)/ASI
C**** Limit RSIX and RSIY so that sea ice is positive at the edges
  610 RSI(I,J) = MAX(0d0,RSI(I,J))
      IF(RSI(I,J)-RSIX(I,J).lt.0.)  RSIX(I,J) =    RSI(I,J)
      IF(RSI(I,J)+RSIX(I,J).lt.0.)  RSIX(I,J) =   -RSI(I,J)
      IF(RSI(I,J)-RSIX(I,J).gt.1d0) RSIX(I,J) =    RSI(I,J)-1d0
      IF(RSI(I,J)+RSIX(I,J).gt.1d0) RSIX(I,J) =1d0-RSI(I,J)
      IF(RSI(I,J)-RSIY(I,J).lt.0.)  RSIY(I,J) =    RSI(I,J)
      IF(RSI(I,J)+RSIY(I,J).lt.0.)  RSIY(I,J) =   -RSI(I,J)
      IF(RSI(I,J)-RSIY(I,J).gt.1d0) RSIY(I,J) =    RSI(I,J)-1d0
      IF(RSI(I,J)+RSIY(I,J).gt.1d0) RSIY(I,J) =1d0-RSI(I,J)
      GO TO 630
C**** Sea ice crunches into itself and completely covers grid box
  620 RSI(I,J)   = 1d0
      RSIX(I,J)  = 0.
      RSIY(I,J)  = 0.
      MHS(1,I,J) = AMSI(1)/ASI
      MHS(2,I,J) =(AMSI(1)+AMSI(2))*BYDXYP(J) - MHS(1,I,J)
      DO K=1,(NTRICE-2)/LMI
        MHS(3+LMI*(K-1),I,J) = AMSI(3+LMI*(K-1)) / ASI
        MHS(4+LMI*(K-1),I,J) = AMSI(4+LMI*(K-1)) / ASI
        DMHSI = (AMSI(3+LMI*(K-1))+AMSI(4+LMI*(K-1))+AMSI(5+LMI*(K-1))
     *       +AMSI(6+LMI*(K-1)))*(BYDXYP(J) -1d0/ ASI)
        MHS(5+LMI*(K-1),I,J) = AMSI(5+LMI*(K-1)) / ASI +
     *       XSI_GLOB(3)*DMHSI
        MHS(6+LMI*(K-1),I,J) = AMSI(6+LMI*(K-1)) / ASI +
     *       XSI_GLOB(4)*DMHSI
      END DO
C**** End of loop over I
  630 IM1=I
C**** End of loop over J
  640 CONTINUE

C**** set global variables from local array
      SNOWI(:,:)= MAX(0d0,MHS(1,:,:) - ACE1I)
      MSI(:,:)  = MHS(2,:,:)
      DO L=1,LMI
        HSI(L,:,:) = MHS(L+2,:,:)
        SSI(L,:,:) = MHS(L+2+LMI,:,:)
#ifdef TRACERS_WATER
        DO ITR=1,NTM
          TRSI(ITR,L,:,:)=MHS(L+2+(1+ITR)*LMI,:,:)
        END DO
#endif
      END DO
C**** Set atmospheric arrays
      DO J=1,JM
        DO I=1,IM
          IF (FOCEAN(I,J).gt.0) THEN
            FTYPE(ITOICE ,I,J)=FOCEAN(I,J)*RSI(I,J)
            FTYPE(ITOCEAN,I,J)=FOCEAN(I,J) - FTYPE(ITOICE ,I,J)
            GTEMP(1:2,2,I,J)=(HSI(1:2,I,J)/(XSI_GLOB(1:2)*(ACE1I
     *           +SNOWI(I,J)))+LHM)*BYSHI
#ifdef TRACERS_WATER
            GTRACER(:,1,I,J)=TRSI(:,1,I,J)/(MHS(1,I,J)-SSI(1,I,J))
#endif
          END IF
        END DO
      END DO
C****
      RETURN
      END

      SUBROUTINE AT2IT(FIELDA,FIELDI,NF,QCONSERV)
!@sum  AT2IT interpolates Atm Tracer grid to Dynamic Ice Tracer grid
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : ima=>im,jma=>jm
      USE ICEDYN, only : nx1,ny1,ratic
      IMPLICIT NONE
!@var QCONSERV true if integrated field must be conserved
      LOGICAL, INTENT(IN) :: QCONSERV
!@var N number of fields
      INTEGER, INTENT(IN) :: NF
!@var FIELDA array on atmospheric tracer grid
      REAL*8, INTENT(IN), DIMENSION(NF,IMA,JMA) :: FIELDA
!@var FIELDI array on dynamic ice tracer grid
      REAL*8, INTENT(OUT), DIMENSION(NF,NX1,NY1) :: FIELDI
      INTEGER I,J

C**** currently no need for interpolation,
C**** just scaling due to area differences for fluxes
C**** and shift one box in longitude
      IF (QCONSERV) THEN
        DO J=1,NY1
          DO I=2,NX1-1
            FIELDI(:,I,J) = FIELDA(:,I-1,J)*RATIC(J)
          END DO
        END DO
      ELSE
        DO J=1,NY1
          DO I=2,NX1-1
            FIELDI(:,I,J) = FIELDA(:,I-1,J)
          END DO
        END DO
      END IF
      FIELDI(:,1,J)=FIELDI(:,NX1-1,J)
      FIELDI(:,NX1,J)=FIELDI(:,2,J)
C****
      RETURN
      END SUBROUTINE AT2IT


      SUBROUTINE IT2AT(FIELDI,FIELDA,NF,QCONSERV)
!@sum  IT2AT interpolates Dynamic Ice Tracer grid to Atm Tracer grid
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : ima=>im,jma=>jm
      USE GEOM, only : imaxj
      USE ICEDYN, only : nx1,ny1,ricat
      IMPLICIT NONE
!@var QCONSERV true if integrated field must be conserved
      LOGICAL, INTENT(IN) :: QCONSERV
!@var N number of fields
      INTEGER, INTENT(IN) :: NF
!@var FIELDA array on atmospheric tracer grid
      REAL*8, INTENT(OUT), DIMENSION(NF,IMA,JMA) :: FIELDA
!@var FIELDI array on dynamic ice tracer grid
      REAL*8, INTENT(IN), DIMENSION(NF,NX1,NY1) :: FIELDI
      INTEGER I,J

C**** currently no need for interpolation,
C**** just scaling due to area differences for fluxes
C**** and shift one box in longitude
      IF (QCONSERV) THEN
        DO J=1,JMA
          DO I=1,IMAXJ(J)
            FIELDA(:,I,J) = FIELDI(:,I+1,J)*RICAT(J)
          END DO
        END DO
      ELSE
        DO J=1,JMA
          DO I=1,IMAXJ(J)
            FIELDA(:,I,J) = FIELDA(:,I+1,J)
          END DO
        END DO
      END IF
      DO I=2,IMA
        FIELDA(:,I,  1)=FIELDA(:,1,JMA)
        FIELDA(:,I,JMA)=FIELDI(:,1,  1)
      END DO
C****
      RETURN
      END SUBROUTINE IT2AT

      SUBROUTINE init_icedyn(iniOCEAN)
!@sum  init_STRAITS initializes strait variables
!@auth Gary Russell/Gavin Schmidt
!@ver  1.0
      USE ICEDYN, only : RSIX,RSIY,USI,VSI
      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: iniOCEAN

C**** Initiallise ice dynamic variables if ocean model starts
      if (iniOCEAN) THEN
        RSIX=0.
        RSIY=0.
        USI=0.
        VSI=0.
      end if

      RETURN
      END SUBROUTINE init_icedyn
