!@sum  Routines to calculates solar zenith angles, weighted by time/sunlight
!@auth Original Development Team
!@auth M. Kelley
!@ver  1.0
!      module RAD_COSZ0
      USE CONSTANT, only : twopi,pi,teeny
      USE MODEL_COM, only : im,jm
      USE GEOM, only : lon2d
      USE RAD_COM, only : cosd,sind,sinj,cosj
      USE DOMAIN_DECOMP, ONLY: grid
      implicit none
      save
      private
      public :: coszs,coszt,cosz_init,update_dusk

c
c local variables common to these routines
c
      REAL*8 :: LT1,LT2,SLT1,SLT2,S2LT1,S2LT2,LT1_SV
C**** ZERO1 HAS TO EQUAL THE CUT-OFF VALUE FOR COSZ USED IN SOLAR
C**** COSZS WORKS CORRECTLY ONLY IF ZERO1 >> 1.D-3
      REAL*8, PARAMETER :: ZERO1=1.D-2
      INTEGER I,J
      REAL*8 ZERO2,DUSK,DAWN,CJCD,SJSD,DROT
      INTEGER :: I_0, I_1, J_0, J_1
      logical :: lwrap

      real*8, dimension(:,:), allocatable :: duskij,sinlatij,coslatij
      contains

      subroutine cosz_init
      implicit none
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP
      J_0 = grid%J_STRT
      J_1 = grid%J_STOP
      allocate(duskij(i_0:i_1,j_0:j_1))
      allocate(sinlatij(i_0:i_1,j_0:j_1))
      allocate(coslatij(i_0:i_1,j_0:j_1))
      do j=j_0,j_1
        sinlatij(:,j) = sinj(j)
        coslatij(:,j) = cosj(j)
      enddo
      return
      end subroutine cosz_init

      subroutine COSZT (ROT1,ROT2,COSZ)
C****
C**** THIS ENTRY COMPUTES THE ZENITH ANGLE WEIGHTED BY DAYTIME
C**** HOURS FROM ROT1 TO ROT2, GREENWICH MEAN TIME IN RADIANS.  ROT1
C**** MUST BE BETWEEN 0 AND 2*PI.  ROT2 MUST BE BETWEEN ROT1 AND
C**** ROT1+2*PI.  I=1 MUST LIE ON THE INTERNATIONAL DATE LINE.
C****

      IMPLICIT NONE
      REAL*8 ROT1,ROT2
      REAL*8, DIMENSION(grid%I_STRT_HALO:grid%I_STOP_HALO,
     &                  grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     &     COSZ

      DROT=ROT2-ROT1

C****
C**** LOOP OVER GRID BOXES
C****
      DO J=J_0,J_1
      DO I=I_0,I_1

      DUSK=DUSKIJ(I,J)
      DAWN=-DUSK

      LT1=ROT1+LON2D(I,J)
      LT2=ROT2+LON2D(I,J)

      call adjust_angles

      if(dawn.gt.lt2 .or. dusk.lt.lt1) then ! constant darkness in this cell
        cosz(i,j) = 0d0
      else
        SJSD=SINLATIJ(I,J)*SIND
        CJCD=COSLATIJ(I,J)*COSD
        lt1 = max(lt1,dawn)
        lt2 = min(lt2,dusk)
        slt1 = sin(lt1)
        slt2 = sin(lt2)
        COSZ(I,J)=SJSD*(LT2-LT1)+CJCD*(SLT2-SLT1)
        if(lwrap) then ! second contribution: lt1 to dusk
          lt1 = lt1_sv
          lt2 = dusk
          slt1 = sin(lt1)
          slt2 = sin(lt2)
          COSZ(I,J)=COSZ(I,J)+SJSD*(LT2-LT1)+CJCD*(SLT2-SLT1)
        endif
        COSZ(I,J)=COSZ(I,J)/DROT
      endif

      ENDDO ! I
      ENDDO ! J

      if(grid%have_south_pole) then
        cosz(:,1) = sum(cosz(:,1))/im
      endif
      if(grid%have_north_pole) then
        cosz(:,jm) = sum(cosz(:,jm))/im
      endif

      RETURN
      END subroutine COSZT

      subroutine COSZS (ROT1,ROT2,COSZ,COSZA)
C****
C**** THIS ENTRY COMPUTES THE ZENITH ANGLE TWICE, FIRST WEIGHTED BY THE
C**** DAYTIME HOURS FROM ROT1 TO ROT2 AND SECONDLY WEIGHTED BY THE
C**** INCIDENT SUN LIGHT FROM ROT1 TO ROT2.
C****
      IMPLICIT NONE
      REAL*8 ROT1,ROT2
      REAL*8, DIMENSION(grid%I_STRT_HALO:grid%I_STOP_HALO,
     &                  grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     &     COSZ,COSZA
      REAL*8 ECOSZ,ECOSQZ,ECOSZ1

      DROT=ROT2-ROT1

C****
C**** LOOP OVER GRID BOXES
C****
      DO J=J_0,J_1
      DO I=I_0,I_1

      DUSK=DUSKIJ(I,J)
      DAWN=-DUSK

      LT1=ROT1+LON2D(I,J)
      LT2=ROT2+LON2D(I,J)

      call adjust_angles

      if(dawn.gt.lt2.or.dusk.lt.lt1) then ! constant darkness in this cell
        cosz(i,j) = 0d0
        cosza(i,j) = 0d0
      else
        SJSD=SINLATIJ(I,J)*SIND
        CJCD=COSLATIJ(I,J)*COSD
        lt1 = max(lt1,dawn)
        lt2 = min(lt2,dusk)
        slt1 = sin(lt1)
        slt2 = sin(lt2)
        S2LT1 = 2.*SLT1*cos(lt1)
        S2LT2 = 2.*SLT2*cos(lt2)
        ECOSZ=SJSD*(LT2-LT1)+CJCD*(SLT2-SLT1)
        ECOSQZ=SJSD*ECOSZ+CJCD*(SJSD*(SLT2-SLT1)+
     *       .5*CJCD*(LT2-LT1+.5*(S2LT2-S2LT1)))
        if(lwrap) then ! second contribution: lt1 to dusk
          lt1 = lt1_sv
          lt2 = dusk
          slt1 = sin(lt1)
          slt2 = sin(lt2)
          S2LT1=2.*SLT1*cos(lt1)
          S2LT2=2.*SLT2*cos(lt2)
          ECOSZ1=SJSD*(LT2-LT1)+CJCD*(SLT2-SLT1)
          ECOSQZ = ECOSQZ + SJSD*ECOSZ1+CJCD*(SJSD*(SLT2-SLT1)+
     *         .5*CJCD*(LT2-LT1+.5*(S2LT2-S2LT1)))
          ECOSZ = ECOSZ + ECOSZ1
        endif
        COSZ(I,J)=ECOSZ/DROT
        COSZA(I,J)=ECOSQZ/ECOSZ
      endif

      ENDDO ! I
      ENDDO ! J

      if(grid%have_south_pole) then
        cosza(1,1) = sum(cosza(:,1)*cosz(:,1))
        cosz (1,1) = sum(cosz(:,1))
        cosza(:,1) = cosza(1,1)/(cosz(1,1)+teeny)
        cosz (:,1) = cosz(1,1)/im
      endif
      if(grid%have_north_pole) then
        cosza(1,jm) = sum(cosza(:,jm)*cosz(:,jm))
        cosz (1,jm) = sum(cosz(:,jm))
        cosza(:,jm) = cosza(1,jm)/(cosz(1,jm)+teeny)
        cosz (:,jm) = cosz(1,jm)/im
      endif

      return
      end subroutine COSZS

      subroutine update_dusk
      implicit none
C**** COMPUTE DUSK FOR EACH GRID CELL
      DO J=J_0,J_1
      DO I=I_0,I_1
        SJSD=SINLATIJ(I,J)*SIND
        CJCD=COSLATIJ(I,J)*COSD
        IF(SJSD-CJCD.GT.0.) THEN ! constant daylight
          DUSKIJ(I,J) = PI
        ELSEIF(SJSD+CJCD.LT.0.) THEN ! constant darkness
          DUSKIJ(I,J) = 0D0
        ELSE
          DUSKIJ(I,J)=ACOS(-SJSD/CJCD)
        ENDIF
      ENDDO
      ENDDO
      return
      end subroutine update_dusk

      subroutine adjust_angles
c
c domain is lt=-pi:+pi, noon is at lt=0
c
      implicit none
      do while(lt1.gt.pi)
        lt1 = lt1 - twopi
      enddo
      do while(lt2.gt.pi)
        lt2 = lt2 - twopi
      enddo
      lwrap = .false.
      if(lt1.gt.lt2) then
        if(lt1.lt.dusk .and. lt2.gt.dawn) then
c
c daylight at lt1 and lt2, brief darkness in between.
c computation will be divided into two parts:
c 1. dawn to lt2
c 2. lt1 to dusk
c
          lt1_sv = lt1
          lt1 = dawn
          lwrap = .true.
        elseif(lt1.lt.dusk .and. lt2.lt.dawn) then
          lt2 = dusk
        elseif(lt2.gt.dawn .and. lt1.gt.dusk) then
          lt1 = dawn
        endif
      endif
      end subroutine adjust_angles

      subroutine calc_asc(lam,phi,n,acell,scell,ccell)
c Calculate Area,Sine,Cosine
c
c calculates the area and area-averaged sin(lat),cos(lat) of
c an area whose n vertices lie at lon,lat=lam(1:n),phi(1:n)
c and whose edges are assumed to be great circles.
c
      implicit none
      integer :: n
      real*8, dimension(n) :: lam,phi
      real*8 :: acell,scell,ccell
      integer :: i,i1,i2
      real*8 :: gam,eps,x1,x2,asum,c2sum,scsum,phisum
      asum = 0d0
      c2sum = 0d0
      scsum = 0d0
      phisum = 0d0
      do i=1,n
        i1 = i
        if(i.lt.n) then
          i2 = i+1
        else
          i2 = 1
        endif
        if(abs(lam(i1)-lam(i2)).lt.1d-6) cycle
        call calc_gameps(lam(i1),lam(i2),phi(i1),phi(i2),
     &       gam,eps)
        x1 = -eps
        x2 = lam(i2)-lam(i1)-eps
        asum = asum + lineint_sinphi(x1,x2,gam)
        c2sum = c2sum + lineint_cosphi2(x1,x2,gam)
        scsum = scsum + lineint_sinphicosphi(x1,x2,gam)
        phisum = phisum + lineint_phi(x1,x2,gam)
      enddo
      acell = -asum
      scell = -.5*c2sum/asum
      ccell = .5*(phisum+scsum)/asum
      return
      end subroutine calc_asc


      subroutine calc_asc_pol(lam,phi,acell,scell,ccell)
c
c Variant of calc_asc for 4-sided cells having
c the north or south pole at one corner.
c

      implicit none
      real*8, dimension(4) :: lam,phi
      real*8 :: acell,scell,ccell
      real*8, dimension(3) :: lamshift,phishift
      real*8, dimension(5) :: lam5,phi5
      integer :: i,ii,ipole
      real*8 :: aeq,atot,ceq,seq
      real*8 :: pi

      pi = acos(-1d0)
c
c place the pole in position 4
c
      ipole = 0
      do i=1,4
        if(abs(phi(i)-.5*pi).lt.1d-6) then
          ipole = i
          exit
        elseif(abs(phi(i)+.5*pi).lt.1d-6) then
          ipole = i
          exit
        endif
      enddo
      if(ipole.eq.0) return
      ii=0
      do i=ipole+1,4
        ii = ii + 1
        lamshift(ii) = lam(i)
        phishift(ii) = phi(i)
      enddo
      do i=1,ipole-1
        ii = ii + 1
        lamshift(ii) = lam(i)
        phishift(ii) = phi(i)
      enddo

c
c calculate the integrals between the equator and the
c equatorial edges of the cell
c
      lam5(1) = lamshift(1)
      lam5(2) = lamshift(3)
      phi5(1:2) = 0d0
      do i=3,5
        lam5(i) = lamshift(6-i)
        phi5(i) = phishift(6-i)
      enddo
      call calc_asc(lam5,phi5,5,aeq,seq,ceq)

c
c obtain the properties of this cell as residuals
c
      atot = abs(lam5(2)-lam5(1))
      acell = atot-aeq
      ccell = (.25*pi*atot-aeq*ceq)/acell
      if(seq.lt.0.) then ! southern hemisphere
        scell = (-.5*atot-aeq*seq)/acell
      else
        scell = (+.5*atot-aeq*seq)/acell
      endif
      return
      end subroutine calc_asc_pol

      subroutine calc_gameps(lam1,lam2,phi1,phi2,gam,eps)
c given two latlon positions on a great circle, calculates
c the functional dependence of latitude upon longitude in
c the form
c tan(latitude) = gam*sin(x-eps)
c where x = longitude minus lam1
      implicit none
      real*8 :: lam1,lam2,phi1,phi2,gam,eps
      real*8 :: dlam,term1
      dlam = lam2-lam1
      term1 = tan(phi1)*cos(dlam)-tan(phi2)
      if(phi1.eq.0. .and. phi2.eq.0.) then
        eps = 0d0
      else
        eps = atan(tan(phi1)*sin(dlam)/term1)
      endif
      gam = term1/(cos(eps)*sin(-dlam))
      return
      end subroutine calc_gameps

      function lineint_sinphi(x1,x2,gam)
c
c integrates the function gam*sin(x)/sqrt(1+(gam*sin(x)**2))
c from x1 to x2
c
      implicit none
      real*8 :: lineint_sinphi
      real*8 :: x1,x2,gam
      real*8 :: a
      a = sqrt(1./(gam*gam)+1.)
      if(gam.ge.0.) then
        lineint_sinphi = asin(cos(x1)/a)-asin(cos(x2)/a)
      else
        lineint_sinphi = asin(cos(x2)/a)-asin(cos(x1)/a)
      endif
      return
      end function lineint_sinphi

      function lineint_cosphi2(x1,x2,gam)
c
c integrates the function 1/(1+(gam*sin(x)**2))
c from x1 to x2
c
      implicit none
      real*8 :: lineint_cosphi2
      real*8 :: x1,x2,gam
      real*8 :: b
      b = sqrt(1.+gam*gam)
      lineint_cosphi2 = (
     &     atan(b*tan(x2))-atan(b*tan(x1)))/b
      return
      end function lineint_cosphi2

      function lineint_sinphicosphi(x1,x2,gam)
c
c integrates the function gam*sin(x)/(1+(gam*sin(x)**2))
c from x1 to x2
c
      implicit none
      real*8 :: lineint_sinphicosphi
      real*8 :: x1,x2,gam
      real*8 :: y1,y2,b
      b = sqrt( (gam*gam)/(gam*gam+1.) )
      if(gam.ge.0.) then
        y1 = cos(x1)
        y2 = cos(x2)
      else
        y1 = cos(x2)
        y2 = cos(x1)
      endif
      lineint_sinphicosphi = .5*log(
     &     ((1d0-b*y2)*(1d0+b*y1))/((1d0+b*y2)*(1d0-b*y1)))/
     &     sqrt(1.+gam*gam)
      return
      end function lineint_sinphicosphi

      function lineint_phi(x1,x2,gam)
c
c numerically integrates the function atan(gam*sin(x))
c from x1 to x2
c
      implicit none
      real*8 :: lineint_phi
      real*8 :: x1,x2,gam
      real*8, parameter :: dxmax=.0025d0
      real*8 :: tmpsum,dx,x,b,bcosx,sinx,cosx,cosx_sv,sindx,cosdx
      integer :: i,nx
      if(gam.eq.0.) then
        lineint_phi = 0.
      else
        b = sqrt(gam*gam/(1d0+gam*gam))
        if(gam.lt.0.) b=-b
        nx = 10
        dx = min(dxmax,abs(x2-x1)/nx)
        nx = max(nx,int(abs(x2-x1)/dx))
        dx = (x2-x1)/nx
        x = x1-.5*dx
        sindx = sin(dx)
        cosdx = cos(dx)
        sinx = sin(x)
        cosx = cos(x)
        tmpsum = 0d0
        do i=1,nx
          x = x+dx
          cosx_sv = cosx
          cosx = cosdx*cosx - sindx*sinx
          sinx = cosdx*sinx + sindx*cosx_sv
          bcosx = b*cosx
          tmpsum = tmpsum + x*bcosx/(1d0-bcosx**2)
        enddo
        lineint_phi = -dx*tmpsum/sqrt(1d0+gam*gam) +
     &       x2*atan(gam*sin(x2))-x1*atan(gam*sin(x1))
      endif
      end function lineint_phi

!      end module RAD_COSZ0
