C**** QUSDEF.f

      MODULE QUSDEF
!@sum  QUSDEF contains definitions for manipulating moments
!@auth Maxwell Kelley
      integer, parameter ::
     &     nmom=9,
     &     mx =1, my =2, mz =3,
     &     mxx=4, myy=5, mzz=6,
     &     mxy=7, mzx=8, myz=9

C**** some useful vector subscripts:

C**** moments with a vertical component
      integer, dimension(4), parameter ::
     &     zmoms=(/MZ,MZZ,MYZ,MZX/)
C**** moments with no vertical component
      integer, dimension(5), parameter ::
     &     xymoms=(/MX,MY,MXX,MXY,MYY/)
C**** moments with a horizontal component
      integer, dimension(7), parameter ::
     &     ihmoms=(/MX,MY,MXX,MYY,MXY,MYZ,MZX/)
C**** moments with no horizontal component
      integer, dimension(2), parameter ::
     &     zomoms=(/MZ,MZZ/)

C**** x-x, x-y, x-z switches
      integer, dimension(nmom), parameter ::
     &     xdir=(/mx,my,mz,mxx,myy,mzz,mxy,myz,mzx/)
      integer, dimension(nmom), parameter ::
     &     ydir=(/my,mx,mz,myy,mxx,mzz,mxy,mzx,myz/)
      integer, dimension(nmom), parameter ::
     &     zdir=(/mz,my,mx,mzz,myy,mxx,myz,mxy,mzx/)

!@dbparam prather_limits forces +ve sub-grid scale profiles (default=0)
      integer :: prather_limits = 0

      Integer, Parameter :: FLUX_NEGATIVE=-1
      Integer, Parameter :: FLUX_NONNEGATIVE=+1

      END MODULE QUSDEF

      subroutine adv1d(s,smom, f,fmom, mass,dm, nx,qlimit,stride,dir
     *     ,ierr,nerr)
!@sum  adv1d implements the quadratic upstream scheme in one dimension
!@auth G. Russell, modified by Maxwell Kelley
c--------------------------------------------------------------
c adv1d advects tracers in x-direction using the qus
c the order of the moments in dir is: x,y,z,xx,yy,zz,xy,yz,zx
c--------------------------------------------------------------
      use QUSDEF, only : nmom,prather_limits
      implicit none
      !
!@var s      mean tracer amount (kg or J)
!@var smom   qus tracer moments (kg or J)
!@var f      tracer flux (diagnostic output) (kg or J)
!@var fmom   tracer moment flux (diagnostic output) (kg or J)
!@var mass   mass field (kg)
!@var dm     mass flux (kg)
!@var nx     length of 1D vector
!@var qlimit true if negative tracer is to be avoided
!@var stride spacing in s array between elements of relevant 1D array
!@var dir    direction switch (equals one of xdir ydir or zdir)
!@var ierr, nerr error codes
      integer, intent(in) :: nx,stride
      logical, intent(in) :: qlimit
      REAL*8, dimension(nx) :: dm, f
      REAL*8, dimension(nx*stride) :: s,mass
      REAL*8, dimension(nmom,nx*stride) :: smom
      REAL*8, dimension(nmom,nx) :: fmom
      integer, dimension(nmom) :: dir
      integer :: mx,my,mz,mxx,myy,mzz,mxy,myz,mzx
      integer :: n,np1,nm1,nn,ns
      integer,intent(out) :: ierr,nerr
      REAL*8 :: fracm,frac1,bymnew,mnew,dm2,tmp
      ! qlimit variables
      REAL*8 :: an, anm1, fn, fnm1, sn, sxn, sxxn
      !
      ierr=0 ; nerr=0
      mx  = dir(1)
      my  = dir(2)
      mz  = dir(3)
      mxx = dir(4)
      myy = dir(5)
      mzz = dir(6)
      mxy = dir(7)
      myz = dir(8)
      mzx = dir(9)
c-----------------------------------------------------------
      ! calculate tracer mass flux f
c-----------------------------------------------------------
      n = nx
      do np1=1,nx

         if(dm(n).lt.0.) then ! air mass flux is negative
            nn=np1
            frac1=+1.
         else                 ! air mass flux is positive
            nn=n
            frac1=-1.
         endif
         ns=1+(nn-1)*stride
         fracm=dm(n)/mass(ns)
         if(mass(ns).le.0.d0) fracm=0.d0
         frac1=fracm+frac1
         f(n)=fracm*(s(ns)-frac1*(smom(mx,ns)-
     &        (frac1+fracm)*smom(mxx,ns)))
      ! temporary storage of fracm in fx, to be used below
         fmom(mx,n)=fracm
      ! temporary storage of frac1 in fxx, to be used below
         fmom(mxx,n)=frac1
      !
        n = np1
      enddo
      if(qlimit) then
         nm1 = nx
         do n=1,nx
            ns=1+(n-1)*stride
            an = fmom(mx,n)      ! reading fracm which was stored in fx
            anm1 = fmom(mx,nm1)
            fn = f(n)
            fnm1 = f(nm1)
            sn = s(ns)
            sxn = smom(mx,ns)
            sxxn = smom(mxx,ns)
            call limitq(anm1,an,fnm1,fn,sn,sxn,sxxn,ierr)
            if (ierr.gt.0) then
              nerr=n
              if (ierr.eq.2) return
            end if
            f(n) = fn
            f(nm1) = fnm1
            smom(mx,ns) = sxn
            smom(mxx,ns) = sxxn
            nm1 = n
         enddo
      endif
c--------------------------------------------------------------------
      ! calculate tracer fluxes of slopes and curvatures
c--------------------------------------------------------------------
      n = nx
      do np1=1,nx
         if(dm(n).lt.0.) then ! air mass flux is negative
            nn=np1
         else                 ! air mass flux is positive
            nn=n
         endif
         ns=1+(nn-1)*stride
      ! retrieving fracm, which was stored in fx
         fracm=fmom(mx,n)
      ! retrieving frac1, which was stored in fxx
         frac1=fmom(mxx,n)
      !
         fmom(mx,n)=dm(n)*(fracm*fracm*(smom(mx,ns)
     &        -3.*frac1*smom(mxx,ns))-3.*f(n))
         fmom(mxx,n)=dm(n)*(dm(n)*fracm**3 *smom(mxx,ns)
     &        -5.*(dm(n)*f(n)+fmom(mx,n)))
      ! cross moments
         fmom(my,n)  = fracm*(smom(my,ns)-frac1*smom(mxy,ns))
         fmom(mxy,n) = dm(n)*(fracm*fracm*smom(mxy,ns)-3.*fmom(my,n))
         fmom(mz,n)  = fracm*(smom(mz,ns)-frac1*smom(mzx,ns))
         fmom(mzx,n) = dm(n)*(fracm*fracm*smom(mzx,ns)-3.*fmom(mz,n))
         fmom(myy,n) = fracm*smom(myy,ns)
         fmom(mzz,n) = fracm*smom(mzz,ns)
         fmom(myz,n) = fracm*smom(myz,ns)
         n = np1
      enddo
c-------------------------------------------------------------------
c update tracer mass, moments of tracer mass, air mass distribution
c-------------------------------------------------------------------
      nm1 = nx
      do n=1,nx
         ns=1+(n-1)*stride
         tmp=mass(ns)+dm(nm1)
         mnew=tmp-dm(n)
      !     mnew=mass(ns)+dm(nm1)-dm(n)
         bymnew = 1./mnew
         dm2=dm(nm1)+dm(n)
         tmp=s(ns)+f(nm1)
         s(ns)=tmp-f(n)
      !     s(ns)=s(ns)+f(nm1)-f(n)
      !
         smom(mx,ns)=(smom(mx,ns)*mass(ns)-3.*(-dm2*s(ns)
     &     +mass(ns)*(f(nm1)+f(n)))+(fmom(mx,nm1)-fmom(mx,n)))*bymnew
         smom(mxx,ns) = (smom(mxx,ns)*mass(ns)*mass(ns)
     &     +2.5*s(ns)*(mass(ns)*mass(ns)-mnew*mnew-3.*dm2*dm2)
     &     +5.*(mass(ns)*(mass(ns)*(f(nm1)-f(n))-fmom(mx,nm1)
     &     -fmom(mx,n))+dm2*smom(mx,ns)*mnew)
     &     +(fmom(mxx,nm1)-fmom(mxx,n))) * (bymnew*bymnew)
      ! cross moments
         smom(my,ns)=smom(my,ns)+fmom(my,nm1)-fmom(my,n)
         smom(mxy,ns)=(smom(mxy,ns)*mass(ns)-3.*(-dm2*smom(my,ns) +
     &        mass(ns)*(fmom(my,nm1)+fmom(my,n))) +
     &        (fmom(mxy,nm1)-fmom(mxy,n)))*bymnew
         smom(mz,ns)=smom(mz,ns)+fmom(mz,nm1)-fmom(mz,n)
         smom(mzx,ns)=(smom(mzx,ns)*mass(ns)-3.*(-dm2*smom(mz,ns) +
     &        mass(ns)*(fmom(mz,nm1)+fmom(mz,n))) +
     &        (fmom(mzx,nm1)-fmom(mzx,n)))*bymnew
      !
         smom(myy,ns)=smom(myy,ns)+fmom(myy,nm1)-fmom(myy,n)
         smom(mzz,ns)=smom(mzz,ns)+fmom(mzz,nm1)-fmom(mzz,n)
         smom(myz,ns)=smom(myz,ns)+fmom(myz,nm1)-fmom(myz,n)
      !
c------------------------------------------------------------------
         mass(ns) = mnew
         if(mass(ns).le.0.) then
            s(ns)=0.
            smom(:,ns)=0.
         endif
         if (qlimit .and. prather_limits.eq.1) then ! force Prather limits
           smom(mx,ns)=min(1.5*s(ns),max(-1.5*s(ns),smom(mx,ns)))
           smom(mxx,ns)=min(2.*s(ns)-abs(smom(mx,ns))/3.,max(abs(smom(mx
     *          ,ns))-s(ns),smom(mxx,ns)))
           smom(mxy,ns)=min(s(ns),max(-s(ns),smom(mxy,ns)))
           smom(mzx,ns)=min(s(ns),max(-s(ns),smom(mzx,ns)))
         end if
c-----------------------------------------------------------------
         nm1 = n
      enddo
      return
      end subroutine adv1d
c************************************************************************

      subroutine advection_1D_custom(s,smom, f,fmom, mass,dm, nx, 
     *     qlimit,stride,dir,ierr,nerr)
!@sum  advection_1d_custom is a parallel variant of adv1d which does not
!@sum  include qlimits.  This increases locality such that global 
!@sum  communication can be significantly reduced for advection along
!@sum  longitudes.
!@auth T. Clune
c--------------------------------------------------------------
c adv1d advects tracers in x-direction using the qus
c the order of the moments in dir is: x,y,z,xx,yy,zz,xy,yz,zx
c--------------------------------------------------------------
      use QUSDEF, only : nmom,prather_limits
      USE DOMAIN_DECOMP, only: grid, GET
      USE DOMAIN_DECOMP, only: NORTH, SOUTH
      USE DOMAIN_DECOMP, only: HALO_UPDATE, HALO_UPDATE_COLUMN
      USE DOMAIN_DECOMP, only: CHECKSUM
      implicit none
      !
!@var s      mean tracer amount (kg or J)
!@var smom   qus tracer moments (kg or J)
!@var f      tracer flux (diagnostic output) (kg or J)
!@var fmom   tracer moment flux (diagnostic output) (kg or J)
!@var mass   mass field (kg)
!@var dm     mass flux (kg)
!@var nx     length of 1D vector
!@var qlimit true if negative tracer is to be avoided
!@var stride spacing in s array between elements of relevant 1D array
!@var dir    direction switch (equals one of xdir ydir or zdir)
!@var ierr, nerr error codes
      integer, intent(in) :: nx,stride
      logical, intent(in) :: qlimit
      REAL*8, dimension(stride,0:nx+1) :: dm, f
      REAL*8, dimension(stride,0:nx+1) :: s,mass
      REAL*8, dimension(nmom,stride,0:nx+1) :: smom
      REAL*8, dimension(nmom,stride,0:nx+1) :: fmom
      integer, dimension(nmom) :: dir
      integer :: mx,my,mz,mxx,myy,mzz,mxy,myz,mzx
      integer :: n,np1,nm1,nn,ns
      integer,intent(out) :: ierr,nerr
      REAL*8 :: fracm,frac1,bymnew,mnew,dm2,tmp
      INTEGER :: i
      INTEGER :: J_0, J_1, J_0S, J_1S
      LOGICAL :: HAVE_SOUTH_POLE
      ! qlimit variables
      REAL*8 :: an, anm1, fn, fnm1, sn, sxn, sxxn
      !
      
      ierr=0 ; nerr=0
      mx  = dir(1)
      my  = dir(2)
      mz  = dir(3)
      mxx = dir(4)
      myy = dir(5)
      mzz = dir(6)
      mxy = dir(7)
      myz = dir(8)
      mzx = dir(9)

      CALL GET(grid, J_STRT = J_0, J_STOP=J_1, 
     &     J_STRT_SKP=J_0S, J_STOP_SKP=J_1S,
     &     HAVE_SOUTH_POLE=HAVE_SOUTH_POLE)
      J_1  = J_1 - J_0 + 1
      J_1S = J_1S- J_0 + 1
      J_0S = J_0S- J_0 + 1
      J_0 =  J_0 - J_0 + 1

      CALL HALO_UPDATE(grid, mass, FROM=NORTH)
      CALL HALO_UPDATE(grid, s,    FROM=NORTH)
      CALL HALO_UPDATE_COLUMN(grid, smom, FROM=NORTH)

      Do i=1, stride
        Call calc_tracer_mass_flux()
      End Do

      ! Limit fluxes to maintain positive mean values?
      CALL HALO_UPDATE_COLUMN(grid, fmom, FROM=NORTH+SOUTH)

      If (qlimit) Then
         
        Call HALO_UPDATE_COLUMN(grid, f, FROM=NORTH+SOUTH)
        Do i=1, stride
          Call apply_limiter()
        End Do
        Call HALO_UPDATE_COLUMN(grid, f, FROM=NORTH+SOUTH)
      End If
c--------------------------------------------------------------------
         ! calculate tracer fluxes of slopes and curvatures
c--------------------------------------------------------------------
      Do i=1, stride
         Call tracer_slopes_and_curvatures()
      End Do
c-------------------------------------------------------------------
c update tracer mass, moments of tracer mass, air mass distribution
c-------------------------------------------------------------------
      CALL HALO_UPDATE(grid,  f,    FROM=SOUTH)
      CALL HALO_UPDATE(grid, dm,    FROM=SOUTH)
      CALL HALO_UPDATE_COLUMN(grid, fmom, FROM=SOUTH)

      DO i = 1, stride
         Call update_tracer_mass()
      END DO

      return

      Contains

      Integer Function CheckFlux(flux) Result (sign_of_flux)
      USE QUSDEF, Only : FLUX_NEGATIVE, FLUX_NONNEGATIVE
      Real*8, Intent(In) :: flux
      If (flux < 0) Then
        sign_of_flux = FLUX_NEGATIVE
      Else
        sign_of_flux = FLUX_NONNEGATIVE
      End If
      End Function CheckFlux

      Integer Function NeighborByFlux(n, dm)
      USE QUSDEF, Only : FLUX_NEGATIVE
        Integer, Intent(In) :: n
        Real*8,  Intent(In) :: dm

        Integer :: nn

        If (CheckFlux(dm) == FLUX_NEGATIVE) Then ! air mass flux is negative
          nn=n+1
        else                                     ! air mass flux is positive
          nn=n
        endif
        NeighborByFlux = nn

      End Function NeighborByFlux

      Function FluxFraction(dm) Result(frac)
      USE QUSDEF, Only : FLUX_NEGATIVE
        Real*8, Intent(In) :: dm
        Real*8 :: frac

        If (CheckFlux(dm) == FLUX_NEGATIVE) Then
          frac = +1.
        Else ! Flux non negative
          frac = -1.
        End If

      End Function FluxFraction
      
      Function MassFraction(dm, mass) Result (fracm)
        Real*8, Intent(In) :: dm
        Real*8, Intent(In) :: mass
        Real*8 :: fracm

        If (mass > 0.0d0) Then
          fracm = dm / mass
        Else
          fracm = 0.d0
        End If
        
      End Function MassFraction
      
      Subroutine calc_tracer_mass_flux()

      Do n = J_0, J_1

        nn = NeighborByFlux(n, dm(i,n))
        fracm = MassFraction(dm(i,n), mass(i,nn))

        frac1 = fracm + FluxFraction(dm(i,n))

        f(i,n)=fracm*(s(i,nn)-frac1*(smom(mx,i,nn)-
     *                               (frac1+fracm)*smom(mxx,i,nn)))
      ! temporary storage of fracm in fx, to be used below
        fmom(mx,i,n)=fracm
      ! temporary storage of frac1 in fxx, to be used below
        fmom(mxx,i,n)=frac1
      !
      enddo
      End Subroutine calc_tracer_mass_flux

      Subroutine tracer_slopes_and_curvatures()

      Do n = J_0, J_1

        nn = NeighborByFlux(n, dm(i,n))
      ! retrieving fracm, which was stored in fx
        fracm=fmom(mx,i,n)
      ! retrieving frac1, which was stored in fxx
        frac1=fmom(mxx,i,n)
      !
        fmom(mx,i,n)=dm(i,n)*(fracm*fracm*(smom(mx,i,nn)
     &       -3.*frac1*smom(mxx,i,nn))-3.*f(i,n))
        fmom(mxx,i,n)=dm(i,n)*(dm(i,n)*fracm**3 *smom(mxx,i,nn)
     &       -5.*(dm(i,n)*f(i,n)+fmom(mx,i,n)))
      ! cross moments
         fmom(my,i,n)  = fracm*(smom(my,i,nn)-frac1*smom(mxy,i,nn))
         fmom(mxy,i,n) = dm(i,n)*(fracm*fracm*smom(mxy,i,nn)
     &                            -3.*fmom(my,i,n))
         fmom(mz,i,n)  = fracm*(smom(mz,i,nn)-frac1*smom(mzx,i,nn))
         fmom(mzx,i,n) = dm(i,n)*(fracm*fracm*smom(mzx,i,nn)
     &                            -3.*fmom(mz,i,n))
         fmom(myy,i,n) = fracm*smom(myy,i,nn)
         fmom(mzz,i,n) = fracm*smom(mzz,i,nn)
         fmom(myz,i,n) = fracm*smom(myz,i,nn)

      enddo
      End Subroutine tracer_slopes_and_curvatures

      Subroutine apply_limiter
      If (HAVE_SOUTH_POLE) Then
         n = J_0
         an = fmom(mx,i,n)      ! reading fracm which was stored in fx
         anm1 = 0
         fn = f(i,n)
         fnm1 = 0
         sn = s(i,n)
         sxn = smom(mx,i,n)
         sxxn = smom(mxx,i,n)
         call limitq(anm1,an,fnm1,fn,sn,sxn,sxxn,ierr)
         if (ierr.gt.0) then
            nerr=n
            if (ierr.eq.2) return
         end if
         f(i,n)   = fn
         f(i,n-1) = fnm1
         smom(mx,i,n) = sxn
         smom(mxx,i,n) = sxxn
      End If

      nm1=J_0S-1
      DO n = J_0S, J_1S+1

         an = fmom(mx,i,n)      ! reading fracm which was stored in fx
         anm1 = fmom(mx,i,nm1)
         fn = f(i,n)
         fnm1 = f(i,nm1)
         sn = s(i,n)
         sxn = smom(mx,i,n)
         sxxn = smom(mxx,i,n)
         call limitq(anm1,an,fnm1,fn,sn,sxn,sxxn,ierr)
         if (ierr.gt.0) then
            nerr=n
            if (ierr.eq.2) return
         end if
         f(i,n)   = fn
         f(i,nm1) = fnm1
         smom(mx,i,n) = sxn
         smom(mxx,i,n) = sxxn
         nm1 = n
      enddo

      End Subroutine apply_limiter

      Subroutine update_tracer_mass()

      If (HAVE_SOUTH_POLE) THEN
         n = J_0
!        mnew=mass(ns)+dm(nm1)-dm(n)
        tmp=mass(i,n)
        mnew=tmp-dm(i,n)

         bymnew = 1./mnew
         dm2=dm(i,n)

!       s(i,n)=s(i,n)+f(nm1)-f(n)
         tmp=s(i,n)
         s(i,n)=tmp-f(i,n)

         smom(mx,i,n)=(smom(mx,i,n)*mass(i,n)-
     &                  3.*(-dm2*s(i,n)
     &                      +mass(i,n)*(f(i,n)))+
     &                  (-fmom(mx,i,n)))*bymnew
         smom(mxx,i,n) = (smom(mxx,i,n)*mass(i,n)*mass(i,n)
     &     +2.5*s(i,n)*(mass(i,n)*mass(i,n)-mnew*mnew-3.*dm2*dm2)
     &     +5.*(mass(i,n)*(mass(i,n)*(-f(i,n))
     &     -fmom(mx,i,n))+dm2*smom(mx,i,n)*mnew)
     &     +(-fmom(mxx,i,n))) * (bymnew*bymnew)
      ! cross moments
         smom(my,i,n)=smom(my,i,n)-fmom(my,i,n)
         smom(mxy,i,n)=(smom(mxy,i,n)*mass(i,n)
     &                   -3.*(-dm2*smom(my,i,n) +
     &                       mass(i,n)*(fmom(my,i,n))) +
     &        (-fmom(mxy,i,n)))*bymnew
         smom(mz,i,n)=smom(mz,i,n)-fmom(mz,i,n)
         smom(mzx,i,n)=(smom(mzx,i,n)*mass(i,n)
     &                   -3.*(-dm2*smom(mz,i,n) +
     &        mass(i,n)*(+fmom(mz,i,n))) +
     &        (-fmom(mzx,i,n)))*bymnew
      !
         smom(myy,i,n)=smom(myy,i,n)-fmom(myy,i,n)
         smom(mzz,i,n)=smom(mzz,i,n)-fmom(mzz,i,n)
         smom(myz,i,n)=smom(myz,i,n)-fmom(myz,i,n)
      !
c------------------------------------------------------------------
         mass(i,n) = mnew
         if(mass(i,n).le.0.) then
            s(i,n)=0.
            smom(:,i,n)=0.
         endif
         if (qlimit .and. prather_limits.eq.1) then ! force Prather limits
           smom(mx,i,n)=min(1.5*s(i,n),max(-1.5*s(i,n),smom(mx,i,n)))
           smom(mxx,i,n)=min(2.*s(i,n)-abs(smom(mx,i,n))/3.,
     *          max(abs(smom(mx,i,n))-s(i,n),smom(mxx,i,n)))
           smom(mxy,i,n)=min(s(i,n),max(-s(i,n),smom(mxy,i,n)))
           smom(mzx,i,n)=min(s(i,n),max(-s(i,n),smom(mzx,i,n)))
         end if

c-----------------------------------------------------------------
         if (qlimit .and. prather_limits.eq.1) then ! force Prather limits
           smom(mx,i,n)=min(1.5*s(i,n),max(-1.5*s(i,n),smom(mx,i,n)))
           smom(mxx,i,n)=min(2.*s(i,n)-abs(smom(mx,i,n))/3.,
     *          max(abs(smom(mx,i,n))-s(i,n),smom(mxx,i,n)))
           smom(mxy,i,n)=min(s(i,n),max(-s(i,n),smom(mxy,i,n)))
           smom(mzx,i,n)=min(s(i,n),max(-s(i,n),smom(mzx,i,n)))
         end if

      End If

      Do n=J_0S,J_1
         nm1 = n-1
!        mnew=mass(ns)+dm(nm1)-dm(n)
        tmp=mass(i,n)+dm(i,nm1)
        mnew=tmp-dm(i,n)

         bymnew = 1./mnew
         dm2=dm(i,nm1)+dm(i,n)

!       s(i,n)=s(i,n)+f(nm1)-f(n)
         tmp=s(i,n)+f(i,nm1)
         s(i,n)=tmp-f(i,n)

         smom(mx,i,n)=(smom(mx,i,n)*mass(i,n)-
     &                  3.*(-dm2*s(i,n)
     &                      +mass(i,n)*(f(i,nm1)+f(i,n)))+
     &                  (fmom(mx,i,nm1)-fmom(mx,i,n)))*bymnew
         smom(mxx,i,n) = (smom(mxx,i,n)*mass(i,n)*mass(i,n)
     &     +2.5*s(i,n)*(mass(i,n)*mass(i,n)-mnew*mnew-3.*dm2*dm2)
     &     +5.*(mass(i,n)*(mass(i,n)*(f(i,nm1)-f(i,n))-fmom(mx,i,nm1)
     &     -fmom(mx,i,n))+dm2*smom(mx,i,n)*mnew)
     &     +(fmom(mxx,i,nm1)-fmom(mxx,i,n))) * (bymnew*bymnew)
      ! cross moments
         smom(my,i,n)=smom(my,i,n)+fmom(my,i,nm1)-fmom(my,i,n)
         smom(mxy,i,n)=(smom(mxy,i,n)*mass(i,n)
     &                   -3.*(-dm2*smom(my,i,n) +
     &                       mass(i,n)*(fmom(my,i,nm1)+fmom(my,i,n))) +
     &        (fmom(mxy,i,nm1)-fmom(mxy,i,n)))*bymnew
         smom(mz,i,n)=smom(mz,i,n)+fmom(mz,i,nm1)-fmom(mz,i,n)
         smom(mzx,i,n)=(smom(mzx,i,n)*mass(i,n)
     &                   -3.*(-dm2*smom(mz,i,n) +
     &        mass(i,n)*(fmom(mz,i,nm1)+fmom(mz,i,n))) +
     &        (fmom(mzx,i,nm1)-fmom(mzx,i,n)))*bymnew
      !
         smom(myy,i,n)=smom(myy,i,n)+fmom(myy,i,nm1)-fmom(myy,i,n)
         smom(mzz,i,n)=smom(mzz,i,n)+fmom(mzz,i,nm1)-fmom(mzz,i,n)
         smom(myz,i,n)=smom(myz,i,n)+fmom(myz,i,nm1)-fmom(myz,i,n)
      !
c------------------------------------------------------------------
         mass(i,n) = mnew
         if(mass(i,n).le.0.) then
            s(i,n)=0.
            smom(:,i,n)=0.
         endif
c-----------------------------------------------------------------
         if (qlimit .and. prather_limits.eq.1) then ! force Prather limits
           smom(mx,i,n)=min(1.5*s(i,n),max(-1.5*s(i,n),smom(mx,i,n)))
           smom(mxx,i,n)=min(2.*s(i,n)-abs(smom(mx,i,n))/3.,
     *          max(abs(smom(mx,i,n))-s(i,n),smom(mxx,i,n)))
           smom(mxy,i,n)=min(s(i,n),max(-s(i,n),smom(mxy,i,n)))
           smom(mzx,i,n)=min(s(i,n),max(-s(i,n),smom(mzx,i,n)))
         end if

      enddo
      End Subroutine Update_Tracer_Mass

      end subroutine advection_1D_custom

c************************************************************************

      subroutine limitq(anm1,an,fnm1,fn,sn,sx,sxx,ierr)
!@sum  limitq adjusts moments to maintain non-neg. tracer means/fluxes
!@auth G. Russell, modified by Maxwell Kelley
        implicit none
        REAL*8 :: anm1,an,fnm1,fn,sn,sx,sxx
c local variables
        REAL*8 :: sl,sc,sr, frl,frl1, frr,frr1, gamma,g13ab,
     &       fr,fr1, fsign,su,sd
        integer, intent(out) ::ierr
        ierr=0
c****
c**** modify the tracer moments so that the tracer mass in each
c**** division is non-negative
c****
c**** no air leaving the box
        if(anm1.ge.0. .and. an.le.0.) return
c**** air is leaving through both the left and right edges
        if(anm1.lt.0. .and. an.gt.0.) then
           sl = -fnm1
           sr = +fn
           sc = sn - sl
           sc = sc - sr
c**** all three divisions are non-negative
           if(sl.ge.0. .and. sr.ge.0. .and. sc.ge.0.) return
c**** first check for the cases when only one of the three is negative
           frl = anm1
           frl1 = frl+1.
           frr = an
           frr1 = frr-1.
           if(sl.ge.0. .and. sr.ge.0.) then ! center division
              gamma = 1.+(frl-frr)
              g13ab = gamma*gamma - 1. + 3.*(frl+frr)**2
              sxx = sxx - sc*10.*g13ab /
     &             (gamma*(12.*(frl+frr)**2 + 5.*g13ab*g13ab))
              sx = sx + sc*12.*(frl+frr) /
     &             (gamma*(12.*(frl+frr)**2 + 5.*g13ab*g13ab))
              sl = -frl*(sn-frl1*(sx-(frl+frl1)*sxx))
              sr = sn-sl
           else if(sr.ge.0.) then           ! leftmost division
              sxx = sxx + sl*(frl+frl1)/(frl*frl1*(.6d0+(frl+frl1)**2))
              sx = sn/frl1 + (frl+frl1)*sxx
              sr = frr*(sn-frr1*(sx-(frr+frr1)*sxx))
              sl = 0.
           else if(sl.ge.0.) then           ! rightmost division
              sxx = sxx - sr*(frr+frr1)/(frr*frr1*(.6d0+(frr+frr1)**2))
              sx = sn/frr1 + (frr+frr1)*sxx
              sl = -frl*(sn-frl1*(sx-(frl+frl1)*sxx))
              sr = 0.
           endif
           sc = sn - sl
           sc = sc - sr
c check for the cases where two of the three divisions are nonpositive
c these cases arise due to adjustments made when only one of the three
c divisions was negative
           if(sl.le.0. .and. sr.le.0.) then ! two outer divisions
              gamma = 1.+(frl-frr)
              sxx = sn*(1.+gamma)/(2.*gamma*frl1*frr1)
              sx  = sn*(frl+frr)*(1.+2.*gamma)/(2.*gamma*frl1*frr1)
              sl = 0.
              sr = 0.
           else if(sl.le.0. .and. sc.le.0.) then ! center/left divs
              sxx = sn/(2.*frr*frl1)
              sx  = sn*(frl+frr+.5)/(frr*frl1)
              sl = 0.
              sr = sn
           else if(sr.le.0. .and. sc.le.0.) then ! center/right divs
              sxx = sn/(2.*frl*frr1)
              sx  = sn*(frl+frr-.5)/(frl*frr1)
              sl = sn
              sr = 0.
           endif
           fnm1 = -sl
           fn   = +sr
        else
c**** air is leaving only through one edge
           if(an.gt.0.)  then ! right edge
              fr=an
              sd=fn
              fsign=-1.
           else                  ! left edge
              fr=anm1
              sd=-fnm1
              fsign=1.
           endif
           if(abs(fr).gt.1.)  then
c**** give warnings if fractional mass loss > 1
             ierr=1
             write(6,*) "limitq warning: abs(a)>1",fr,sd
             write(6,*) "limitq input: anm1,an,fnm1,fn,sn,sx,sxx",anm1
     *            ,an,fnm1,fn,sn,sx,sxx
c**** only stop if net tracer mass is negative
             if (sn+fnm1-fn.lt.0) then
               ierr=2
               write(6,*) "limitq error: new sn < 0",sn,sn+fnm1-fn
               return
             end if
           end if
           su = sn-sd
           if(sd.ge.0. .and. su.ge.0.) return
           fr1=fr+fsign
           if(sd.lt.0.)  then
c**** downstream division is negative
              sxx = sxx +fsign*sd*(fr+fr1)/(fr*fr1*(.6d0+(fr+fr1)**2))
              sx = sn/fr1 + (fr+fr1)*sxx
              su = sn
           else
c**** upstream division is negative
              sxx = sxx -fsign*su*(fr+fr1)/(fr*fr1*(.6d0+(fr+fr1)**2))
              sx = sn/fr + (fr+fr1)*sxx
              su = 0.
           endif
           sd = sn - su
           if(an.gt.0.) then
              fn=sd
           else
              fnm1=-sd
           endif
        endif
        return
      end subroutine limitq


      SUBROUTINE CTMIX (RM,RMOM,FMAIR,FMIX,FRAT)
!@sum CTMIX  Cloud top mixing of tracer moments (incl. q,t) from CONDSE
!@auth Gary Russell, Jean Lerner, Gavin Schmidt
!@ver 1.0
!@var RM,RX,RY,RZ,RXX,RYY,RZZ,RXY,RYZ,RZX mean and moments of tracer
      USE QUSDEF
      IMPLICIT NONE
      REAL*8, DIMENSION(2) :: RM
      REAL*8, DIMENSION(NMOM,2) :: RMOM
      REAL*8 FMIX      !@var FMIX  fraction of lower box mixed up
      REAL*8 FRAT      !@var FRAT  fraction of upper box mixed down
      REAL*8 FMAIR     !@var FMAIR mass of air mixed
      REAL*8 RTEMP     !@var RTEMP dummy variable

      RTEMP = RM(1)*(1.-FMIX)+FRAT*RM(2)
      RM(2) = RM(2)*(1.-FRAT)+FMIX*RM(1)
      RM(1) = RTEMP
C X
      RTEMP       = RMOM(MX ,1)*(1.-FMIX)+FRAT*RMOM(MX ,2)
      RMOM(MX ,2) = RMOM(MX ,2)*(1.-FRAT)+FMIX*RMOM(MX ,1)
      RMOM(MX ,1) = RTEMP
C Y
      RTEMP       = RMOM(MY ,1)*(1.-FMIX)+FRAT*RMOM(MY ,2)
      RMOM(MY ,2) = RMOM(MY ,2)*(1.-FRAT)+FMIX*RMOM(MY ,1)
      RMOM(MY ,1) = RTEMP
C XX
      RTEMP       = RMOM(MXX,1)*(1.-FMIX)+FRAT*RMOM(MXX,2)
      RMOM(MXX,2) = RMOM(MXX,2)*(1.-FRAT)+FMIX*RMOM(MXX,1)
      RMOM(MXX,1) = RTEMP
C YY
      RTEMP       = RMOM(MYY,1)*(1.-FMIX)+FRAT*RMOM(MYY,2)
      RMOM(MYY,2) = RMOM(MYY,2)*(1.-FRAT)+FMIX*RMOM(MYY,1)
      RMOM(MYY,1) = RTEMP
C XY
      RTEMP       = RMOM(MXY,1)*(1.-FMIX)+FRAT*RMOM(MXY,2)
      RMOM(MXY,2) = RMOM(MXY,2)*(1.-FRAT)+FMIX*RMOM(MXY,1)
      RMOM(MXY,1) = RTEMP
C Z MOMENTS
      RMOM(ZMOMS,:) = RMOM(ZMOMS,:)*(1.-FMAIR)

      END SUBROUTINE CTMIX
