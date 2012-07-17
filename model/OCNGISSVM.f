C****  
C**** GISSVM.f  GISS ocean vertical mixing scheme    2012/03/21
C****
#include "rundeck_opts.h"

      MODULE GISSMIX_COM
!@sum  GISSMIX_COM holds variables related to the GISS mixing scheme
!@auth Ye Cheng
#ifdef TRACERS_OCEAN
c     USE OCN_TRACER_COM, only : ntm
#endif
      USE OCEAN, only : im,jm,lmo
      USE SW2OCEAN, only : lsrpd
      IMPLICIT NONE
      SAVE
!@var otke turbulent kinetic energy in ocean (m/s)**2
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: otke
!@var otke_init_max maximum initial vaule of otke
      real*8, parameter :: otke_init_max=0.5d0/800. 
      real*8, parameter :: emin=1d-6,emax=1000. ! (m/s)^2
#ifdef TRACERS_OCEAN
c     REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: TRMO1,TXMO1,TYMO1
#endif
      ! bottom drag over lon,lat used by bottom drag routine
      !@ref C2010 section 7.2 equation (72)
!@var rhobot in-situ density at ocean bottom (kg/m^3)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: rhobot
!@var taubx x component of tau_b
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: taubx
!@var tauby y component of tau_b
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: tauby
!@var exya internal tidal energy (w/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: exya
!@var ut2a unresolved bottom velocity squared (m/s)^2
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: ut2a
      integer, parameter :: idrag=1    !@var idrag 1: tides are not explicit
!@var nonlocal 0: local; 1: non-local; 2: non-local with counter-gradient terms
      integer, parameter :: nonlocal=0
      END MODULE GISSMIX_COM


      MODULE GISS_OTURB
!@sum GISS_OTURB contains variables and routines for GISS mixing scheme
!@ref Canuto et al. 2010, Ocean Modelling, 34, 70-91 (C2010)
!@+   Canuto et al. 2011, Ocean Modelling, 36, 198-207 (C2011)
!@+   Cheng et al. 2002, JAS, 59,1550-1565 (C2002)
!@auth AHOWARD and YCHENG
!@date May 2011

      USE OCEAN, only : im,jm,lmo
      USE DOMAIN_DECOMP_1d, Only: AM_I_ROOT
      USE CONSTANT, only : omega,by3,grav
      USE GISSMIX_COM, only : otke,otke_init_max,emin,emax
      USE GISSMIX_COM, only : rhobot,taubx,tauby,exya,ut2a,idrag
     &    ,nonlocal

      implicit none

      integer, parameter :: mt=214      !@var mt dim of ri in table
c     integer, parameter :: nt=162      !@var mt dim of rr in table
      integer, parameter :: nt=108      !@var mt dim of rr in table
      real*8, parameter :: kmin=1d-3    !@var kmin min of diffusivities, (m^2/s)
      real*8, parameter :: kmax=100.   !@var kmax max of diffusivities, (m^2/s)
      real*8, parameter :: lrmax=1.    !@var lrmax max of length/l2
      real*8, parameter :: osocb1=21.6,kappa=0.4
c     real*8, parameter :: gammam_bg=.46
      real*8, save ::
     &    ria(mt)       !@var ria ri 1d array of richardson #, for 2d tables
     &   ,rra(nt)       !@var rra rr 1d array of density ratio,for 2d tables
     &   ,gma(mt,nt)    !@var gma 2d table for gm=(tau*shear)^2
     &   ,sma(mt,nt)    !@var sma 2d table for sm=structure fuction for momentum
     &   ,sha(mt,nt)    !@var sha 2d table for sh=structure fuction for heat
     &   ,ssa(mt,nt)    !@var ssa 2d table for ss=structure fuction for salinity
     &   ,rfa(mt,nt)    !@var rfa 2d table for rf=flux richardson #
     &   ,phim2a(mt,nt) !@var phim2a 2d table for phim2, used for bottom shear
     &   ,lr1a(mt)      !@var lr1a 1d table for length/(blackadar length)
c    &   ,lra(mt,nt)    !@var lra 2d table for length/(blackadar length)
c     real*8, save ::
c    &    gma1(mt)    !@var gma 1d table for gm=(tau*shear)^2
c    &   ,sma1(mt)    !@var sma 1d table for sm=structure fuction for momentum
c    &   ,sha1(mt)    !@var sha 1d table for sh=structure fuction for heat
c    &   ,ssa1(mt)    !@var ssa 1d table for ss=structure fuction for salinity
c    &   ,rfa1(mt)    !@var rfa 1d table for rf=flux richardson #
c     real*8, save ::
c    &    ribga(nt)     !@var ribg 1d array of backgound richardson number ri 
c    &   ,smbga(nt)     !@var smbga 1d table for background sm
c    &   ,shbga(nt)     !@var shbga 1d table for background sh
c    &   ,ssbga(nt)     !@var ssbga 1d table for background ss
c    &   ,rfbga(nt)     !@var rfbga 1d table for background rf
c    &   ,gammbga(nt)   !@var gammbga 1d table for background gamma_m
      real*8, save :: 
     &    rimin         !@var rimin min of richardson # ri
     &   ,rimax         !@var rimax max of richardson # ri
     &   ,rrmin         !@var rrmin min of density ratio rr
     &   ,rrmax         !@var rrmax max of density ratio rr

c     real*8, save :: rtwi_rr     ! Rrho, passed to fct in rtwi routine
      real*8, save ::  pi10,pi20,pi30,pi40,pi50
      real*8, save ::  p1,p2,p3,p4
      real*8, save ::  rf1,rf2

      CONTAINS

      subroutine gissmix_init(iniOCEAN)
!@sum creates tables for the giss diffusivities(ri,rr) 
!@Author AHoward/YCheng
!@ref Canuto et al., Ocean Modelling, 2010
!@+   Canuto et al., Ocean Modelling, 2010
!@+   Cheng et al., JAS, 2002
!@date May 2011    
      USE FILEMANAGER
      USE DOMAIN_DECOMP_1D, only : READT_PARALLEL
      USE OCEANR_DIM, only : ogrid

      implicit none

      ! in:
      logical, intent(in) :: iniOCEAN
      ! local:

      integer m,n,m0,m1,n0,j,k,i,l
      real*8 tmpa(mt)
      real*8 rini,rend,ff,tmp,pow
      real*8 gm1,sm1,sh1,ss1,sr

      real*8 ri,rr,gmbg
      integer iend,ier
      real*8 :: val,rest,eps
      character*80 path 
      INTEGER :: iu_TIDES
c     integer i_0,i_1,j_0,j_1

      ! establish ri grids for look-up table
      ! the grids extend from -rend to +rend, with more grids near zero
      rini=1d-4
      rend=1d4
      ff=4.
      m0=nint(log(rend/rini)/dlog(2.d0)*ff+1)
      m=2*m0
      if(m.ne.mt) then
         write(*,*) "m=",m, "    mt=",mt
         write(*,*) "Stop, m must equal mt"
         stop
      endif
      do j=1,m0-1
         tmpa(j)=rini*2**(dfloat(j-1)/ff)
      end do
      tmpa(m0)=rend
      do j=1,m
         if(j.le.m0) then
            ria(j)=-tmpa(m0-j+1)
         else
            ria(j)=tmpa(j-m0)
         endif
      end do

      ! establish rr grids for look-up table
      ! the grids extend from -rend to +rend, with more grids near 0
      ! abs(rr)>=rini
       rini=1d-2
       rend=1.-1d-2
       ff=8.
       n0=nint(log(rend/rini)/dlog(2.d0)*ff+1)
       n=2*n0
       if(n.ne.nt) then
          write(*,*) "n=",n, "    nt=",nt
          write(*,*) "Stop,  n must equal nt"
          stop
       endif
       do k=1,n0-1
          tmpa(k)=rini*2**(dfloat(k-1)/ff)
       end do
       tmpa(n0)=rend
       do k=1,n
          if(k.le.n0) then
             rra(k)=-tmpa(n0-k+1)
          else
             rra(k)=tmpa(k-n0)
          endif
       end do

      ! 2d tables for gm,sm,sh,ss,rf,phim2
      ! phim2 is defined in C2002, (36), with l=kz; for bottom shear

      rimin=ria(1) ! -1.d4
      rimax=ria(m) !  1.d4
      rrmin=rra(1) ! -0.99
      rrmax=rra(n) !  0.99
c     if( AM_I_ROOT() ) then
c        write(67,'(a,1e14.4)') "rimin=",rimin
c        write(67,'(a,1e14.4)') "rimax=",rimax
c        write(67,'(a,1e14.4)') "rrmin=",rrmin
c        write(67,'(a,1e14.4)') "rrmax=",rrmax
c     endif
      do k=1,n
         do j=1,m
            call gissdiffus(
            ! In:
     &         ria(j),rra(k)
            ! Out:
     &        ,gma(j,k),sma(j,k),sha(j,k),ssa(j,k),rfa(j,k))
            phim2a(j,k)=2./osocb1**2*gma(j,k)**.5d0/sma(j,k)
c           if( AM_I_ROOT() ) then
c              write(61,'(9e14.4)') ria(j),rra(k)
c    &        ,gma(j,k),sma(j,k),sha(j,k),ssa(j,k),rfa(j,k)
c           endif
         end do ! loop j
      end do ! loop k

      ! 1d tables for gm,sm,sh,ss,rf (rf used in length scale formula)

      do j=m,1,-1
         call gissdiffus(
         ! In:
     &      ria(j),rrmin
         ! Out:
     &     ,gm1,sm1,sh1,ss1,rf1)
         if(j.eq.m) rf2=rf1
         tmp=min(max(1-rf1/rf2,1d-30),1.d0)
         lr1a(j)=tmp**(4.*by3)
c        if( AM_I_ROOT() ) then
c           write(61,'(9e14.4)') ria(j),lr1a(j),rf2
c        endif
      end do ! loop j

c     ! 2d table for length/(blackadar leng)
c     call gissdiffus(
c     ! In:
c    &   rimax,rrmin
c     ! Out:
c    &  ,gm1,sm1,sh1,ss1,rf1)
c        if( AM_I_ROOT() ) then
c           write(61,'(9e14.4)') rf1,rf2,rf1/rf2
c        endif
c     do j=1,m
c        if( AM_I_ROOT() ) then
c           write(61,'(9e14.4)') ria(j),rf1,rfa1(j)/rf1
c        endif
c     end do

c     do k=1,n
c        do j=1,m
c           if(rfa(j,k)/rf1.gt.1.) then
c             write(*,*) "rfa(j,k)/rf1.gt.1., stop"
c             write(*,*) rfa(j,k)/rf1
c             stop
c           endif
c           if(rfa(j,k).gt.0.) then
c              pow=4./3.
c           else
c              pow=.15
c           endif
c           tmp=1.-4./(3.*pow)*rfa(j,k)/rf1
c           lra(j,k)=tmp**pow
c           lra(j,k)=min(tmp**pow,lrmax)
c        end do ! loop j
c     end do ! loop k
c     ! make the 1-d tables for background sm,sh,ss
c     eps=1.d-6
c     iend=40
c     ! rtwi finds the root of ri=fct(ri)
c     do k=1,n
c        rr=rra(k)
c        ! to estimate rest
c        if(k.eq.1) then
c            rest=.4
c        else
c            rest=ri
c        endif
c        rtwi_rr=rr
c        ! in rtwi, rest is the input (estimate of ri), ri is the output
c        ! fct is the function ri=fct(ri) from which ri is solved
c        call rtwi(ri,val,fct,rest,eps,iend,ier)
c        ribga(k)=ri
c        ! make 1d table for sm,sh,ss of rr
c        ribga(k)=0.4
c        call gissdiffus(
c        ! In:
c    &      ribga(k),rra(k)
c        ! Out:
c    &     ,gmbg,smbga(k),shbga(k),ssbga(k),rfbga(k))
c        gammbga(k)=.5*ribga(k)*gmbg*smbga(k)
c        if( AM_I_ROOT() ) then
c           write(62,'(9e14.4)') ribga(k),rra(k)
c    &      ,gammbga(k),smbga(k),shbga(k)/smbga(k),ssbga(k)/smbga(k)
c        endif
c     end do

C**** initialize exya, ut2a
      ! tidally induced diffusivities: C2010, (69), (75)-(77)
      call openunit("TIDES",iu_TIDES,.true.,.true.)
      CALL READT_PARALLEL(ogrid,iu_TIDES,NAMEUNIT(iu_TIDES),exya ,1)
      CALL READT_PARALLEL(ogrid,iu_TIDES,NAMEUNIT(iu_TIDES),ut2a ,1)
c     ! i_0=ogrid%I_STRT
c     ! i_1=ogrid%I_STOP
c     j_0=ogrid%J_STRT
c     j_1=ogrid%J_STOP
c     do i=1,im
c        do j=j_0,j_1
c        end do
c     end do

C**** initialize rhobot, taubx, tauby
      rhobot(:,:)=0.
      taubx(:,:)=0.
      tauby(:,:)=0.

C**** initialize otke
      if(iniOCEAN) then
         do l=1,lmo
            otke(l,:,:)=min(max(otke_init_max/(float(l)**2),emin),emax)
         end do
      endif

      return
      end subroutine gissmix_init

c     real*8 function fct(ri)
c!@sum rtwi finds the root of ri=fct(ri)
c      implicit none
c      real*8 ri,rr,gm,sm,sh,ss,rf
c      rr=rtwi_rr
c      call gissdiffus(
c      ! In:
c     &    ri,rr
c      ! Out:
c     &   ,gm,sm,sh,ss,rf)
c      fct=2.*gammam_bg/(gm*sm)
c      return
c      end function fct


      subroutine gissdiffus(
      ! In:
     &    ri,rr
      ! Out:
     &   ,gm,sm,sh,ss,rf)
!@sum finds structure functions sm,sh,ss and rf using giss model
!@ref Canuto et al., Ocean Modelling, 2010
!@auth aHoward/YCheng
!@date May, 2011
      
      implicit none

      ! in:
      real*8 :: ri,rr
      ! out:
      real*8 :: gm,sm,sh,ss,rf
      ! local:
      real*8 :: pi1,pi2,pi3,pi4,pi5,a3,a2,a1,a0
      real*8 :: tmp,a,b,c
      real*8 :: x1r,x2r,x3r,x1i,x2i,x3i,omrr
      real*8 :: x,q,p,bygam,xx,am,ah,as,w2byk,sr
      complex*16 x1,x2,x3

      ! prepare for solving follwoing cubic eqn for gm
      ! a3*gm**3+a2*gm**2+a1*gm+a0=0

      call cubic_coeff(
      ! In:
     &    ri,rr
      ! Out:
     &   ,pi1,pi2,pi3,pi4,pi5
     &   ,a3,a2,a1,a0)
      if(abs(a3).lt.1d-10.and.abs(a2).lt.1d-10) then
         gm=-a0/a1
      elseif(abs(a3).lt.1d-10) then
         ! solve quadratic a2*x**2+a1*x+a0=0
         ! x=gm
         tmp=sqrt(a1**2-4*a2*a0)
         gm=(-a1-tmp)/(2*a2)
         if(gm.lt.0.) then
            gm=(-a1+tmp)/(2*a2)
         endif
      else
         ! solve cubic x**3+a*x**2+b*x+c=0
         ! x=gm
         a=a2/a3
         b=a1/a3
         c=a0/a3
         call cubic_r(a,b,c,x1,x2,x3)
         x1r=dreal(x1); x1i=imag(x1);
         x2r=dreal(x2); x2i=imag(x2);
         x3r=dreal(x3); x3i=imag(x3);
         if(x1r.lt.0.) x1r=1d40
         if(x2r.lt.0.) x2r=2d40
         if(x3r.lt.0.) x3r=3d40
         if(abs(x1i).gt.1d-30) x1r=1d30
         if(abs(x2i).gt.1d-30) x2r=2d30
         if(abs(x3i).gt.1d-30) x3r=3d30
         gm=min(x1r,min(x2r,x3r))
         if(gm.ge.1d30) then
            write(*,*) "gm=",gm
            write(*,*) "why gm >= 1d30; stop and check"
            stop
         endif
      endif
      if(gm.lt.0.) then
         write(*,*) "gm < 0", gm
         stop
      endif
      ! below, x=gh, notation change (it was used to represent gm)
      omrr=1-rr
      x=ri*gm/omrr
      q=pi1*(pi2*(1+rr)-pi3*rr)
      p=pi4*(pi5-pi2*(1+rr))
      bygam=rr/( (pi4/pi1)*(1+q*x)/(1+p*x) ) ! lower case gamma
      ah=pi4/(1+p*x+pi2*pi4*x*(1-bygam))
      xx=(1-bygam)*x*ah
      am=2/gm*(15./7.+xx)
      w2byk=2. / (30./7.+xx)
      as=ah/( (pi4/pi1)*(1+q*x)/(1+p*x) )
      sh=ah*w2byk
      ss=as*w2byk
      sm=am*w2byk
      if((gm.le.0.).or.(sm.le.0.).or.(sh.le.0.).or.(ss.le.0.)) then
         write(100,*) "gm,sm,sh,ss <= 0., stop and check"
         write(100,*)  gm,sm,sh,ss
         stop
      endif
      sr=(sh-rr*ss)/(1-rr)
      rf=ri*sr/sm

      return
      end subroutine gissdiffus

      subroutine cubic_coeff(
         ! In:
     &      ri,rr
         ! Out:
     &     ,pi1,pi2,pi3,pi4,pi5
     &     ,a3,a2,a1,a0
     &     )

      implicit none

      real*8 ri,rr
      real*8 pi1,pi2,pi3,pi4,pi5
      real*8 a3,a2,a1,a0

      real*8 aa1,aa2,aa3,aa4,aa5,aa6
      real*8 sgmt0
      real*8 a,b,b2,rrsy,rrs,rrm,rrm1,rrm2,rrm1s,rr2,rrss
      real*8 Rrs2,bs2,Rrns,tmp
      real*8 pi2_,pi1_,pi4_,b2su,rr2su,pi1_1,pi4_1,pi2_1,w
      integer test

      !@ pi0's:

      sgmt0=.72d0
      pi20=1./3.d0
      pi40=1./5.*(1./(1.+1/sgmt0))
      pi10=pi40
      pi50=sgmt0
      pi30=pi50
      pi3=pi30
      pi5=pi50
      if(ri.lt.0.) then
         pi1=pi10
         pi2=pi20
         pi4=pi40
      else ! ri >= 0
         if(rr.gt.0.) then
c           a=5.
            b=2.*ri/(.1+ri)
            rrsy = 2./(rr + 1/rr)
            ! pi1 = pi10/(1. + ri/(1. + (a*rrsy**2)))
            tmp=atan2(.5*rrsy**2*(1-rrsy**2),2*rrsy**2-1)/3.1415927
            pi1 = pi10/(1.+ri**tmp*tmp)   !ok
            pi4 = pi1
            pi2=pi20*( 1-b*rrsy*(1-rrsy) )
         else ! rr < 0
            pi1=pi10/(1+ri)
            pi2=pi20
            pi4=pi40/(1+ri)
         endif
      endif

!@    a3*gm**3+a2*gm**2+a1*gm+a0=0

      aa1 = pi1 * pi4 *
     &  (    pi2 * (15 * pi3 + 7) *(rr**2 +1) + (14 * pi2
     &  - 14 * pi3 - 15 * pi3 ** 2) * rr    )
     & * (pi1 *  rr - pi4) / (-1 + rr) ** 3 / 150

      aa2 = pi1 * pi4 * 
     &  (    pi2 * (-150 * pi3 + 210 * pi1 + 7) *( rr **2+1)
     &    +( 14*(pi2-pi3)*(1+15*(pi1+pi4))+150*pi3**2 )*rr
     &    + 210* pi2 * ( pi4-pi1)    )
     &  / (-1 + rr) ** 2  / 9000

      aa3 = (pi1 * ( 5*pi2*pi4*(30 * pi3  + 17 )
     &      + pi1*(15 * pi3 + 7) ) * (rr ** 2+1)
     &   - ( 10*pi1*pi3*pi4*(15*pi3+17) + 15*pi2*(pi4**2+pi1**2)
     &    + 14*pi1*pi4*(1-10*pi2) ) * rr
     &       - (15*pi3+7)*(pi1**2-pi4**2))
     &      / (-1 + rr) ** 2 / 150

      aa4 = (( 150 *(pi1*pi3+pi4*pi2)-7 *pi1*(1+30*pi1) )*rr
     &    - 150*(pi1*pi2+pi4*pi3)+7 *pi4*(1 + 30 * pi4)  )
     &      / (1 - rr) / 9000

      aa5 =  ((-17 * pi1 - 30 *( pi1 * pi3 +  pi4 * pi2)) * rr  
     &          +17 * pi4 + 30 *( pi1 * pi2 +  pi4 * pi3))
     &    / (1 - rr) / 30

      aa6 = -1. / 60.

      a3 = aa1*ri**3+aa2*ri**2
      a2 = aa3*ri**2+aa4*ri
      a1 = aa5*ri+aa6
      a0 = 1.

      ! for the turbulence soc model at level 2.5:
      p1=pi40*pi50
      p2=pi40*(pi50-pi20)
      p3=pi10*pi20
      p4=pi40-.02d0*by3

      return
 1001 format(15(1pe14.4))
      end subroutine cubic_coeff
      
      subroutine cubic_r(a,b,c,x1,x2,x3)
!@sum solves the cubic eqn x**3+a*x**2+b*x+c=0
!@+   for a, b, c are real
      implicit none
      real*8 a,b,c
      complex*16 x1,x2,x3
      real*8 pi,q,r,the,temp,aa,bb
      pi=acos(-1.d0)
      q=(a**2-3.d0*b)/9.d0
      r=(2.d0*a**3-9.d0*a*b+27.d0*c)/54.d0
      if (r**2.lt.q**3) then
          the=acos(r/sqrt(q**3))
          temp=2.d0*sqrt(q)
          x1=-temp*cos(the/3.d0)-a/3.d0
          x2=-temp*cos((the+2.d0*pi)/3.d0)-a/3.d0
          x3=-temp*cos((the-2.d0*pi)/3.d0)-a/3.d0
      else
          aa=-sign(1.d0,r)*(abs(r)+sqrt(r**2-q**3))**(1.d0/3.d0)
          if(aa.ne.0.d0) then
              bb=q/aa
          else
              bb=0.d0
          endif
          x1=(aa+bb)-a/3.d0
          x2=cmplx(-0.5d0*(aa+bb)-a/3.d0, sqrt(3.d0)/2.d0*(aa-bb))
          x3=conjg(x2)
      endif
      return
      end subroutine cubic_r


      SUBROUTINE locate(n,xa,x,klo,khi,a,b)
!@sum locate finds the grids that embrace x
!@+   after call locate, the interpreted value can be calculated as
!@+   y=a*ya(klo)+b*ya(khi)
      implicit none
      ! in:
      INTEGER n
      REAL*8 xa(n),x
      ! out:
      INTEGER klo
      INTEGER khi
      REAL*8 a ! a=(xa(khi)-x)/h
      REAL*8 b ! b=(x-xa(klo))/h
      ! local:
      INTEGER k
      REAL*8 h
      klo=1
      khi=n
      do while (khi-klo.gt.1)
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      end do
      h=xa(khi)-xa(klo)
      if (h.eq.0.) then
         write(*,*) 'bad xa input in locate'
         stop
      endif
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
c     y=a*ya(klo)+b*ya(khi) ! calc. outside locate for efficiency
      return
      END SUBROUTINE locate


c      subroutine rtwi(x,val,fct,xst,eps,iend,ier)
cc     to solve general nonlinear equations of the form x=fct(x)
cc     by means of wegsteins iteration method
cc     prepare iteration
c
c      implicit none
c
c      real*8 x,val,fct,xst,eps
c      integer iend,ier
c      real*8 tol,a,b,d
c      integer i
c
c      ier=0
c      tol=xst
c      x=fct(tol)
c      a=x-xst
c      b=-a
c      tol=x
c      val=x-fct(tol)
cc     start iteration loop
c      do 6 i=1,iend
c      if(val) 1,7,1
cc     equation is not satisfied by x
c 1    b=b/val-1.
c      if(b) 2,8,2
cc     iteration is possible
c 2    a=a/b
c      x=x+a
c      b=val
c      tol=x
c      val=x-fct(tol)
cc     test on satisfactory accuracy
c      tol=eps
c      d=abs(x)
c      if(d-1.) 4,4,3
c 3    tol=tol*d
c 4    if(abs(a)-tol) 5,5,6
c 5    if(abs(val)-10.*tol) 7,7,6
c 6    continue
cc     end of iteration loop
cc     no convergence after iend iteration steps. error return.
c      ier=1
c 7    return
cc     error return in case of zero divisor
c 8    ier=2
c      return
c      end subroutine rtwi

      END MODULE GISS_OTURB


      subroutine gissmix( 
      ! in:
     &    n,ze,zg,db,dv2,adt,bds,to,s,rho,uob,vob,u2b
     &   ,fc,ustar,bf,hbl,kbl,strait,ilon,jlat
      ! inout:
     &   ,e
      ! out:
     &   ,ri,rr,km,kh,ks,buoynl) 

!@sum giss turbulence model for ocean
!@ref Canuto et al. 2004, GRL, 31, L16305 (C2004)
!@+   Canuto et al. 2010, Ocean Modelling, 34, 70-91 (C2010)
!@+   Canuto et al. 2011, Ocean Modelling, 36, 198-207 (C2011)
!@auth AHoward and YCheng
!@date May 2011

      USE GISS_OTURB
      !@var grav acceleration of gravity m/s^2
      !@var lmo max. number of vertical layers
      !@var taubx(im,jm) x component of tau_b = kinematic bottom drag (m/s)^2
      !@var tauby(im,jm) y component of tau_b = kinematic bottom drag (m/s)^2
      !@var rhobot(im,jm) in-situ density at ocean bottom (kg/m^3)	
      !@var exya internal tidal energy (w/m^2)
      !@var ut2a unresolved bottom shear squared (m/s)^2

      implicit none

      ! in:
      integer n          !@var n number of vert. layers on this column
      real*8 ze(0:lmo)   !@var ze vertical grid-edge depth (m), > 0
      real*8 zg(0:lmo+1) !@var zg vertical grid depth (m), < 0
      real*8 db(lmo)     !@var db -grav/rho*d(rho) (m/s^2)
      real*8 dv2(lmo)    !@var dv2 vel. diff. squared btw layers (m/s)^2
      real*8 adt(lmo)    !@var adt rho*alpha*DT (kg/m^3)
      real*8 bds(lmo)    !@var bds rho*beta*DS  (kg/m^3)
      real*8 to(lmo)     !@var to ocean potential tmperature (K)
      real*8 s(lmo)      !@var s ocean salinity (1)
      real*8 rho(lmo)    !@var rho density
      REAL*8 uob         !@var uob x component of velocity at zg(n)
      REAL*8 vob         !@var vob x component of velocity at zg(n)
      REAL*8 u2b         !@var u2b velocity squared at zg(n)
      real*8 fc          !@var fc Coriolis parameter=2*omega*sin(lat) (1/s)
      real*8 ustar       !@var ustar surface friction velocity (m/s)
      real*8 bf          !@var bf surface buoyancy forcing (m^2/s^3)
      real*8 hbl         !@var hbl pbl depth (m)
      integer kbl        !@var kbl first grid level below hbl
      integer strait,ilon,jlat
      intent (in) n,ze,zg,db,dv2,adt,bds,to,s,rho
     &  ,uob,vob,u2b,fc,ustar,bf,hbl,kbl,strait,ilon,jlat

      ! inout:
      REAL*8 e(lmo)      !@var e ocean turbulent kinetic energy (m/s)**2
      intent (inout) e
      ! out:

      real*8 ri(0:lmo+1) !@var ri local richardson number
      real*8 rr(0:lmo+1) !@var rr local alpha*dTdz/(beta*dSdz)
      real*8 km(0:lmo+1) !@var km vertical momentun diffusivity (m**2/s)
      real*8 kh(0:lmo+1) !@var kh vertical heat diffusivity (m**2/s)
      real*8 ks(0:lmo+1) !@var ks vertical salinity diffusivity (m**2/s)
      real*8 buoynl(lmo) !@var buoynl non-local part of buoyancy flux (m^2/s^3)
      intent (out) ri,rr,km,kh,ks,buoynl

      integer :: flag !@var flag =0 if abs(rr)<=1; =1 if abs(rr)>1
c     integer, parameter :: num_smooth=1
      real*8, parameter :: l0min=3.d0,l2min=.05d0
      integer iter, mr
c     real*8 bylenr ! rotation effect
      real*8 vs2(0:lmo+1)!@var vs2 velocity shear squared (1/s**2)
      real*8 bv2(0:lmo+1)!@var bv2 Brunt Vaisala frequency squared (1/s**2)
      real*8 len         !@var len turbulence length scale (m)

      integer l,jlo,jhi,klo,khi
      real*8 a1,a2,b1,b2,c1,c2,c3,c4
      real*8 ril,rrl,gm,sm,sh,ss,kml,khl,ksl,lr,etau
      real*8 l0,l1,l2,kz,zbyh,bydz,zl,tmp,tmp1,tmp2

      ! for nonlocal diffusivities
      real*8, parameter :: ghmin=-5.
      real*8 wstar,wstar3,ustar1,ustar2,ustar3,bylmonin
     &      ,zeta,phi_m,eps,tau,lb,byn,qturb
     &      ,ah,as,am,gh,w2byk,zili,buoy,krl
     &      ,omrr,x,pp,qq,bygam,xx

      ! for background diffusivities
      ! consts appeared in C2010, (65a)-(66)
      !@var f30 2*omega*sin(30 degrees)=omega
      real*8, parameter :: bv0=5.24d-3    ! (1/s), below (65b)
     &   ,cd=3.d-3                        !@var cd dry drag coeff.
     &   ,f30=omega                       ! (1/s), below (65b)
     &   ,bv0byf30=bv0/f30                ! (1), (65b)
     &   ,byden=1./(f30*acosh(bv0byf30))  ! (1), (65b)
     &   ,epsbyn2=.288d-4                 ! (m^2/s), (66)
     &   ,q=.7d0
     &   ,byzet=1./500.d0                   ! (1/m)
      real*8 fbyden,afc,ltn
      real*8 bvbyf,fac
c     real*8 gammbg,smbg,shbg,ssbg
      real*8 kmbg,khbg,ksbg

      ! for tidally-induced diffusivities, exy from table in GISS_OTURB
      real*8 exy,den,fz,epstd_byn2,kmtd,khtd,kstd
      ! for unresolved bottom shear, ut2 from table in GISS_OTURB
      real*8 ut2,ustarb2,phim2,zb,unr20

      !             grid levels                       interface levels
      !
      !                   -----------------------------  surf
      !                1    - - - - - - - - - - - - -
      !                   -----------------------------  1
      !                2    - - - - - - - - - - - - -
      !                   -----------------------------  2
      !               l-1   - - - - - - - - - - - - -
      !                   -----------------------------  l-1
      !zg(l)<0,tracer  l    - - - - - - - - - - - - -
      !                   -----------------------------  l  ze(l)>=0,edge
      !               l+1   - - - - - - - - - - - - -
      !                   -----------------------------  l+1
      !               lm    - - - - - - - - - - - - -
      !                   -----------------------------  lm  
      !             lm+1    - - - - - - - - - - - - -
      !                   -----------------------------  lm+1

      !-----------------------------------------------------------------
      ! more details of giss model are included in subroutine gissdiffus
      ! which is only called from within the subroutine gissmix_init,
      ! the latter is called only once. the lookup uses bisection.
      !-----------------------------------------------------------------

      ! have moved "call gissmix_init" to init_OCEAN in OCNDYN_ye.f

      do l=1,n-1
         bydz=1./(zg(l)-zg(l+1))
         bv2(l)=db(l)*bydz ! N^2
         vs2(l)=dv2(l)*bydz*bydz
      end do
c     do mr = 1,num_smooth
c        call z121(bv2,n-1,lmo)
c        call z121(vs2,n-1,lmo)
c     end do
      do l=1,n-1
         ! ri(l)=db(l)*(zg(l)-zg(l+1))/(dv2(l)+1d-30)
         ri(l)=bv2(l)/(vs2(l)+1d-30)
         rr(l)=bds(l)/(adt(l)+1d-30)
      end do

      ! vertically smooth Rr num_smooth times
c     do mr = 1,num_smooth
c        call z121(rr,n-1,lmo)
c     end do
      l0=.15*hbl

      if(nonlocal.ne.0) then ! nonlocal
         ustar1=max(ustar,3.5d-3)
         ustar2=ustar1*ustar1
         ustar3=ustar1*ustar2
         bylmonin=kappa*bf/ustar3 !@var bylmonin 1/Lmonin
         if(bf.lt.0.) then
            wstar3=-bf*hbl
            wstar=wstar3**by3
            if(nonlocal.eq.2) then
               l=max(kbl,2)
               zili=.2/(1.+7.*(-bf/(hbl*hbl))**(2.*by3)
     &             /max(bv2(l),1d-30))
            endif
         else
            wstar3=0.
            wstar=0.
         endif
      endif

      ! modify ri at the interface nearest to ocean bottom
      ! due to unresolved bottom shear, C2010, eqs,(73)-(77)
      ! using C2002, (35)-(36), find phim(ri), make a table for phim2(ri)
      ! (generalized to include rr dependence)
      ! phim=f/(1-a*f*rf), f=sqrt(2.)/b1*(gm/sm**2)**.25, a=0 or 2.7
      ! iterate for ri=n2/(sig2+unr2)

      if(strait.eq.1) then ! in the strait area
         exy=0.
         ut2=0.
      else                 ! not in the strait area
         exy=exya(ilon,jlat)  ! for tide-induced diffusivities 
         ut2=ut2a(ilon,jlat)  ! for bottom shear effects on Ri 
         ! tidally enhanced bottom drag and bottom density
         ! kinematic drag tau_b given by C2010. eq.(72),
         !@var tau_b momentum flux into rock beneath bottom ocean grid box
         !@+ divided by mass of bottom ocean grid box
         tmp1=cd*(u2b+ut2)**.5d0
         taubx(ilon,jlat)=tmp1*uob
         tauby(ilon,jlat)=tmp1*vob
         rhobot(ilon,jlat)=rho(n)
         ustarb2=cd*(ut2*(u2b+ut2))**.5d0
         l=n-1
         zb=ze(n)-ze(l)
         unr20=ustarb2/(kappa*zb)**2
         rrl=rr(l)
         if(abs(rrl).gt.1.) rrl=1/(rr(l)+1d-30)
         rrl=min(max(rrl,rrmin),rrmax)
         call locate(nt,rra,rrl,klo,khi,a2,b2)
         ril=ri(l)
         do iter=1,10
            ril=min(max(ril,rimin),rimax)
            call locate(mt,ria,ril,jlo,jhi,a1,b1)
            phim2=a2*(a1*phim2a(jlo,klo)+b1*phim2a(jhi,klo))
     &           +b2*(b1*phim2a(jhi,khi)+a1*phim2a(jlo,khi))
            tmp=ril
            ril=bv2(l)/(vs2(l)+unr20*phim2+1d-30)
            ril=.9d0*ril+.1d0*tmp ! best
            if(abs((ril-tmp)/(ril+tmp)).le.1d-2) exit
         end do
         ri(l)=ril
      endif
      
      den=1.-exp(-ze(n)*byzet)
      afc=abs(fc)
      fbyden=afc*byden
c     if(bf.lt.0.) then ! rotation effects
c        bylenr=(afc**3/max(-bf,1d-20))**.5
c     else
c        bylenr=0.
c     endif

      do l=1,n-1

         ril=ri(l)
         rrl=rr(l)
         flag=0
         if(abs(rrl).gt.1.) then
            rrl=1/(rr(l)+1d-30)
            flag=1
         endif

         ! find foreground diffusivities from 2d lookup tables

         ril=min(max(ril,rimin),rimax)
         rrl=min(max(rrl,rrmin),rrmax)
         call locate(mt,ria,ril,jlo,jhi,a1,b1)
         call locate(nt,rra,rrl,klo,khi,a2,b2)
         c1=a1*a2
         c2=b1*a2
         c3=b1*b2
         c4=a1*b2
         gm=c1*gma(jlo,klo)+c2*gma(jhi,klo)
     &     +c3*gma(jhi,khi)+c4*gma(jlo,khi)
         sm=c1*sma(jlo,klo)+c2*sma(jhi,klo)
     &     +c3*sma(jhi,khi)+c4*sma(jlo,khi)
         sh=c1*sha(jlo,klo)+c2*sha(jhi,klo)
     &     +c3*sha(jhi,khi)+c4*sha(jlo,khi)
         ss=c1*ssa(jlo,klo)+c2*ssa(jhi,klo)
     &     +c3*ssa(jhi,khi)+c4*ssa(jlo,khi)
         ! symmetry of giss model: sh <-> ss if rr<->1/rr
         if(flag.eq.1) then
            tmp=sh
            sh=ss
            ss=tmp
         endif

         ! length scale:
         zl=ze(l)
         zbyh=zl/hbl
         kz=kappa*zl
c        lr=c1*lra(jlo,klo)+c2*lra(jhi,klo)
c    &     +c3*lra(jhi,khi)+c4*lra(jlo,khi)
c        if(ril.ge.0.) then 
c           lr=a1*lr1a(jlo)+b1*lr1a(jhi)
c        else
c           lr=1.
c        endif
         lr=1./(1.+max(ril,0.d0))
         if(zl.le.hbl) then         ! within obl
            l1=l0
         else
            l1=l0min+max(l0-l0min,0.d0)*exp(1.-zbyh)
         endif
         l2=max(l1*lr,l2min)
         len=l2*kz/(l2+kz)
         tmp=(osocb1*len)**2*vs2(l)
         e(l)=.5*tmp/(gm+1.d-20)
         e(l)=min(max(e(l),emin),emax)
         etau=.5*osocb1*sqrt(2.*e(l))*len
         kml=etau*sm
         khl=etau*sh
         ksl=etau*ss
         buoynl(l)=0.
         if(nonlocal.ne.0.and.zl.le.hbl) then ! nonlocal, in obl
            if (bv2(l).gt.0.) then
               byn=1./(sqrt(bv2(l))+1d-30)
               qturb=.75*sqrt(max(e(l),emin))
               lb=qturb*byn
            else
               lb=1.d30
            endif
            l1=l1*lb/(l1+lb)
            len=l1*kz/(l1+kz)
            ! find e and update sm,sh,ss within obl
            zeta=zl*bylmonin
            if(zeta.lt.0.) then
               phi_m=(1.-15.*zeta)**(-.25)
            else
               phi_m=1.+4.7*zeta
            endif
            ! eps: modified Moeng and Sullivan 1994
            eps=.4*wstar3/hbl+ustar3*(1.-zbyh)*phi_m/kz
            e(l)=.5*(osocb1*len*eps)**(2.*by3)
            e(l)=min(max(e(l),emin),emax)
            tau=2.*e(l)/(eps+1d-20)
            etau=e(l)*tau
            ! Level 2.5 of OIII model
            rrl=-1.
            gm=tau*tau*vs2(l)
            omrr=1-rrl
            !x=ril*gm/omrr ! x=gh
            x=tau*tau*bv2(l)/omrr ! x=gh
            !if(x.lt.-10.) x=-10. ! if rrl=0
            if(x.lt.ghmin) x=ghmin  ! if rrl=-1
            qq=pi10*(pi20*(1+rrl)-pi30*rrl)
            pp=pi40*(pi50-pi20*(1+rrl))
            bygam=rrl/((pi40/pi10)*(1+qq*x)/(1+pp*x))
            ah=pi40/(1+pp*x+pi20*pi40*x*(1-bygam))
            as=ah/((pi40/pi10)*(1+qq*x)/(1+pp*x))
            xx=(1-bygam)*x*ah
            am=(.8-(pi40-pi10+(pi10-.0066667)*(1-bygam))*x*ah)
     &        /(10.+(pi40-pi10*rrl)*x+.02*gm)
            w2byk=2./(3.+.4*xx+.3*gm*am)
            sh=ah*w2byk
            ss=as*w2byk
            sm=am*w2byk
            kml=max(etau*sm,kml)
            khl=max(etau*sh,khl)
            ksl=max(etau*ss,ksl)
            if(nonlocal.eq.2.and.bf.lt.0.) then ! find cg terms
               buoy=-bf*(1.-(1.+zili)*zbyh)
     &              -ustar3/hbl*zbyh
               krl=(khl-ksl*rrl)/omrr
               buoynl(l)=buoy+krl*bv2(l)
            endif
         endif ! end nonlocal and in obl

         ! background and tidally induced diffusivities

         if(ril.gt.0.) then
            ! C2010, eqs.(65a)-(66); C2011, near end of Sec 4, p.203
            ! afc = abs(Coriol)
            bvbyf=sqrt(max(bv2(l),1d-8))/(afc+1.d-30) 
            if(bvbyf.gt.1.) then
               ltn=acosh(bvbyf)*fbyden  ! dimensionless
               ltn=max(ltn,7.d-2)       ! limit described in C2004
            else
               ltn=7.d-2                ! dimensionless
            endif
            fac=epsbyn2*ltn      ! in m^2/s, Km=GAMMAm*fac

c           ! background diffusivities from 1d lookup tables
c           gammbg=a2*gammbga(klo)+b2*gammbga(khi)
c           smbg=a2*smbga(klo)+b2*smbga(khi)
c           shbg=a2*shbga(klo)+b2*shbga(khi)
c           ssbg=a2*ssbga(klo)+b2*ssbga(khi)
c           ! symmetry of giss model: sh <-> ss if rr<->1/rr
c           if(flag.eq.1) then
c              tmp=shbg
c              shbg=ssbg
c              ssbg=tmp
c           endif
cc          !kmbg=.46d0*fac
c           ! gammam_bg=.5*ribg*gmbg*smbg
c           ! gammah_bg=.5*ribg*gmbg*shbg
            kmbg=fac
            khbg=by3*fac
            ksbg=khbg
c           tmp1=shbg/(smbg+1.d-30)
c           tmp2=ssbg/(smbg+1.d-30)
c           khbg=kmbg*tmp1
c           ksbg=kmbg*tmp2

            ! tidally induced diffusivities, C2010, eqs.(69)-(71)

            ! use the following fz instead of (70) of C10:
            fz=(exp((-zg(l+1)-ze(n))*byzet)
     &         -exp((-zg(l)  -ze(n))*byzet))/(den*(zg(l)-zg(l+1)))
            epstd_byn2=q*exy*fz*2./(rho(l)+rho(l+1))/max(bv2(l),1d-8)
            kmtd=epstd_byn2
            khtd=by3*epstd_byn2
            kstd=khtd
c           khtd=kmtd*tmp1
c           kstd=kmtd*tmp2
         else
            kmbg=0.
            khbg=0.
            ksbg=0.
            kmtd=0.
            khtd=0.
            kstd=0.
         endif      

         ! foreground, background and tidal diffusivities added up
         ! C2010, eqs.(69)-(71)

         km(l)=min(kml+kmbg+kmtd,kmax)
         kh(l)=min(khl+khbg+khtd,kmax)
         ks(l)=min(ksl+ksbg+kstd,kmax)
         
      end do  
      km(0)=0.; kh(0)=0.; ks(0)=0.
      km(1)=max(km(1),kmin);kh(1)=max(kh(1),kmin);ks(1)=max(ks(1),kmin)
      km(n:lmo+1)=0.; kh(n:lmo+1)=0.; ks(n:lmo+1)=0.
      buoynl(n:lmo)=0.
      ri(0)=0.; rr(0)=0.
      ri(n:lmo+1)=0.; rr(n:lmo+1)=0.; e(n:lmo)=emin
      
      return
      end subroutine gissmix


#ifdef OCN_GISSMIX /*alloc_gissmix_com,def_rsf_gissmix,new_io_gissmix*/

      SUBROUTINE alloc_gissmix_com(grid)
!@sum  To allocate arrays who sizes now need to be determined at
!@+    run-time
!@auth Reto Ruedy/Ye Cheng

      USE DOMAIN_DECOMP_1D, only : dist_grid,get
!      USE OCEANR_DIM

      USE GISSMIX_COM

      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(IN) :: grid

      INTEGER :: J_1H, J_0H, L
      INTEGER :: IER

      call getDomainBounds(grid, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)
      !ALLOCATE( otke(LSRPD,IM,J_0H:J_1H) , STAT = IER)
      ALLOCATE( otke(lmo,IM,J_0H:J_1H) , STAT = IER)
      ALLOCATE( rhobot(IM,J_0H:J_1H) , STAT = IER)
      ALLOCATE( taubx(IM,J_0H:J_1H) , STAT = IER)
      ALLOCATE( tauby(IM,J_0H:J_1H) , STAT = IER)
      ALLOCATE( exya(IM,J_0H:J_1H) , STAT = IER)
      ALLOCATE( ut2a(IM,J_0H:J_1H) , STAT = IER)

#ifdef TRACERS_OCEAN
c     ALLOCATE( TRMO1(NTM,IM,J_0H:J_1H),
c    *          TXMO1(NTM,IM,J_0H:J_1H),
c    *          TYMO1(NTM,IM,J_0H:J_1H),
c    *   STAT = IER)
#endif

      END SUBROUTINE alloc_gissmix_com

      subroutine def_rsf_gissmix(fid)
!@sum
!@+    this subroutine is called at the end of
!@+    subroutine def_rsf_ocean, the latter is in OCNDYNtn.f
!@+    remember to put #ifdef OCN_GISSMIX around the call 
!@date Nov 29,2011
      use pario, only : defvar
      use oceanr_dim, only : grid=>ogrid
      use gissmix_com, only : otke
      use straits, only : otkest
      implicit none
      integer fid   !@var fid file id
      call defvar(grid,fid,otke,'otke(lmo,dist_imo,dist_jmo)')
      call defvar(grid,fid,otkest,'otkest(lmo,nmst)')
      return
      end subroutine def_rsf_gissmix

      subroutine new_io_gissmix(fid,iaction)
!@sum
!@+    this subroutine is called at the end of
!@+    subroutine new_io_ocean, the latter is in OCNDYNtn.f
!@+    remember to put #ifdef OCN_GISSMIX around the call 
!@date Nov 29,2011
      use model_com, only : ioread,iowrite
      use pario, only : read_dist_data,write_dist_data
     &                 ,read_data,write_data
      use oceanr_dim, only : grid=>ogrid
      use gissmix_com, only : otke
      use straits, only : otkest
      implicit none
      integer fid     !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
      select case (iaction)
      case (iowrite)            ! output to restart file
        call write_dist_data(grid,fid,'otke',otke,jdim=3)
        call write_data(grid,fid,'otkest',otkest)
      case (ioread)            ! input from restart file
        call read_dist_data(grid,fid,'otke',otke,jdim=3)
        call read_data(grid,fid,'otkest',otkest,bcast_all=.true.)
      end select
      return
      end subroutine new_io_gissmix

#endif /*alloc_gissmix_com,def_rsf_gissmix,new_io_gissmix*/
