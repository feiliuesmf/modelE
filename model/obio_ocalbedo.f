#include "rundeck_opts.h"      

      subroutine obio_ocalbedo(wind,solz,bocvn,xocvn,chl,vrbos)
c  Computes ocean surface albedo from solar zenith angle (solz)
c  and wind speed (wind, m/s).
c  Albedo is provided as direct (albd) and diffuse (albs).
c  Derive surface reflectance as a function of solz and wind
c  Includes spectral dependence of foam reflectance derived from Frouin
c  et al., 1996 (JGR)
      USE FILEMANAGER
      implicit none

      integer nl
      real*8 cn,rof,rosps,rospd,rtheta
      real*8 sintr,rthetar,rmin,rpls,sinrmin,sinrpls,tanrmin
      real*8 tanrpls,sinp,tanp,a,b

      real*8, intent(in)  :: wind, solz, chl
      real*8 :: sunz
 
      real*8, dimension(6), intent(out) :: bocvn
      real*8, dimension(6), intent(out) :: xocvn

      real*8 :: sum1, sum2, part_sum
      real*8 :: refl, ref

!@sum Those are the weights which were obtained by normalizing solar flux
!@sum for the Lamda's range 0 - 4000 nm. We used Landau fitting function for it.
!@sum They are used to be used for getting 6 band approximation based on 33
!@sum Watson Gregg band calculations
!@sum band_6(j) = Sum(band_33(i)*weight(i))/Sum(part_sum(i))
C**** WHY IS WEIGHT ONLY DECLARED TO BE 31 AND NOT 33?
      real*8 :: weight(31) = (/0.0158378,0.0201205,0.0241885,0.0277778,
     .                       0.0307124,0.0329082,0.0343586,0.0351143,
     .                       0.0352609,0.0349008,0.0341389,0.0330742,
     .                       0.0317941,0.0303725,0.0288696,0.0273329,
     .                       0.0500921,0.0643897,0.0686573,0.0532013,
     .                       0.0416379,0.0330341,0.0265929,0.0217156,
     .                       0.0179725,0.0150596,0.0127618,0.0158128,
     .                       0.0232875,0.0313132,0.0184843/)

C**** gband is the distribution of the 31 bands for the 6 band GISS code
      integer :: gband(7) = (/ 1, 18, 19, 23, 26, 30, 32 /)
c      real*8 :: part_sum(6) = (/0.526854,0.0643897,0.196531,0.066281,
c     .                        0.066922,0.0497974/)
      logical vrbos

C**** this should be initialised somewhere else, but is here for now to
C**** decouple from obio.

      integer, parameter :: nlt=33
      real*8 :: pi, rad, roair, rn, rod(nlt),ros(nlt), wfac(nlt),
     *     aw(nlt), bw(nlt), saw, sbw
      real*8 :: b0, b1, b2, b3, a0, a1, a2, a3, t, tlog, fac, rlam
      integer :: ic , iu_bio, ngiss, lambda, lam(nlt)
      character title*50
      data a0,a1,a2,a3 /0.9976d0, 0.2194d0,  5.554d-2,  6.7d-3 /
      data b0,b1,b2,b3 /5.026d0, -0.01138d0, 9.552d-6, -2.698d-9/

      pi = dacos(-1.0D0)
      rad = 180.0D0/pi
      rn = 1.341d0  ! index of refraction of pure seawater
      roair = 1.2D3 ! density of air g/m3  SHOULD BE INTERACTIVE?

      call openunit('cfle1',iu_bio)
      do ic = 1,6
        read(iu_bio,'(a50)')title
      enddo
      do nl = 1,nlt
        read(iu_bio,20) lambda,saw,sbw
        lam(nl) = lambda
        aw(nl) = saw
        bw(nl) = sbw
        if (lam(nl) .lt. 900) then
          t = exp(-(aw(nl)+0.5*bw(nl)))
          tlog = alog(1.0D-36+t)
          fac = a0 + a1*tlog + a2*tlog*tlog + a3*tlog*tlog*tlog
          wfac(nl) = max(0d0,min(fac,1d0))
        else
          rlam = float(lam(nl))
          fac = b0 + b1*rlam + b2*rlam*rlam + b3*rlam*rlam*rlam
          wfac(nl) = max(fac,0d0)
        endif
      enddo
      call closeunit(iu_bio)
 20   format(i5,f15.4,f10.4)

C**** end initialisation
      
      sunz=acos(solz)*rad  !in degs

c  Foam and diffuse reflectance
      if (wind .gt. 4.0) then
        if (wind .le. 7.0) then
          cn = 6.2D-4 + 1.56D-3/wind
          rof = roair*cn*2.2D-5*wind*wind - 4.0D-4
        else
          cn = 0.49D-3 + 0.065D-3*wind
          rof = (roair*cn*4.5D-5 - 4.0D-5)*wind*wind
        endif
        rosps = 0.057d0
      else
        rof = 0.0
        rosps = 0.066d0
      endif
      
      if(vrbos)write(*,*)'ocalbedo: ', wind,rof,rosps,sunz

c  Direct
c   Fresnel reflectance for sunz < 40, wind < 2 m/s
      if (sunz .lt. 40.0 .or. wind .lt. 2.0) then
        if (sunz .eq. 0.0) then
          rospd = 0.0211d0
        else
          rtheta = sunz/rad
          sintr = sin(rtheta)/rn
          rthetar = asin(sintr)
          rmin = rtheta - rthetar
          rpls = rtheta + rthetar
          sinrmin = sin(rmin)
          sinrpls = sin(rpls)
          tanrmin = tan(rmin)
          tanrpls = tan(rpls)
          sinp = (sinrmin*sinrmin)/(sinrpls*sinrpls)
          tanp = (tanrmin*tanrmin)/(tanrpls*tanrpls)
          rospd = 0.5*(sinp + tanp)
        endif
      else
       !Empirical fit otherwise
        a = 0.0253d0
        b = -7.14D-4*wind + 0.0618d0
        rospd = a*exp(b*(sunz-40.0))
      endif
      if(vrbos)write(*,*)'ocalbedo: ', rospd

c  Reflectance totals
      do nl = 1,nlt
        rod(nl) = rospd + rof*wfac(nl)
        ros(nl) = rosps + rof*wfac(nl)
        if(vrbos)write(*,*)'ocalbedo: ', nl,rod(nl),ros(nl),rof,wfac(nl)
      enddo
      
!diffuse spectral reflectance average
! c=0.03   0.0404007
! c=0.1   0.0271808
! c=0.3    0.0195185
! c=1.0    0.015129
! c=3.0    0.0136362
! c=10.0   0.012881
! Average : 0.0214577 

!  transition between band33 and band6 approximation

! loop over giss radiation bands
      do ngiss=1,6
        sum1 = 0.0
        sum2 = 0.0
        do nl=gband(ngiss), gband(ngiss+1)-1 
          sum1 = sum1+weight(nl)*rod(nl)
          sum2 = sum2+weight(nl)*ros(nl)
        enddo
        part_sum=sum(weight(gband(ngiss):gband(ngiss+1)-1))
        xocvn(ngiss) = sum1/part_sum
        bocvn(ngiss) = sum2/part_sum
      end do

C**** add chlorophyll term
      refl = ref(chl)
      bocvn(1) = bocvn(1) + refl
      if (vrbos) print*,"ocalbedo:",refl,chl

      return
      end subroutine obio_ocalbedo


      real*8 function ref(chl)
!@sum calculate reflectance as a function of chlorophyll calculation
      implicit none
      real*8, intent(in) :: chl  !@var chl Chlorohpyll concentration

C**** formula from ??????
      if (chl .le. 0.03) then 
        ref = 0.004d0
      else      
        ref = exp(-4.13846d0-0.0246239d0*chl) 
     .       + exp(-3.44072d0-9.54688d0*chl)
      endif
     
      return
      end function ref

