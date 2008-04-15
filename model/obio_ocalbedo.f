#include "rundeck_opts.h"      
       subroutine obio_ocalbedo(vrbos)
 
c  Computes ocean surface albedo from solar zenith angle (solz)
c  and wind speed (wind, m/s).
c  Albedo is provided as direct (albd) and diffuse (albs).
  
      USE obio_dim

      USE hycom_dim_glob
      USE hycom_arrays_glob
      implicit none

!!#include "dimensions.h"
#include "dimension2.h"
!!#include "common_blocks.h"

!     parameter(nltgmao=8)
!     real albd(nltgmao),albs(nltgmao)
 
      logical vrbos

c  Derive surface reflectance as a function of solz and wind
c  for OASIM Bands
      call sfcrfl(vrbos)
 
!the rest is not done yet--do be completed later with the coupling?
c  Albedo at GMAO Bands
!     albd = 0.0
!     albs = 0.0
!     albd(1) = rod(1)
!     albs(1) = ros(1)
!     albd(2) = rod(1)
!     albs(2) = ros(1)
!     albd(3) = rod(2)
!     albs(3) = ros(2)
!     albd(4) = rod(2)*0.52
!    *                  + rod(3)
!    *                  + rod(4)
!     albs(4) = ros(2)*0.52
!    *                  + ros(3)
!    *                  + ros(4)
!     albd(4) = albd(4)/2.52
!     albs(4) = albs(4)/2.52
!     n = 0
!     do nl = 5,17
!      albd(5) = albd(5) + rod(nl)
!      albs(5) = albs(5) + ros(nl)
!      n = n+1
!     enddo
!     albd(5) = albd(5)/float(n)
!     albs(5) = albs(5)/float(n)
!     n = 0
!     do nl = 18,23
!      albd(6) = albd(6) + rod(nl)
!      albs(6) = albs(6) + ros(nl)
!      n = n+1
!     enddo
!     albd(6) = albd(6)/float(n)
!     albs(6) = albs(6)/float(n)
!     n = 0
!     do nl = 24,31
!      albd(7) = albd(7) + rod(nl)
!      albs(7) = albs(7) + ros(nl)
!      n = n+1
!     enddo
!     albd(7) = albd(7)/float(n)
!     albs(7) = albs(7)/float(n)
!     n = 0
!     do nl = 32,nlt
!      albd(8) = albd(8) + rod(nl)
!      albs(8) = albs(8) + ros(nl)
!      n = n+1
!     enddo
!     albd(8) = albd(8)/float(n)
!     albs(8) = albs(8)/float(n)
c
      return
      end

c------------------------------------------------------------------
      subroutine sfcrfl(vrbos)
 
c  Computes surface reflectance for direct (rod) and diffuse (ros)
c  components separately, as a function of sunz, wind speed or
c  stress.
c  Includes spectral dependence of foam reflectance derived from Frouin
c  et al., 1996 (JGR)
 
      USE obio_dim
      USE obio_incom,only : rad,rn,wfac,roair
      USE obio_forc, only : sunz,rod,ros,wind
#ifdef OBIO_RAD_coupling
     .          ,chl,obio_bocvn,obio_xocvn
#endif

      USE hycom_dim_glob
      USE hycom_arrays_glob
      implicit none


!!      include 'dimensions.h'
!!      include 'common_blocks.h'

      integer nl
      real cn,rof,rosps,rospd,rtheta
      real sintr,rthetar,rmin,rpls,sinrmin,sinrpls,tanrmin
      real tanrpls,sinp,tanp,a,b

#ifdef OBIO_RAD_coupling
!@sum Those are the weights which were obtained by normalizing solar flux
!@sum for the Lamda's range 0 - 4000 nm. We used Landau fitting function for it.
!@sum They are used to be used for getting 6 band approximation based on 33
!@sum Watson Gregg band calculations
!@sum band_6(j) = Sum(band_33(i)*weight(i))/Sum(part_sum(i))

      real*8 :: sum1, sum2
      real*8 :: refl, ref
      real*8 :: weight(31) = (/0.0158378,0.0201205,0.0241885,0.0277778,
     .                       0.0307124,0.0329082,0.0343586,0.0351143,
     .                       0.0352609,0.0349008,0.0341389,0.0330742,
     .                       0.0317941,0.0303725,0.0288696,0.0273329,
     .                       0.0500921,0.0643897,0.0686573,0.0532013,
     .                       0.0416379,0.0330341,0.0265929,0.0217156,
     .                       0.0179725,0.0150596,0.0127618,0.0158128,
     .                       0.0232875,0.0313132,0.0184843/)
      real*8 :: part_sum(6) = (/0.526854,0.0643897,0.196531,0.066281,
     .                        0.066922,0.0497974/)
#endif
      logical vrbos



c  Foam and diffuse reflectance
      if (wind .gt. 4.0)then
       if (wind .le. 7.0)then
        cn = 6.2E-4 + 1.56E-3/wind
        rof = roair*cn*2.2E-5*wind*wind - 4.0E-4
       else
        cn = 0.49E-3 + 0.065E-3*wind
        rof = (roair*cn*4.5E-5 - 4.0E-5)*wind*wind
       endif
       rosps = 0.057
      else
       rof = 0.0
       rosps = 0.066
      endif
 
c     if(vrbos)write(*,*)'ocalbedo: ',
c    .         nstep,wind,rof,rosps,sunz

c  Direct
c   Fresnel reflectance for sunz < 40, wind < 2 m/s
      if (sunz .lt. 40.0 .or. wind .lt. 2.0)then
       if (sunz .eq. 0.0)then
        rospd = 0.0211
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
       a = 0.0253
       b = -7.14E-4*wind + 0.0618
       rospd = a*exp(b*(sunz-40.0))
      endif

 
c  Reflectance totals
      do nl = 1,nlt
       rod(nl) = rospd + rof*wfac(nl)
       ros(nl) = rosps + rof*wfac(nl)
c     if(vrbos)write(*,*)'ocalbedo: ',
c    .   nstep,nl,rod(nl),ros(nl)
      enddo
 
#ifdef OBIO_RAD_coupling

!  transition between band33 and band6 approximation
      sum1 = 0.0
      sum2 = 0.0
      do nl=1,17 
       sum1 = sum1+weight(nl)*rod(nl)
       sum2 = sum2+weight(nl)*ros(nl)
      enddo

      refl = ref(chl)

!      obio_bocvn(1) = sum2/part_sum(1) + 0.02
!      obio_bocvn(1) = sum2/part_sum(1) + 0.0214577
!      obio_bocvn(1) = sum2/part_sum(1)

       obio_bocvn(1) = sum2/part_sum(1) + refl
       obio_xocvn(1) = sum1/part_sum(1)

!diffuse spectral reflectance average
! c=0.03   0.0404007
! c= 0.1   0.0271808
! c=0.3    0.0195185
! c=1.0    0.015129
! c=3.0    0.0136362
! c=10.0   0.012881
! Average : 0.0214577 

      obio_bocvn(2) = weight(18)*ros(18)/part_sum(2)
      obio_xocvn(2) = weight(18)*rod(18)/part_sum(2)

      sum1 = 0.0
      sum2 = 0.0
      do nl=19,22 
       sum1 = sum1+weight(nl)*rod(nl)
       sum2 = sum2+weight(nl)*ros(nl)
      enddo
      obio_bocvn(3) = sum2/part_sum(3)
      obio_xocvn(3) = sum1/part_sum(3)


      sum1 = 0.0
      sum2 = 0.0
      do nl=23,25
       sum1 = sum1+weight(nl)*rod(nl)
       sum2 = sum2+weight(nl)*ros(nl)
      enddo
      obio_bocvn(4) = sum2/part_sum(4)
      obio_xocvn(4) = sum1/part_sum(4)

      sum1 = 0.0
      sum2 = 0.0
      do nl=26,29
       sum1 = sum1+weight(nl)*rod(nl)
       sum2 = sum2+weight(nl)*ros(nl)
      enddo
      obio_bocvn(5) = sum2/part_sum(5)
      obio_xocvn(5) = sum1/part_sum(5)


      sum1 = 0.0
      sum2 = 0.0
      do nl=30,31 
       sum1 = sum1+weight(nl)*rod(nl)
       sum2 = sum2+weight(nl)*ros(nl)
      enddo
      obio_bocvn(6) = sum2/part_sum(6)
      obio_xocvn(6) = sum1/part_sum(6)

#endif

      return
      end

!_____________________________________________________________
      real*8 function ref(value)

      implicit none

      real*8 :: value

      if (value .le. 0.03)then
        ref = 0.004
      else      
        ref = exp(-4.13846-0.0246239*value) 
     .      + exp(-3.44072-9.54688*value)
      endif

     
      return
      end function ref
!_____________________________________________________________ 
