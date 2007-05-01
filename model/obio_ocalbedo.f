       subroutine obio_ocalbedo(vrbos)
 
c  Computes ocean surface albedo from solar zenith angle (solz)
c  and wind speed (wind, m/s).
c  Albedo is provided as direct (albd) and diffuse (albs).
  
      USE obio_dim

      implicit none

#include "dimensions.h"
#include "dimension2.h"
#include "common_blocks.h"

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
      USE obio_incom,only : rad,lam,aw,bw,rn,wfac,roair
      USE obio_forc, only : sunz,rod,ros,wind

      implicit none


      include 'dimensions.h'
      include 'common_blocks.h'

      integer nl
      real cn,rof,rosps,rospd,rtheta
      real sintr,rthetar,rmin,rpls,sinrmin,sinrpls,tanrmin
      real tanrpls,sinp,tanp,a,b

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
 
cdiag if(vrbos)write(888,'(i5,4e12.4)')
cdiag.         nstep,wind,rof,rosps,sunz

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
cdiag if(vrbos)write(889,'(i5,2e12.4)')nstep,rod(nl),ros(nl)
      enddo
 
      return
      end
