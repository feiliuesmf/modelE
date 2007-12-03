      subroutine obio_daysetrad
c
c  Sets daily parameters for ocean irradiance.
c
      USE obio_dim
      USE obio_incom, only : lam,ac,bc,aw,nl450,excdom
      USE obio_com, only : npst,npnd,WtoQ,obio_P,avgq1d,ihra_ij
     .                    ,acdom

      implicit none

#include "dimensions.h"
#include "dimension2.h"
#include "common_blocks.h"


      integer :: nl,nt
      real :: actot450,atot450

c  Compute acdom
      if (nl450.eq.0) stop 'obio_daysetrad: nl450=0'
      !!m = indext2
      do k = 1,kdm
        actot450 = 0.0
        atot450  = 0.0
        do nt = 1,nchl
         actot450 = actot450  + obio_P(k,nnut+nt)*ac(nt,nl450)
        enddo
        atot450 = aw(nl450) + actot450
        do nl = npst,npnd
         acdom(k,nl) = 0.2*atot450*excdom(nl)
        enddo
      enddo
 
c  Compute average quanta; Initialize light history arrays
c    (ihra set to 1 first time thru to preserve restart file values)
      do k = 1,kdm
       if (ihra_ij .gt. 0)then
        avgq1d(k) = avgq1d(k)/float(ihra_ij)
       endif
      enddo

      return
      end
