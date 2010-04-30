#include "rundeck_opts.h"

      subroutine obio_alkalinity(vrbos,kmax,i,j)

!@sum  online computation of alkalinity
!@auth Natassa Romanou

! compute source/sink term for alkalinity
! zc     compensation depth: photosynthesis = respiration
! J_PO4  source/sinks of phosphate
! J_Ca   source/sinks of calcium carbonate
! Jprod  production of organic phosphorus which is linked to particulate
!        organic carbon POC production by 
!        J_POC=(1-sigma_Ca)*J_DOC=(1-sigma_Ca)*rC:P*Jprod
!        to estimate J_DOC use total productivity pp
!        Really, production of particulate organic carbon (or phosphorus)
!        above the compensation depth
! F_Ca   downward flux of CaCO3
! Fc     instantaneous downward flux of particulate organic phosphorus at 
!        compensation depth
! sigma_Ca fraction of P converted to DOP, or fraction of C converted to DOC
! 1-sigma_Ca fraction of P NOT converted to DOP, or fraction of C NOT converted to DOC
! kappa_Ca semi-labile DOP consumption rate constant
! alk tendency is also computed under ice.  this is because 
! alk changes due to sal and temp changes.
! alkalinity units should be umol/kg (uE/kg to go into obio_carbon)
! alk tendency computed here is in umol/m3/s, convert to umol/kg
!-------------------------------------------------------------------------
! J_ALK = -rN:P* J_PO4 + 2*J_Ca
! J_PO4 = J_NO3/rN:P, tendency of phosphates
! J_NO3 : tendency of nitrates, P_tend(:,1)
! J_Ca = R* rC:P* (1-sigma_Ca)* Jprod, for  z<zc
! J_Ca = - dF_Ca/dz,   for z>zc
! Jprod = kappa_Ca*[DOP]-J_PO4, for z<zc 
! Jprod = 0, for z>zc
! [DOP] concentration of dissolved organic phosphorus=[DOC]/rC:P
! [DOC] concentration of dissolved organic C, this is based on phytoplankton
! Distinguish two cases: DOC based on all species
!                        DOC based on coccolithophores only
! F_Ca = R * rC:P * Fc * exp(-(z-zc)/d)
!   Fc = (1-sigma_Ca) * integral_0_zc (J_DOC*dz)
!
! note sign errors in equations 28a and 31 of notes.
!-------------------------------------------------------------------------

      USE obio_dim
      USE obio_incom, only: rain_ratio,cpratio,sigma_Ca,d_Ca,
     .      npratio,uMtomgm3,cnratio,kappa_Ca,bn
      USE obio_com, only: P_tend,p1d,pp2_1d,dp1d,A_tend,
     .      rhs,zc,alk1d,pnoice

#ifdef OBIO_ON_GARYocean
      USE MODEL_COM, only: nstep=> itime
      USE ODIAG, only: ij_fca,oij=>oij_loc
      USE OCEAN, only: dxypo,lmm
#else
      USE hycom_scalars, only: nstep
#endif

      implicit none

      integer nt,k,kmax,nchl1,nchl2,i,j,kzc
      real*8 J_PO4(kmax),pp,Jprod(kmax),Jprod_sum,Fc,zz,F_Ca(kmax+1),
     .       J_Ca(kmax),term,term1,term2,DOP
      logical vrbos
!--------------------------------------------------------------------------
!only compute tendency terms if total depth greater than conpensation depth
      if (p1d(kmax+1) .lt. zc) then
         A_tend(:) = 0.
         go to 100
      endif

!pp and therefore Jprod are defined at dp points (mid level)
!F_Ca is defined at pressure interfaces
!J_Ca is therefore defined mid-layer (dp points)

!compute sources/sinks of phosphate
!J_PO4 units uM/hr
      do k=1,kmax
      J_PO4(k) = 0.1d0 * P_tend(k,1)    !approximate by nitrate conc tendency
                                        !NO3/PO4 ratio from Conkright et al, 1994
      term = -1.d0*npratio * J_PO4(k) 
      rhs(k,15,1) = term
      A_tend(k)= term 
      enddo

      term1=A_tend(1)

!distinguish two cases: OCMIP uses total pp for Jprod, whereas here we also
!consider the case where Jprod only includes coccolithophores
#ifdef Jprod_based_on_pp
      nchl1=1
      nchl2=nchl
#endif
#ifdef Jprod_based_on_cocc
      nchl1=4
      nchl2=4   
#endif

      if (nchl1.le.0 .or. nchl2.le.0) then
          print*, nchl1, nchl2, 
     .   'MUST SET Jprod_based_on_pp or Jprod_based_on_cocc in rundeck'
          stop  
      endif

!compute total primary production
!and the Jprod term
      do k=1,kmax
        if (p1d(k+1) .le. zc)  then
          pp=0.
          do nt=nchl1,nchl2
             !pp2_1d is in mg,C/m3/hr
             !pp is in mg,C/m2/hr     
             pp=pp+pp2_1d(k,nt)/max(p1d(k+1),1.d-3)/bn/cnratio   

cdiag        write(*,'(a,6i5,4e12.4)')'obio_alkalinity, pp:',
cdiag.        nstep,i,j,k,nchl1,nchl2,pp2_1d(k,nt),
cdiag.                    p1d(k+1),bn,cnratio
          enddo
!         DOP = pp/uMtomgm3/cpratio     !convert to uM C
!                                       !Jprod units uM/hr
!         Jprod(k) = kappa_Ca*DOP-J_PO4(k)

          Jprod(k) = pp/uMtomgm3/cpratio     !convert to uM C
                                             !Jprod units uM/hr
        else

          Jprod(k)=0.
        endif
      enddo

!integrate net primary production down to zc
      Jprod_sum = 0.
      do k=1,kmax
      if (p1d(k+1) .le. zc)  then
       Jprod_sum = Jprod_sum + Jprod(k)*dp1d(k)
      endif
      enddo

      Fc = (1.d0-sigma_Ca)* Jprod_sum

!compute downward flux of CaCO3
!only below the euphotic zone (the compensation layer)
      do k=1,kmax+1
         if (p1d(k) .gt. zc)  then
           F_Ca(k) = rain_ratio*cpratio*Fc*exp(-1.d0*(p1d(k)-zc)/d_Ca)
         endif
      enddo

      !find layer index for zc
      do k=kmax+1,1,-1
           if (p1d(k).gt.zc) kzc = k
      enddo
      if (kzc.lt.1) kzc=1
      if (kzc.gt.lmm(i,j)) kzc=lmm(i,j)

      !p1d(kzc) is really the compensation depth
      !F_Ca(kzc) is the CaCO3 export
!     write(*,'(a,i8,3i5,7e12.4)')'CaCO3 downward flux:',
!    .        nstep,i,j,kzc,p1d(kzc),zc,rain_ratio,cpratio,Fc,
!    .        exp(-1.d0*(p1d(kzc)-zc)/d_Ca),F_Ca(kzc)

#ifdef OBIO_ON_GARYocean
      OIJ(I,J,IJ_fca) = OIJ(I,J,IJ_fca) + F_Ca(kzc)
     .                *24.d0*365.d0*1.d-3*12.d0
     .                *dxypo(J)*75.d0*1.d-15
!to convert this into units of Pgr,C/yr: 
! * 24 *365          ! uM/hr -> uM/yr = mili-mol,C/m3/hr
! * 1e-3             ! mol,C/m3/yr
! * 361098046745487. ! mol,C/m/yr, better multiply by dxypo(J)
! * p1d(kzc)         ! mol,C/yr
! * 12.d0            ! gr,C/yr
! * 1.d-15           ! Pgr,C/yr
#endif

!compute sources/sinks of CaCO3
      do k=1,kmax
         if (p1d(k) .le. zc)  then
             !formation of calcium carbonate above compensation depth
             J_Ca(k) = -1.d0*rain_ratio*cpratio*(1.-sigma_Ca)*Jprod(k)
         else
             !dissolution of calcium carbonate below compensation depth
             J_Ca(k) = -1.d0* (F_Ca(k+1)-F_Ca(k)) / dp1d(k)
         endif
       term = 2.d0* J_Ca(k) 
       rhs(k,15,5) = term
       A_tend(k) = A_tend(k) + term      
      enddo

      term2=A_tend(1)-term1

      !for consistency, keep term that goes into rhs table in uM/hr = mili-mol,C/m3/hr
      !convert A_tend terms into uE/kg/hr, the actual units of alkalinity 
      A_tend = A_tend /1024.5d0 *1.d3     ! mili-mol,C/m3/hr -> umol/m3/hr -> umol/kg/hr

! surface boundary condition 
! we probably do not need one here since the effects of buoyancy 
! changes at the surface i.e. the E-P term are included in the 
! total mass changes of the surface model boxes. YES/NO?

!!!!!!!!!! NEED TO ADD BOTTOM BOUNDARY CONDITIONS 

!     if (vrbos) then
!     do k=1,kmax
!     k=1
!     write(*,'(a,4i5,12e12.4)')'obio_alkalinity; ',
!    .    nstep,i,j,k
!    .   ,p1d(k),zc,J_PO4(k),pp,Jprod_sum,Fc
!    .   ,F_Ca(k),J_Ca(k),alk1d(k),A_tend(k),term1,term2
!     enddo
!     endif

 100  continue
      end subroutine obio_alkalinity
