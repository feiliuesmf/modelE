#include "rundeck_opts.h"

      subroutine calc_flammability(t,p,r,v,flam)
!@sum calculated the flammability of vegetation based on model
!@+ variables for temperature, precipitation, relative humidity,
!@+ and an index of vegetation density read from an input file.
!
!@auth Greg Faluvegi based on information from Olga Pechony
!@ver  1.0 (based on Olga's Word document Flammability.doc)
!
!@var T the surface air temperature passed in Kelvin for current I,J
!@var P the precipitation rate passed in mm/day for current I,J
!@var R the relative humidity passed as fraction for current I,J
!@var V vegetative density (unitless 0 to 1) for current I,J
!@var Z component of the Goff-Gratch saturation vapor pressure equation
!@var tsbyt = reciprocal of the temperature times ts
!@var flam the flammability returned for current I,J
!@param a,b,s,d,f,h,ts,cr coefficients for the parameterization
!@+ e.g. from Goff-Gratch saturation vapor pressure equation
!
      implicit none

      real*8, parameter :: a=-7.90298d0,d=11.344d0,c=-1.3816d-7,
     & b=5.02808,f=8.1328d-3,h=-3.49149d0,ts=373.16d0,cr=-2.d0
      real*8, intent(in) :: t,p,r,v
      real*8, intent(out) :: flam
      real*8 :: z,tsbyt

      tsbyt=ts/t

      z= a*(tsbyt-1.d0) + b*log10(tsbyt) + 
     &   c*(10.d0**(d*(1.d0-tsbyt))-1.d0) +
     &   f*(10.d0**(h*(tsbyt-1.d0))-1.d0)

      flam=(10.d0**z*(1.d0-r)) * exp(cr*p) * v

      return
      end subroutine calc_flammability


      subroutine prec_running_average(p,avg,iH,iD,i0,first,HRA,DRA,PRS)
!@sum prec_running_average keeps a running average of the model
!@+ precipitation variable for use in the flammability model.
!@+ In practive, this does hourly and daily running averages and
!@+ uses those to get the period-long running average, to avoid
!@+ saving a huge array.
!@auth Greg Faluvegi
!@ver 1.0
      use flammability_com, only: nday=>nday_prec,nmax=>maxHR_prec
 
      implicit none
      
!@var p model variable which will be used in the average (precip)
!@var temp just for holding current day average for use in avg_prec
!@var nmax number of accumulations in one day (for checking)
!@var bynmax reciprocal of nmax
!@var bynday reciprocal of nday
!@var i0 see i0fl(i,j) (day in period pointer)
!@var iH see iHfl(i,j) (hour in day index)
!@var iD see iDfl(i,j) (day in period index)
!@var first see first_prec(i,j) whether in first period
!@var HRA see HRAfl(i,j,time) (hourly average)
!@var DRA see DRAfl(i,j,time) (daily average)
!@var PRS see PRSfl(i,j) (period running sum)
!@var avg see ravg_prec(i,j) the running average prec returned
      real*8, intent(IN) :: p
      real*8, dimension(nday) :: DRA
      real*8, dimension(nmax) :: HRA
      real*8 :: temp, bynmax, PRS, avg, iH, iD, i0, first, bynday
      integer :: n

      if(nint(iH) < 0 .or. nint(iH) > nmax) then
        write(6,*) 'iH maxHR_prec=',iH,nint(iH),nmax
        call stop_model('iHfl or maxHR_prec problem',255)
      endif
      bynmax=1.d0/real(nmax)
      bynday=1.d0/real(nday)
      iH = iH + 1.d0
      HRA(nint(iH)) = p
      ! do no more, unless it is the end of the day:

      if(nint(iH) == nmax) then ! end of "day":
        iH = 0.d0
        if(nint(first) == 1) then ! first averaging period only
          iD = iD + 1.d0
          do n=1,nmax
            DRA(nint(iD)) = DRA(nint(iD)) + HRA(n)
          end do
          DRA(nint(iD)) = DRA(nint(iD))*bynmax
          if(nint(iD) == nday) then ! end first period
            PRS = 0.d0
            do n=1,nday
              PRS = PRS + DRA(n)
            end do
            avg = PRS * bynday
            first=0.d0
            iD=0.d0
            i0=0.d0
          end if
        else ! not first averaging period: update the running average
          i0 = i0 + 1.d0 ! move pointer
          if(nint(i0) == nday+1) i0=1.d0 ! reset pointer
          temp=0.d0
          do n=1,nmax
            temp = temp + HRA(n)
          end do
          temp = temp * bynmax ! i.e. today's average
          PRS = PRS - DRA(nint(i0))
          DRA(nint(i0)) = temp
          PRS = PRS + DRA(nint(i0))
          avg = PRS * bynday
        end if
      end if

      end subroutine prec_running_average


#if defined DYNAMIC_BIOMASS_BURNING && defined CALCULATE_FLAMMABILITY
      subroutine dynamic_biomass_burning(n,ns)
!@sum dynamic_biomass_burning fills in the surface source ns for
!@+ tracer n with biomass burning based on flammability and offline
!@+ correlations of emissions with observed fire counts.
!@auth Greg Faluvegi based on direction from Olga Pechony
!@ver 1.0

      use domain_decomp_atm,only: grid, get
      use flammability_com, only: epfc,mfcc,flammability,first_prec,
     &                            missing
      use tracer_com, only: sfc_src
      use model_com, only: nday,DTsrc
      use ghy_com, only: fearth

      implicit none
   
      integer :: J_0S, J_1S, I_0H, I_1H, i, j
      integer, intent(in) :: n,ns

      call get(grid, J_STRT_SKP=J_0S, J_STOP_SKP=J_1S,
     &              I_STRT_HALO=I_0H,I_STOP_HALO=I_1H)

      ! units: sfc_src = kg/m2/s
      ! units: epfc = kg/m2/#fire
      ! units: flammability*mfcc = #fire/day
      ! units: DTsrc*nday = s/day

      do j=J_0S,J_1S
        do i=I_0H,I_1H
          if(fearth(i,j) > 0. .and. first_prec(i,j)==0 
     &      .and. flammability(i,j)/=missing)then
            sfc_src(i,j,n,ns)=
     &      epfc(i,j,n)*flammability(i,j)*mfcc/(nday*DTsrc)
          else
            sfc_src(i,j,n,ns)=0.d0
          endif
        enddo
      enddo
    
      end subroutine dynamic_biomass_burning
#endif

