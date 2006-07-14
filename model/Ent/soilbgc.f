      module soilbgc

!@sum Routines to simulate soil biogeochemistry:
!@sum microbial dynamics, C & N pools, respirationa and N fluxes.


      use ent_const
      use ent_types
      use ent_pfts
      implicit none


      contains
      !*********************************************************************
      
      subroutine soil_bgc(dtsec, pp)
      real*8, intent(in) :: dtsec  !main ent time step (s)
!      real*8 livefr(nlive)   !fractions of plant C
      real*8 soilmoist       !soil moisture avg over top 30 cm
      real*8 fnpp            !npp (input kgC/m2/d --> convert to gC/m2/s) 
      integer ivt            !ivt = pft (see above)
      real*8 lai             !total leaf area, unadjusted for burying by snow
      real*8 soiltemp        !soil temperature avg over top 30 cm
      real*8 clayfrac        !fractional clay content in soil
      real*8 sandfrac        !fractional clay content in soil
      real*8 siltfrac        !fractional clay content in soil
      real*8 Tpool(ptrace,NPOOLS)    !total plant and soil C,N pools (gC/m2)
      real*8 Cflux           !total soil C flux to atm (gC/m2/s) **NOT CO2!**
      real*8 fco2                    !NEE = net ecosys co2 flux (umol co2/m2 /s) [+ = to atm]
      real*8 soil_resp                  !soil resp flux (soil CO2 flux, umol co2/m2/s)
      type(entcelltype),pointer :: ecp
      type(patch),pointer :: pp
      type(cohort),pointer :: cop, scop

      scop => pp%sumcohort
      ! pass NPP and LAI to CASA as fnpp, lai (using patch summ. values) 
      lai = scop%LAI                !total lai
      fnpp = scop%NPP*1.d3/sday     !total npp !convert npp from kgC/m2/d to gC/m2/s -PK

      ivt = scop%pft     
      soilmoist = pp%Soilmoist    !soil moist varies by patch
      soiltemp = pp%cellptr%Soiltemp  !soil temp, texture vary by cell
      clayfrac = pp%cellptr%soil_texture(3) !in GHY.f, texture order is sand,loam,clay,peat(+bedrock) -PK 7/13/06
      sandfrac = pp%cellptr%soil_texture(1)
! use siltfrac = 1 - (clayfrac + sandfrac) or 0.4*"loam" for now -PK 6/14/06
      siltfrac = 0.4*pp%cellptr%soil_texture(2)
      
      print *,'casa inputs: dt= ',dtsec,'soilmoist=',soilmoist     !-PK 7/13/06
     &       ,'npp=',fnpp,'ivt=',ivt,'lai=',lai,'soiltemp=',soiltemp
     &       ,'clay=',clayfrac,'sand=',sandfrac,'silt=',siltfrac
     &       ,'Tpool(in)=',Tpool
     
!      print*, 'calling casa (soil_bgc)'
      call casa (dtsec, soilmoist, fnpp, ivt
     &                 ,lai, soiltemp, clayfrac, sandfrac, siltfrac  
     &                 ,Tpool,Cflux,fco2)
      !convert Cflux from gC/m2/s to umol C/m2/s -PK 6/14/06
      soil_resp = Cflux*1.d6/12.  
      !assign soil_resp per patch
      pp%soil_resp = soil_resp   
      !assign Tpool, CO2flux per patch to Tpool, fco2 from casa -PK 7/7/06
      pp%Tpool = Tpool
      pp%CO2flux = fco2 
      
      print *,'casa outputs: soil_resp(~Cflux)=',soil_resp,'Tpool(out)='
     &       ,Tpool,'CO2flux(~fco2)=',fco2 

      end subroutine soil_bgc


      !*********************************************************************
  

      end module soilbgc
