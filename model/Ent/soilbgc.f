      module soilbgc

!@sum Routines to simulate soil biogeochemistry:
!@sum microbial dynamics, C & N pools, respiration and N fluxes.

      use ent_const
      use ent_types
      use ent_pfts
      implicit none

      public soil_bgc

      contains
      
!***********************************************************************      
      subroutine soil_bgc(dtsec, pp)
      use patches, only : print_Tpool
      real*8, intent(in) :: dtsec  !main ent time step (s)
      type(patch),pointer :: pp
      !----Local----------
      real*8 soilmoist       !soil moisture avg over top 30 cm
      integer ivt            !ivt = pft (see above)
      real*8 soiltemp        !soil temperature avg over top 30 cm
      real*8 clayfrac        !fractional clay content in soil
      real*8 sandfrac        !fractional clay content in soil
      real*8 siltfrac        !fractional clay content in soil
      real*8 Tpool(PTRACE,NPOOLS)    !total plant and soil C,N pools (gC/m2)
      real*8 Cflux           !total soil C flux to atm (gC/m2/s) **C flux, NOT CO2!**

      type(cohort),pointer :: scop

      scop => pp%sumcohort
      
      ivt = scop%pft     
      soilmoist = pp%Soilmoist    !soil moist varies by patch
      soiltemp = pp%cellptr%Soiltemp  !soil temp, texture vary by cell
      clayfrac = pp%cellptr%soil_texture(3) !in GHY.f, texture order is sand,loam,clay,peat(+bedrock) -PK 7/13/06
      sandfrac = pp%cellptr%soil_texture(1)
! use siltfrac = 1 - (clayfrac + sandfrac) or 0.4*"loam" for now -PK 6/14/06
      siltfrac = 0.4*pp%cellptr%soil_texture(2)  !**hack** -PK
      Tpool = pp%Tpool !Added - NYK 7/27/06
      
#ifdef DEBUG
      print *,'casa inputs: dt= ',dtsec,'soilmoist=',soilmoist     
     &       ,'ivt=',ivt,'soiltemp=',soiltemp 
     &       ,'clay=',clayfrac,'sand=',sandfrac,'silt=',siltfrac
      call print_Tpool(Tpool)
#endif
     
      call casa_bgfluxes(dtsec, soilmoist,  ivt    
     &                 , soiltemp, clayfrac, sandfrac, siltfrac
     &                 , Tpool, Cflux)
     
      !* Convert Cflux from gC/m2/s to kgC/m2/s 
      !* and assign patch soil_resp, Tpool 
      pp%Soil_resp = Cflux*1.d-3 
      pp%Tpool = Tpool

#ifdef DEBUG
      call print_Tpool(Tpool)
      print *,'casa outputs: soil_resp(kgC/m2/s)=',pp%soil_resp
#endif

      end subroutine soil_bgc

!***********************************************************************
      subroutine casa_bgfluxes(dtsec, soilmoist, ivt
     &                 ,soiltemp, clayfrac, sandfrac, siltfrac 
     &                 ,Tpool,Cflux)

      implicit none

! ------------------------ input/output variables -----------------
! input 
      real*8,intent(in) :: dtsec           !main ent time step (s)
      real*8,intent(in) :: soilmoist       !soil moisture avg over top 30 cm
      integer,intent(in) :: ivt            !ivt = pp%sumcohort%pft (assigned in soilbgc)
      real*8,intent(in) :: soiltemp        !soil temperature avg over top 30 cm
      real*8,intent(in) :: clayfrac        !fractional clay content in soil
      real*8,intent(in) :: sandfrac        !fractional sand content in soil
      real*8,intent(in) :: siltfrac        !fractional silt content in soil
      
! i/o
      real*8,intent(inout) :: Tpool(PTRACE,NPOOLS)  !total plant and soil C,N pools

! output
      real*8,intent(out) :: Cflux      !total respiration flux to atm (gC/m2/s) !main output from this routine -PK 6/15/06

! ------------------------ local variables ------------------------
      integer n,m
      integer ipool, ifirst
      real*8 Cdead
!decomp coefs      
      real*8 fact_soilmic(N_PFT)
      real*8 fact_slow(N_PFT)
      real*8 fact_passive(N_PFT)
      real*8 eff(NRESP_POOLS)          
      real*8 frac_donor(NRESP_POOLS)   
      real*8 kdt(N_PFT,NPOOLS)
!met variables and functions
      real*8 atmp, bgmoist, bgtemp
      real*8 Wlim    !water limitation function for CASA (func of next 2)
      real*8 watopt  !"optimal" water content for ET
      real*8 watdry  !water content when ET stops (~wilting point)
      !watopt, watdry are funcs of next 3 (which are funcs of clayfrac,sandfrac -- see lsmtci.F) 
      real*8 watsat  !saturated volumetric soil water content (porosity)
      real*8 smpsat  !soil matric potential at saturation (mm)
      real*8 bch     !clapp and hornberger "b"

      real*8 Closs(PTRACE,NPOOLS)    !amt C lost per pool (gC/m2)
      real*8 Resp(PTRACE,NPOOLS)     !amt C lost to atm per pool (gC/m2), used to calculate Cflux
      real*8 poolsum                 !accumulates Resp

      
! define rate coefs fact_soilmic, fact_slow, fact_passive (from casatci.F) -PK 5/25/06
        if(ivt.eq.CROPS)then      
          fact_soilmic(ivt) = 1.25d0
          fact_slow(ivt)    = 1.5d0
          fact_passive(ivt) = 1.5d0
        else
          fact_soilmic(ivt) = 1.d0
          fact_slow(ivt)    = 1.d0
          fact_passive(ivt) = 1.d0
        end if
        if(ivt.eq.0)then
          fact_soilmic(ivt) = 0.d0
          fact_slow(ivt)    = 0.d0
          fact_passive(ivt) = 0.d0
        end if

* Maximum RATE CONSTANTS FOR EACH POOL SCALED TO LENGTH OF TIME STEP
* For small delta_t, kdt   is the same as annK*delta_t
*iyf:  Consider dM/dt =  - M/tau
*iyf:  Analytic solution:  M(t) = M0 exp (-t/tau)
*iyf:  Integrate Flux*dt from t=(n-1)*dt to t=n*dt:
*iyf:  integral = M(n*dt) - M[(n-1)*dt]
*iyf:           = - M[(n-1)*dt] {1 - exp [-dt/tau]}
*iyf:    approx = M[(n-1)*dt] * [dt/tau]
*iyf:  variable kdt was previously Krate in CASA

** NOTE: kdt will be used for WOOD and dead pools only
**       so no need to worry about adding stressCD to annK here
      do n = 1, NPOOLS
        do m = 1, N_PFT    
          kdt(m,n)= 1.d0 - (exp(-annK(m,n)*dtsec))
        enddo
      enddo  


*---Step 1: heterotrophic (microbial) respiration -------------------------
C.. Initialize respiration fluxes each timestep
C.. Note: these should be over dead pools only (see resp_pool_index)
*iyf:  Resp in unit of gC/m2/timestep

      Resp(:,:) = 0.d0

*---Step 1a: TEMPERATURE AND MOISTURE CONSTRAINTS ON DECOMP 
! temperature dependence
            bgtemp = (Q10 ** ((soiltemp - 30.d0) / 10.d0))

! moisture dependence 
* mimic calculation of bevap in surphy.F to get Wlim
* but use smoist,soiltemp instead of h2osoi,tsoi 
*   watdry = water content when evapotranspiration stops = wp
            watsat = 0.489d0 - 0.00126d0*sandfrac 
            smpsat = -10.d0 * ( 10.d0**(1.88d0-0.0131d0*sandfrac) )
            bch = 2.91d0 + 0.159d0*clayfrac
            watdry = watsat * (-316230.d0/smpsat) ** (-1.d0/bch)
            watopt = watsat * (-158490.d0/smpsat) ** (-1.d0/bch)
!! 03/11/21 note: there are no limits on Wlim, except if soiltemp < 0
!!Wlim ultimately used to get total C loss per pool,  
!!to both atm and other pools (see step 1b below) -PK 6/8/06
            if (soiltemp .gt. 0.d0) then
               Wlim = min( max(soilmoist-watdry,0.d0) /
     &                   (watopt-watdry), 1.d0)
            else
               Wlim = 0.01d0
            end if

            bgmoist = 0.25d0 + 0.75d0*Wlim
            atmp = bgtemp * bgmoist

*---Step 1b: DETERMINE loss of C FROM EACH DEAD POOL (donor) PER TIMESTEP
*iyf:  Closs is the amount of carbon each pool loses in gC/m2/timestep.
*iyf:  A fraction of Closs is transferred to another pool (receiver pool), 
*iyf:  the remainder to the atm (Resp).
*iyf:  The distribution of Closs is done in subroutine casa_respire, Step 1c. 
            do n = 1, NDEAD
               ipool = NLIVE + n
               Cdead   = Tpool(Carbon,ipool) * kdt(ivt,ipool) * atmp
               Closs(Carbon,ipool) = Cdead
            enddo

** adjust pools
            Closs(Carbon,SURFSTR) = Closs(Carbon,SURFSTR)
     &                  * lignineffect(ivt)
            Closs(Carbon,SOILSTR) = Closs(Carbon,SOILSTR)
     &                  * lignineffect(ivt)
            Closs(Carbon,SOILMIC) = Closs(Carbon,SOILMIC)
     &                  * (1.d0-(0.75d0*(siltfrac+clayfrac)))  !see Potter et al 1993 for ref -PK 6/8/06
     &                  * fact_soilmic(ivt)
            Closs(Carbon,SLOW) = Closs(Carbon,SLOW)
     &                  * fact_slow(ivt)
            Closs(Carbon,PASSIVE) = Closs(Carbon,PASSIVE)
     &                  * fact_passive(ivt)

*iyf:  limits on loss from dead pools.
ciyf - No need to track limits on loss rate
ciyf (only need to track limits on inventories)
            do n = NLIVE+1,NPOOLS
              Closs(Carbon,n)=MIN(Closs(Carbon,n),Tpool(Carbon,n))
            enddo     
 
*---Step 1c:  SOM C DECOMPOSITION
*iyf:  update Tpool's:  C is transferred between donor and receiver pools 
*iyf:  Resp is the amount of C to atm, in gC/m2/timestep

* MICROBIAL EFFICIENCIES FOR PARTICULAR FLOWS
        eff( 1) =  0.45    ! SLOW,PASSIVE
        eff( 2) =  0.45    ! SLOW,SOILMIC
        eff( 3) =  0.40    ! SURFMET,SURFMIC
        eff( 4) =  0.40    ! SURFSTR,SURFMIC
        eff( 5) =  0.70    ! SURFSTR,SLOW
        eff( 6) =  0.45    ! SOILMET,SOILMIC
        eff( 7) =  0.45    ! SOILSTR,SOILMIC
        eff( 8) =  0.70    ! SOILSTR,SLOW
        eff( 9) =  0.40    ! CWD,SURFMIC
        eff(10) =  0.70    ! CWD,SLOW
        eff(11) =  0.40    ! SURFMIC,SLOW
        eff(12) =  0.85 - (0.68 * (siltfrac+clayfrac))  ! SOILMIC,PASSIVE 
        eff(13) =  0.85 - (0.68 * (siltfrac+clayfrac))  ! SOILMIC,SLOW 
        eff(14) =  0.45    ! PASSIVE,SOILMIC
* EXTRA RESPIRATION TRANSFER EFFICIENCIES
        frac_donor( 1) =  0.003 + (0.009*clayfrac)
        frac_donor( 2) =  1.0 - frac_donor(1)
        frac_donor( 3) =  1.0
        frac_donor( 4) =  1.0 - structurallignin(ivt)
        frac_donor( 5) =  structurallignin(ivt)
        frac_donor( 6) =  1.0
        frac_donor( 7) =  1.0 - structurallignin(ivt)
        frac_donor( 8) =  structurallignin(ivt)
        frac_donor( 9) =  1.0 - woodligninfract
        frac_donor(10) =  woodligninfract
        frac_donor(11) =  1.0
        frac_donor(12) =  0.003 + (0.032*clayfrac)
        frac_donor(13) =  1.0 - frac_donor(12)
        frac_donor(14) =  1.0

! casa_respire determines amount of C transferred from each "donor" pool and 
! where it goes ("receiver" pool), in addition to total C respired (Resp-->Cflux) -PK 6/15/06
      call casa_respire(ivt ,eff ,frac_donor
     &                 ,Closs, Resp, Tpool)
     
*---Step 2-----------------------------------------------------------
* CALCULATE NITROGEN POOLS !***keep this for future use*** -PK 8/23/06                                 
!            do n = 1,NPOOLS
!              Tpool(Nitrogen,n)=Tpool(Carbon,n)/CNratio(n)
!            enddo   

*--- Step 3 --------------------------------------------------------
*  Get Total C Fluxes to atm = Sum over respiring (dead) pools/dtsec
*  Note: Resp for live pools should be zero
*iyf:  Cflux in gC/m2/s
           poolsum = 0.0
           do n = NLIVE+1,NPOOLS
              poolsum=poolsum+Resp(Carbon,n)
           enddo  !over dead pools -PK
           Cflux=poolsum/dtsec           !total C flux to atm (gC/m2/s)

C.. mod 03/03/11 change closs/nloss to be g/m2/sec instead of g/m2/timestep
           do n = 1,NPOOLS
              Closs(Carbon,n)=Closs(Carbon,n)/dtsec
!             Closs(Nitrogen,n)=Closs(Nitrogen,n)/dtsec  !not used (yet) -PK  
           enddo  !over all pools

      return
      end subroutine casa_bgfluxes
      
!***********************************************************************
      subroutine casa_respire(ivt ,eff ,frac_donor
     &                       ,Closs, Resp, Tpool)

      implicit none

! ------------------------ input/output variables -----------------
! input
      integer,intent(in) :: ivt            !ivt = pp%sumcohort%pft (assigned in soilbgc)
      real*8,intent(in) :: eff(NRESP_POOLS)        !decomp coefs
      real*8,intent(in) :: frac_donor(NRESP_POOLS) !decomp coefs
      real*8,intent(in) :: Closs(PTRACE,NPOOLS)    !amt C lost per pool (gC/m2)
      
! i/o
      real*8,intent(inout) :: Tpool(PTRACE,NPOOLS)  !total plant and soil C,N pools
      
! output
      real*8,intent(out) :: Resp(PTRACE,NPOOLS)  !C lost to atm per pool (gC/m2)-->Cflux in casa_bgfluxes (step 3)
      
! ------------------------ local variables ------------------------
      integer resp_pool_index(2,NRESP_POOLS)   
      integer irtype
      integer donor_pool
      integer recvr_pool
      real*8 Out

! ----------------------------------------------------------
! set values of parameters/constants used in CASA Respiration
! these are in order in which respiration is called
      resp_pool_index = reshape ( 
     1              (/SLOW     ,PASSIVE,
     2              SLOW     ,SOILMIC,
     3              SURFMET  ,SURFMIC,
     4              SURFSTR  ,SURFMIC,
     5              SURFSTR  ,SLOW   ,
     6              SOILMET  ,SOILMIC,
     7              SOILSTR  ,SOILMIC,
     8              SOILSTR  ,SLOW   ,
     9              CWD      ,SURFMIC,
     a              CWD      ,SLOW   ,
     b              SURFMIC  ,SLOW   ,
     c              SOILMIC  ,PASSIVE,
     d              SOILMIC  ,SLOW   ,
     e              PASSIVE  ,SOILMIC/),
     &              (/2,NRESP_POOLS/)   )

! Loop over all respiring pools
      do irtype = 1, NRESP_POOLS
         donor_pool = resp_pool_index(1,irtype)
         recvr_pool = resp_pool_index(2,irtype)
         Out  = Closs(Carbon,donor_pool) * frac_donor(irtype)
         Tpool(Carbon,donor_pool) = Tpool(Carbon,donor_pool)
     &                                   - Out
         Tpool(Carbon,recvr_pool) = Tpool(Carbon,recvr_pool)  
     &                                   + (Out * eff(irtype))
         Resp(Carbon,donor_pool) =  Resp(Carbon,donor_pool) 
     &                                   + Out * (1.d0 - eff(irtype))
!** 01/10/01 make sure donor pool does not fall below zero
        if(Tpool(Carbon,donor_pool).le.0.0d0)then
            Tpool(Carbon,donor_pool)=0.d0
        end if        
      enddo  !irtype - over all 14 respiring pools

      return
      end subroutine casa_respire
      
!***********************************************************************      
      end module soilbgc
