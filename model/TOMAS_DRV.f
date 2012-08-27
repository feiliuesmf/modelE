#include "rundeck_opts.h"

      MODULE TOMAS_AEROSOL

C-----INCLUDE FILES--------------------------------------------------
      USE RESOLUTION, only : im,jm,lm     ! dimensions
      USE TRACER_COM, only : ntm
      IMPLICIT NONE 

C-----VARIABLE DECLARATIONS------------------------------------------
      integer ibins, icomp !ibins should be matched to nbins in TRACER_COM.f
      integer idiag ! number of diagnostic aerosol species
      parameter(ibins=12, icomp=9, idiag=2)

      integer srtso4, srtna, srth2o, srtecob, srtecil, srtocob,
     &	      srtocil, srtnh4, srtdust
      parameter (srtso4=1,
     &           srtna =2,
     &           srtecob=3,
     &           srtecil=4,
     &           srtocob=5,
     &           srtocil=6,
     &		 srtdust=7,
     &           srtnh4=8,
     &           srth2o=9)

C Nk and Mk contain the number and mass size distributions of the
C aerosol.  Units are #/grid cell or kg/grid cell, respectively.
C Gc are gas phase concentrations (kg/grid cell) of species
C corresponding to the aerosol species (e.g. H2SO4 for sulfate).
C Nkd and Mkd store values of Nk and Mk for diagnostic purposes.

      real*8 Nk(ibins), Mk(ibins,icomp), Gc(icomp-1)
      real*8 Nkd(ibins), Mkd(ibins,icomp), Gcd(icomp-1) 

C The following variables describe the grid cell in which the
C microphysics is operating.

      real*8 boxvol     !volume of grid cell (cm3)
      real*8 boxmass    !volume of grid cell (kg)
      real*8 temp       !temperature (K) of grid cell
      real*8 pres       !air pressure (Pa) of grid cell
      real*8 rh         !relative humidity (0-1)

C Physical properties of aerosol components

      real molwt(icomp)
      data molwt/96., 58.45, 200., 200., 200., 200., 100.,18.,18./

C Flag for which nucleation parameterizations to use
      integer bin_nuc, tern_nuc, ion_nuc, actv_nuc   !flags for binary and ternary nuc
      parameter(bin_nuc=1, tern_nuc=0, ion_nuc=0, actv_nuc=1) ! 1 = on

C soa_amp is the mass growth amplification factor (determined by the amount of
C soa that needs to be condensed
      real*8 soa_amp

C tau_soa is the 1st order timescale in which SOA condenses (days)
      real*8 tau_soa
      parameter(tau_soa=0.5d0)

C SOArate is the rate of SOA condensing (kg/s)
!      real*8, ALLOCATABLE, DIMENSION(:,:)  :: SOA_chem !SOA formation [kg]
      real*8  SOArate

      real*8 surf_area ! aerosol surface area [micon^2 cm^-3]
      real*8 ionrate ! ion pair formation rate [ion pairs cm^-3 s^-1]

      integer, dimension(101,101,101) :: binact10,binact02
      real*8,dimension(101,101,101) :: fraction10,fraction02
      integer, parameter :: ptype=7 ! number of microphysics process 
      real*8, ALLOCATABLE,dimension(:,:,:,:,:) :: AEROD
      real*8, ALLOCATABLE,DIMENSION(:,:,:) :: AQSO4oxid_mc,AQSO4oxid_ls  !1 for Convective and 2 for large-scale
      real*8, ALLOCATABLE,DIMENSION(:,:,:)  ::  h2so4_chem  !h2so4 formation rate from so2+oh [kg of H2SO4/sec\
      real*8, ALLOCATABLE,DIMENSION(:,:,:,:,:) :: N_subgridcg !number changed by subgrid coagulation
      real*8, ALLOCATABLE,DIMENSION(:,:,:,:,:,:)  :: M_subgridcg !mass changed by subgrid coagulation
      real*8, ALLOCATABLE,DIMENSION(:,:,:,:)  :: trm_emis !trm before emission and uses for subgrid coagulation
!      real*8, ALLOCATABLE,DIMENSION(:,:,:,:)  :: tomas_emis !trm before interactive emission and uses for subgrid coagulation
      integer ncomp
      parameter(ncomp=8)

! Mie lookup tables
      REAL*8, DIMENSION(124,101,91)      :: TOMAS_QEXT, TOMAS_QSCA,
     &     TOMAS_QABS,TOMAS_GSCA !,TOMAS_QBACK 

      ! FALSE : one Radiation call
      ! TRUE  : nmodes Radiation calls
      INTEGER                            :: TOMAS_DIAG_FC = 2 !debug 
! 2=external mixing (=icomp-2) radiation calls  |
! 1=internal mixing (but AECOB is externally mixed) (ANUM_01) radiation call
! TOMAS_DIAG_FC=2 is only available now.  
#if (defined TOMAS_12_10NM) 
      real*8, parameter :: scalesizeSO4(ibins)=(/
     &     4.3760E-02, 6.2140E-02, 3.6990E-02, 1.8270E-02,
     &     4.2720E-02, 1.1251E-01, 1.9552E-01, 2.2060E-01,
     &     1.6158E-01, 7.6810E-02, 2.8884E-02, 2.0027E-04/)

      real*8, parameter :: scalesizeCARBO30(ibins)=(/
     &     1.1291E-03,4.9302E-03,1.2714E-02,3.6431E-02,
     &     1.0846E-01,2.1994E-01,2.7402E-01,2.0750E-01,
     &     9.5304E-02,2.6504E-02,1.2925E-02,1.6069E-05/) ! use for fossil fuel (bimodal)
      real*8, parameter :: scalesizeCARBO100(ibins)=(/!0.0,0.0,0.0,
     &     1.9827E-06,3.9249E-05,5.0202E-04,4.1538E-03,
     &     2.2253E-02,7.7269E-02,1.7402E-01,2.5432E-01,
     &     2.4126E-01,1.4856E-01,7.6641E-02,9.8120E-04/) ! use for biomass burning
#elif (defined TOMAS_12_3NM) 
      real*8, parameter :: scalesizeSO4(ibins)=(/0.d0,0.d0,0.d0,
     &     4.3760E-02, 6.2140E-02, 3.6990E-02, 1.8270E-02,
     &     4.2720E-02, 1.1251E-01, 1.9552E-01, 2.2060E-01,
     &     1.6158E-01, 7.6810E-02, 2.8884E-02, 2.0027E-04/)

      real*8, parameter :: scalesizeCARBO30(ibins)=(/0.d0,0.d0,0.d0,
     &     1.1291E-03,4.9302E-03,1.2714E-02,3.6431E-02,
     &     1.0846E-01,2.1994E-01,2.7402E-01,2.0750E-01,
     &     9.5304E-02,2.6504E-02,1.2925E-02,1.6069E-05/) ! use for fossil fuel (bimodal)
      real*8, parameter :: scalesizeCARBO100(ibins)=(/0.d0,0.d0,0.d0,
     &     1.9827E-06,3.9249E-05,5.0202E-04,4.1538E-03,
     &     2.2253E-02,7.7269E-02,1.7402E-01,2.5432E-01,
     &     2.4126E-01,1.4856E-01,7.6641E-02,9.8120E-04/) ! use for biomass burning
#endif

#if (defined TOMAS_15_10NM)
      real*8, parameter :: scalesizeSO4(ibins)=(/
     &     4.3760E-02, 6.2140E-02, 3.6990E-02, 1.8270E-02,
     &     4.2720E-02, 1.1251E-01, 1.9552E-01, 2.2060E-01,
     &     1.6158E-01, 7.6810E-02, 2.8884E-02, 2.0027E-04/)

      real*8, parameter :: scalesizeCARBO30(ibins)=(/
     &     3.8100E-03,2.0700E-02,7.2900E-02,1.6750E-01,
     &     2.5000E-01,2.4500E-01,1.5510E-01,6.4200E-02,1.7320E-02,
     &     3.0320E-03,3.7079E-04,2.3265E-07/) ! use for fossil fuel

      real*8, parameter :: scalesizeCARBO100(ibins)=(/
     &     1.9827E-06,3.9249E-05,5.0202E-04,4.1538E-03,
     &     2.2253E-02,7.7269E-02,1.7402E-01,2.5432E-01,2.4126E-01,
     &     1.4856E-01,7.6641E-02,9.8120E-04/) ! use for biomass burning
#elif (defined TOMAS_15_3NM)
      real*8, parameter :: scalesizeSO4(ibins)=(/
     &     0.d0,0.d0,0.d0,
     &     4.3760E-02, 6.2140E-02, 3.6990E-02, 1.8270E-02,
     &     4.2720E-02, 1.1251E-01, 1.9552E-01, 2.2060E-01,
     &     1.6158E-01, 7.6810E-02, 2.8884E-02, 2.0027E-04/)

      real*8, parameter :: scalesizeCARBO30(ibins)=(/
     &     0.d0,0.d0,0.d0,
     &     3.8100E-03,2.0700E-02,7.2900E-02,1.6750E-01,
     &     2.5000E-01,2.4500E-01,1.5510E-01,6.4200E-02,1.7320E-02,
     &     3.0320E-03,3.7079E-04,2.3265E-07/) ! use for fossil fuel

      real*8, parameter :: scalesizeCARBO100(ibins)=(/
     &     0.d0,0.d0,0.d0,
     &     1.9827E-06,3.9249E-05,5.0202E-04,4.1538E-03,
     &     2.2253E-02,7.7269E-02,1.7402E-01,2.5432E-01,2.4126E-01,
     &     1.4856E-01,7.6641E-02,9.8120E-04/) ! use for biomass burning
#endif

!#if (defined TOMAS_30_10NM)

!#elif (defined TOMAS_30_3NM)
           
!#endif
      END MODULE TOMAS_AEROSOL


C     **************************************************
C     *  aerophys                                      *
C     **************************************************

C     WRITTEN BY Peter Adams, August 2000

C     This subroutine does aerosol microphysics (nucleation,
C     coagulation, condensation).  Gas-phase chemistry, aqueous
C     chemistry, dry and wet deposition are handled elsewhere in
C     the GCM.


      SUBROUTINE TOMAS_DRV 

C-----INCLUDE FILES--------------------------------------------------
      USE DOMAIN_DECOMP_ATM, only : GRID, getDomainBounds, 
     &   write_parallel ,am_i_root
      USE TOMAS_AEROSOL 
      USE TRACER_COM

      USE TRDIAG_COM, only : taijs=>taijs_loc,taijls=>taijls_loc
     *     ,ijts_TOMAS,itcon_TOMAS
!      USE AEROSOL_SOURCES, only: off_HNO3
      USE FLUXES, only: tr3Dsource
      USE RESOLUTION, only : im,jm,lm     ! dimensions
      USE ATM_COM, only :   t            ! potential temperature (C)
     $                     ,q            ! saturated pressure
      USE MODEL_COM, only : dtsrc
      USE GEOM, only: axyp,imaxj,BYAXYP
      USE CONSTANT,   only:  lhe,mair,gasc   
      USE ATM_COM,   only: pmid,pk,byam,gz, am   ! midpoint pressure in hPa (mb)
!                                           and pk is t mess up factor
!                                           BYAM  1/Air mass (m^2/kg)
      IMPLICIT NONE

C-----VARIABLE DECLARATIONS------------------------------------------

      INTEGER J_0, J_1, I_0, I_1

      integer i,j,l,n,jc,mt,k,np  !counters
      integer mpnum       !microphysical process id #
      real adt            !aerosol microphysics time step (seconds)
      real*8 qsat         !used in RH calculation
      integer tracnum
      integer flag
      real*8 frac
      
      real*8 Nkout(NBINS), Mkout(NBINS,icomp)
      real*8 Gcout(icomp-1)
      real*8 tot_aam  ! total aerosol ammonia per grid cell across all bins
      real*8 Gcavg ! average h2so4 concentration during timstep
             
      real*8 H2SO4rate_o, SOAmass ! H2SO4rate for the specific gridcell, SOA mass to be condensed
      integer num_iter
      real*8 fn                 ! nucleation rate of clusters cm-3 s-1
      real*8 fn1                ! formation rate of particles to first size bin cm-3 s-1
      real*8 tot_n_1, tot_n_1a, tot_n_2, tot_n_i ! used for nitrogen mass checks
      real*8 tot_s_1, tot_s_1a,tot_s_1b, tot_s_2 ! used for sulfur mass checks
      real*8 Nknuc(ibins), Mknuc(ibins, icomp)
      real*8 Nkcond(ibins),Mkcond(ibins,icomp)
      real*8 INIT_Nk(ibins),INIT_Mk(ibins,icomp)
      real*8 INIT_H2SO4,INIT_NH3,INIT_NH4,INIT_SOA
      real*8 TSUM(2)

c$$$      real*4, dimension(GRID%J_STRT_HALO:GRID%J_STOP_HALO,lm) :: 
c$$$     &     nucrate,nucrate1

C-----CODE-----------------------------------------------------------
    
C****
C**** Extract useful local domain parameters from "grid"
C****
      call getDomainBounds(grid, J_STRT=J_0,       J_STOP=J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

!debug      nucrate(J_0:J_1,I_0:I_1)        = 0.d0  !DIAG: nucleation rate diagnotics (Jnuc)
!debug      nucrate1(J_0:J_1,I_0:I_1)       = 0.d0  !DIAG: new particle formation rate at lowest boundary in TOMAS


C     Loop over all grid cells
      DO L=1,LM                            
         DO J=J_0,J_1                          
            DO I=I_0,IMAXJ(J)


               temp = pk(l,i,j)*t(i,j,l) !should be in [K]
               rh = MIN(1.,q(i,j,l)/QSAT(temp,lhe,pmid(l,i,j))) ! rH [0-1]
               pres= pmid(l,i,j)*100. ! pmid in [hPa]
               boxmass=am(l,i,j)*axyp(i,j) !kg of air
               boxvol=boxmass/mair*1000.d0
     &              *gasc*temp/pres*1e6 !cm3
 
Cjrp  Initialize all components condensible gas values to zero      
Cjrp  Gc(srtso4) will remain zero until within cond_nuc where the
Cjrp  pseudo steady state H2SO4 concentration will be put in this place.

               Gc(:)=0.d0
               
C     Swap T0M into Nk, Mk, Gc arrays

               do n=1,ibins
                  Nk(n)=TRM(i,j,l,IDTNUMD-1+n)
                  Mk(n,srtso4)=TRM(i,j,l,IDTSO4-1+n)
                  Mk(n,srtna) =TRM(i,j,l,IDTNA -1+n)
                  MK(n,srtecob)=TRM(i,j,l,IDTECOB -1+n)
                  MK(n,srtecil)=TRM(i,j,l,IDTECIL -1+n)
                  MK(n,srtocob)=TRM(i,j,l,IDTOCOB -1+n)
                  MK(n,srtocil)=TRM(i,j,l,IDTOCIL -1+n)      
                  Mk(n,srtdust)=TRM(i,j,l,IDTDUST -1+n)            
                  Mk(n,srth2o)=TRM(i,j,l,IDTH2O-1+n)
                  Mk(n,srtnh4)=0.
               enddo

               INIT_NK(:) = NK(:)
               INIT_MK(:,:)=MK(:,:)
               INIT_H2SO4 = H2SO4_chem(I,J,L)*dtsrc
               INIT_NH3=TRM(I,J,L,n_NH3)
               INIT_NH4=TRM(I,J,L,n_NH4)
               INIT_SOA=TRM(I,J,L,n_SOAgas)

! swap NH3 from giss to tomas               
               tot_n_i = TRM(i,j,l,n_NH3)*14.d0/17.d0 + 
     &              TRM(i,j,l,n_NH4)*14.d0/18.d0

               call NH3_GISStoTOMAS(TRM(i,j,l,n_NH3), !YHL - This might not be needed. 
     &              TRM(i,j,l,n_NH4),Gc,Mk)
              
                                ! nitrogen and sulfur mass checks
                                ! get the total mass of N
               tot_n_1 = Gc(srtnh4)*14.d0/17.d0
               do k=1,ibins
                  tot_n_1 = tot_n_1 + Mk(k,srtnh4)*14.d0/18.d0
               enddo

               H2SO4rate_o = H2SO4_chem(i,j,l) !kg of h2so4/sec  (from SO2+OH)
               SOAmass=TRM(i,j,l,n_SOAgas) !kg of SOA 
               
! Do water eqm at appropriate times
               call ezwatereqm(Mk)

               call storenm()
               call mnfix(Nk,Mk)
               mpnum=6 
               call aerodiag(mpnum,i,j,l)

! in-cloud oxidation is not affected by time step.  
! in-cloud oxidation! move from clouds2.f to here        

               call storenm()
                  call aqoxid(i,j,l,.TRUE.) ! Moist Convective clouds

               mpnum=4 
               call aerodiag(mpnum,i,j,l)

! Do water eqm at appropriate times
               call ezwatereqm(Mk)

               call storenm()
               call mnfix(Nk,Mk)
               mpnum=6 
               call aerodiag(mpnum,i,j,l)


               call storenm()
                  call aqoxid(i,j,l,.false.) ! Large scale clouds

               mpnum=5 
               call aerodiag(mpnum,i,j,l) 


! Do water eqm at appropriate times
               call ezwatereqm(Mk)

               call storenm()
               call mnfix(Nk,Mk)
               mpnum=6 
               call aerodiag(mpnum,i,j,l)

! timestep is shorten to be 10 mins 

               adt=dtsrc/3. 

               do mt=1,3        ! 3 * 10 min inside

                                !calculate SOA to condense
               SOArate = SOAmass*(1.d0-
     &              exp(-adt/(tau_soa*3600.*24.)))/adt

               SOAmass=SOAmass-SOArate*adt ! updated SOA mass within the time step loop

               if(H2SO4rate_o.le.0.) THEN
                  if(am_i_root())
     &                 print*,'problem in soa_amp',i,j,l,
     &                 h2so4rate_o,SOArate

                  soa_amp=0. !no SOA condensing! 
               else

               soa_amp = SOArate/H2SO4rate_o  ! TOMAS- floating invalid due to zero division??
               endif

C     ****************
C     Aerosol dynamics
C     ****************                                  

! get the total mass of S
               tot_s_1 = H2SO4rate_o*adt*32.d0/98.d0
               do k=1,ibins
                  tot_s_1 = tot_s_1 + Mk(k,srtso4)*32.d0/96.d0
               enddo

               Gcavg = 0.0

               Gc(:)=0.d0  ! back to zero
!YUNHA LEE - Gc=0 is moved here due to shorter time step.
 
! WHAT TO DO WITH Gc initialization!! ? 

Cjrp  Initialize all components condensible gas values to zero      
Cjrp  Gc(srtso4) will remain zero until within cond_nuc where the
Cjrp  pseudo steady state H2SO4 concentration will be put in this place.

!               Gc(srtso4)=h2so4rate_o*adt ! this is for condensation diagnostics 

               call storenm()

C If any Nk are zero, then set them to a small value to avoid division by zero
               call cond_nuc(Nk,Mk,Gc,Nkout,Mkout,Gcout,fn,fn1,
     &             H2SO4rate_o,adt,num_iter,Nknuc,Mknuc,Nkcond,Mkcond,l)
                                !get nucleation diagnostic
 
               Mk(:,:)=Mknuc(:,:)
               Nk(:)=Nknuc(:)

!               Gc(srtso4)=h2so4rate_o*adt !to make zero nucleation diag
            
               mpnum=3 
               call aerodiag(mpnum,i,j,l)

               Mk(:,:)=Mkcond(:,:)
               Nk(:)=Nkcond(:)

               
               Gc(srtnh4)=Gcout(srtnh4)
               Gc(srtso4)=Gcout(srtso4)

               TRM(I,J,L,n_H2SO4)=Gc(srtso4)

               mpnum=1
               call aerodiag(mpnum,i,j,l)
               Mk(:,:)=Mkout(:,:)
               Nk(:)=Nkout(:)

!               nucrate(j,l)=nucrate(j,l)+fn
!               nucrate1(j,l)=nucrate1(j,l)+fn1
               
                                ! accumulate nucleation rate diagnostics
                                ! first sum for JL
               TSUM(1)=TSUM(1)+fn*boxvol*adt ! number of particles generated per kg of air in timestep
               TSUM(2)=TSUM(2)+fn1*boxvol*adt
                                ! IJ
c$$$               T3DC(I,J,L,1)=T3DC(I,J,L,1)+fn*boxvol*adt
c$$$               T3DC(I,J,L,2)=T3DC(I,J,L,2)+fn1*boxvol*adt      
               
                                ! nitrogen and sulfur mass checks
                                ! get the total mass of N
               tot_n_1a = Gc(srtnh4)*14.d0/17.d0
               do k=1,ibins
                  tot_n_1a = tot_n_1a + Mk(k,srtnh4)*14.d0/18.d0
               enddo
               
                                ! get the total mass of S
               tot_s_1b = 0.d0
               do k=1,ibins
                  tot_s_1b = tot_s_1b + Mk(k,srtso4)*32.d0/96.d0
               enddo               
        
               
               call storenm()
               call mnfix(Nk,Mk)
               mpnum=6 
               call aerodiag(mpnum,i,j,l)


C YHL - This should be called for next multicoag.  This does not need to be called in last timestep but do it anyway. 

                                !Coagulation
               call storenm()
               call multicoag(adt)

               mpnum=2
               call aerodiag(mpnum,i,j,l)


C     Do water eqm at appropriate times
               call eznh3eqm(Gc,Mk)
               call ezwatereqm(Mk)
                  
C     ***********************
C     End of aerosol dynamics
C     ***********************
               
                                ! do nitrogen and sulfur mass checks
                                ! get the total mass of N
               tot_n_2 = Gc(srtnh4)*14.d0/17.d0
               do k=1,ibins
                  tot_n_2 = tot_n_2 + Mk(k,srtnh4)*14.d0/18.d0
               enddo
               
                                ! get the total mass of S
               tot_s_2 = 0.d0
               do k=1,ibins
                  tot_s_2 = tot_s_2 + Mk(k,srtso4)*32.d0/96.d0
               enddo
               

               if(mt.eq.1)then  !only firs time step
               if(am_i_root())then

                  if (abs(tot_n_2-tot_n_1)/tot_n_1.gt.1.0D-4)then
                     print*,'Nitrogen not conserved in aerophys'
! 1                   print*,'i',i,'j',j,'l',l
!                     print*,'Init,Init1,Intm,Final',tot_n_i,tot_n_1,
!     *                    tot_n_1a,tot_n_2
                  endif
                  
                  if (abs(tot_s_2-tot_s_1)/tot_s_1.gt.1.0D-4)then !TOMAS - increase from 1.0D-4 
                     print*,'Sulfur not conserved in aerophys'
!                     print*,'i',i,'j',j,'l',l
!                     print*,'Init,Init1,Intm,Final',tot_s_1,
!     *                    tot_s_1b,tot_s_2
                  endif
               endif
               endif

             enddo              ! timestep
              
C     Check for negative tracer problems
               flag=0
               do n=1,NBINS
                  if (Nk(n) .lt. 0.0) then
                     write(*,*) 'Nk < 0 for i,j,l,bin:',i,j,l,n
                     flag=1
                  endif
                  do jc=1,icomp
                     if (Mk(n,jc) .lt. 0.0) then
                    write(*,*) 'Mk < 0 for i,j,l,bin,comp:',i,j,l,n,jc
                        flag=1
                     endif
                  enddo
               enddo
               if (flag .eq. 1) 
     &              CALL STOP_MODEL('- tracer in TOMAS_DRV',255)   

!Save diagnostics! 
               do n=1,ibins       
!     Aerosol number             
                  tracnum=IDTNUMD-1+n 
                  tr3Dsource(i,j,l,nOther,tracnum)=
     &                 (NK(N)-INIT_NK(N))/dtsrc
                  
                  do np=1,ptype
                     if (ijts_TOMAS(np,tracnum).gt.0) 
     &                taijs(i,j,ijts_TOMAS(np,tracnum)) 
     &                    =taijs(i,j,ijts_TOMAS(np,tracnum))
     &                    +AEROD(i,j,l,tracnum,np) ! /adt
                     if (itcon_TOMAS(np,tracnum).gt.0) 
     &                    call inc_diagtcb(i,j,AEROD(i,j,l,tracnum,np),
     &                    itcon_TOMAS(np,tracnum),tracnum)
                  enddo

                  do jc=1,icomp-idiag
                     tracnum=IDTSO4-1+n+ibins*(jc-1)
                     tr3Dsource(i,j,l,nOther,tracnum)=
     &                    (MK(n,jc)-INIT_Mk(n,jc))/dtsrc

                  do np=1,ptype
                     if (ijts_TOMAS(np,tracnum).gt.0) 
     &                   taijs(i,j,ijts_TOMAS(np,tracnum)) 
     &                    =taijs(i,j,ijts_TOMAS(np,tracnum))
     &                    +AEROD(i,j,l,tracnum,np) ! /adt

                     if (itcon_TOMAS(np,tracnum).gt.0) 
     &                    call inc_diagtcb(i,j,AEROD(i,j,l,tracnum,np),
     &                    itcon_TOMAS(np,tracnum),tracnum)
                  enddo

                  enddo  
                 
                  tracnum=IDTH2O-1+n 
                  tr3Dsource(i,j,l,nOther,tracnum)=
     &                 (MK(N,SRTH2O)-INIT_MK(N,SRTH2O))/dtsrc
               enddo
               

               tr3Dsource(i,j,l,nOther,n_H2SO4) =
     *              (Gc(srtSO4)-INIT_H2SO4)/dtsrc
               
               tr3Dsource(i,j,l,nChemistry,n_NH3)=
     *              (Gc(srtNH4)-INIT_NH3)/dtsrc

                                ! aerosol ammonia
               tot_aam = 0.d0
               do n=1,NBINS
                  tot_aam = tot_aam + Mk(n,srtnh4)
               enddo
               
               tr3Dsource(i,j,l,nChemistry,n_NH4)=
     *              (tot_aam-INIT_NH4)/dtsrc

               tr3Dsource(i,j,l,nChemistry,n_SOAgas)=
     *              (SOAmass-INIT_SOA)/dtsrc !total SOA rate in 30 min
               
               AEROD(i,j,l,:,:)=0.0 !for aeroupdate

C     End of loop over grid cells
            enddo               !I loop
c$$$  do K=1,2
c$$$  TAJLS(J,L,K+9) = TAJLS(J,L,K+9) + TSUM(K)
c$$$            enddo
         enddo                  !J loop
      enddo                     !L loop
      return

      END SUBROUTINE TOMAS_DRV                     !of main
      

      logical function is_nan(value)
      real*8 value
      if (abs(value).ge.0) then
         is_nan=.false.
      else
         is_nan=.true.
         endif
      return
      end function is_nan

      subroutine nanstop(value, line, var1, var2)
      real*8 value
      integer line, var1, var2

      if (abs(value).ge.0) then
      else
      write (*,*) 'line',line, var1, var2, value
      call stop_model('nan in nanstop',255)
      endif
      return
      end subroutine nanstop

C     **************************************************
C     *  drydep_getdp                                  *
C     **************************************************

C     WRITTEN BY Yunha Lee

C     This function calculates the average diameter of aerosol
C     particles in a given GCM grid cell and size bin.

      subroutine dep_getdp(i,j,l,getdp,size_density)                                            
C-----INCLUDE FILES-----------------------------------------------------
      USE TRACER_COM, only : nbins,IDTSO4,IDTNA,IDTECIL,
     &     IDTECOB,IDTOCIL,IDTOCOB,IDTDUST,IDTH2O,
     &     IDTNUMD,ntm,xk,trm
      USE CONSTANT,   only : pi,lhe,mair,gasc  
      USE ATM_COM, only :   t            ! potential temperature (C)
     $     ,q                   ! saturated pressure
     $     ,pmid                ! midpoint pressure in hPa (mb)
     $     ,pk                  ! pk is t mess up factor
      USE TOMAS_AEROSOL

      IMPLICIT NONE

C-----VARIABLE DECLARATIONS------------------------------------------

      integer i,j,l  !coordinate of GCM grid cell
      integer n,k      !tracer index and size bin index
      real density                 !density (kg/m3) of current size bin
      real mso4, mh2o, mno3, mnh4  !mass of each component (kg/grid box)
      real mecil,mecob,mocil,mocob
      real mdust,mtot,mnacl             
      real*8 mp          !particle mass (kg)
      real*8 mu          !air viscosity (kg/m s)
      real*8 qsat
      real aerodens
      real,intent(out),DIMENSION(nbins) :: getdp,size_density
      real*8 Neps  !a small number of particles (#/box)
      parameter (Neps=1.d-20)

C-----CODE-----------------------------------------------------------

!     Compute particle diameter for each bins - YUNHA LEE 
      do k=1,nbins

C     Swap GCM variables into aerosol algorithm variables
        Nk(k)=trm(i,j,l,IDTNUMD-1+k)
        Mk(k,srtso4)=trm(i,j,l,IDTSO4-1+k)
        Mk(k,srtna )=trm(i,j,l,IDTNA -1+k)
        Mk(k,srtnh4)=0.1875*Mk(k,srtso4) ! artificial for now.. 0.0!t0m(i,j,l,IDTNH4-1+n)
        MK(k,srtecob)=trm(i,j,l,IDTECOB -1+k)
        MK(k,srtecil)=trm(i,j,l,IDTECIL -1+k)
        MK(k,srtocob)=trm(i,j,l,IDTOCOB -1+k)
        MK(k,srtocil)=trm(i,j,l,IDTOCIL -1+k) 
        MK(k,srtdust)=trm(i,j,l,IDTDUST -1+k) 
        Mk(k,srth2o)= trm(i,j,l,IDTH2O-1+k) !I don't think this is necessary!
      enddo

      temp = pk(l,i,j)*t(i,j,l) !should be in [K]
      rh = MIN(1.,q(i,j,l)/QSAT(temp,lhe,pmid(l,i,j))) ! rH [0-1]
            
      call mnfix(Nk,Mk)  
      call ezwatereqm(Mk)

      do k=1,nbins 
         if (Nk(k) .eq. 0.0) then
            if (Mk(k,srtso4) .gt. 1.) then
               print*, 'ERROR in getdp - # = but mass > 0',i,j
               print*, 'bin=',k
               print*, 'TRM(#)=',Nk(k)
               print*, 'TRM(SO4)=',Mk(k,srtso4)
               print*, 'TRM(NACL)=',Mk(k,srtna)
               print*, 'TRM(OCIL)=',Mk(k,srtocil)
               call stop_model('ERROR IN getdp',255)
            endif
         endif

         mso4=Mk(k,srtso4) 
         mnacl=Mk(k,srtna)
         mno3=0.e0
         if ((mso4+mno3) .lt. 1.e-8) mso4=1.e-8
         mnh4=0.1875*mso4  !assume ammonium bisulfate
         mecob=Mk(k,srtecob)
         mecil=Mk(k,srtecil)
         mocil=Mk(k,srtocil)
         mocob=Mk(k,srtocob)
         mdust=Mk(k,srtdust)          
         mh2o=Mk(k,srth2o)   

!in CLOUDS2.f - some tracers goes negative..
!So, set to zero to prevent a problem in density calculation
         if(mnacl.lt.0) mnacl=0.
         if(mecob.lt.0) mecob=0.
         if(mecil.lt.0) mecil=0.
         if(mocob.lt.0) mocob=0.
         if(mocil.lt.0) mocil=0.
         if(mdust.lt.0) mdust=0.
         if(mh2o.lt.0) mh2o=0.

         density=aerodens(mso4,mno3,mnh4 !mno3 taken off!
     *        ,mnacl,mecil,mecob,mocil,mocob,mdust,mh2o) !assume bisulfate   

         mtot= 1.1875*Mk(k,srtso4)+mnacl+mecil+mecob+
     *        mocil+mocob+mdust+mh2o
   
         size_density(k)=density      

         if (Nk(k) .gt. Neps.and.mtot.gt.0.) then
            mp=mtot/Nk(k)
         else
            mp=sqrt(xk(k+1)*xk(k))
!            if(Nk(k) .gt. Neps) 
!     &           print*,'Warning in getdp:#>Neps but mtot=0',
!     &           k,mtot,Nk(k)
         endif
         
!     fix unrealistically large mp for low aerosol conc.
         if (mp .gt. 1.d3*xk(NBINS+1)) then
            
            if ((Nk(k) .lt. 1.d5) .and. !negligible amount of aerosol - fudge mp
     &           (Mk(k,srtso4) .lt. 3.)) then
               mp=sqrt(xk(k+1)*xk(k))
            else
               if (Nk(k) .gt. 1.d12) then
!MODELE-TOMAS: during CONDSE, TM(H2O) is so large that causes too big mp. 
!MODELE-TOMAS: So, if dry mass is less than the max size boundary, just take the max mp. 
                  if((mtot-Mk(k,srth2o)).lt.
     &                 1.d1*xk(nbins+1)*Nk(k))then
                     print*,'Fudge mp in getdp: large mp by AH2O'
                     mp=sqrt(xk(nbins+1)*xk(nbins))                     
                  else
                  print*,'ERROR in getdp: mp too large'
                  print*, 'bin=',k
                  print*, 'TM(#)=', Nk(k)
                  print*, 'TM(SO4)=', mso4, mh2o
                  print*, 'TM(NACL)=', mnacl, mdust
                  print*, 'TM(OC)=',mocob,mocil
                  print*, 'TM(EC)=',mecob,mecil
                  call stop_model('mp too large getdp',255)
                  endif
               endif
            endif
         endif
         getdp(k)=(6.d0*mp/(pi*size_density(k)))**(1.d0/3.d0)
      enddo

      RETURN
      END subroutine dep_getdp

      subroutine readfraction(infile,fraction2)
      
      IMPLICIT NONE
C-----VARIABLE DECLARATIONS---------------------------------------------
      
      character*17 infile
      integer innum, ii, jj, kk
      
      real*8,intent(out),dimension(101,101,101):: fraction2
      parameter (innum=580)
 1    format(f6.5)
      open(unit=innum,file=infile,FORM='FORMATTED',STATUS='OLD')
      do ii=1,101
        do jj=1,101
          do kk=1,101
            read(innum,1) fraction2(kk,jj,ii)
            if (fraction2(kk,jj,ii).gt.1.) fraction2(kk,jj,ii)=0.
          enddo
        enddo
      enddo
!     print*,'fraction last',fraction(101,101,101)
      close(innum)
      return 
      end subroutine readfraction
			 
      subroutine readbinact(infile,binact)
C-----INCLUDE FILES-----------------------------------------------------
      USE TOMAS_AEROSOL, ONLY : ibins      
      IMPLICIT NONE
C-----VARIABLE DECLARATIONS---------------------------------------------
      
      character*15 infile
      integer innum, ii, jj, kk
      integer,intent(out),dimension(101,101,101):: binact
      parameter (innum=590)
 1    format(I2)
      open(unit=innum,file=infile,FORM='FORMATTED',STATUS='OLD')
      do ii=1,101
        do jj=1,101
          do kk=1,101
            read(innum,1) binact(kk,jj,ii)
            if (binact(kk,jj,ii).eq.0) binact(kk,jj,ii)=ibins+1
          enddo
        enddo
      enddo
      close(innum)
      return
      end subroutine readbinact
      
      subroutine getfraction(tr_conv,tm,fract)
      USE TOMAS_AEROSOL, ONLY : binact02,binact10,
     &     fraction02,fraction10 
      USE TRACER_COM, only : nbins,ntm,IDTECIL,
     &     IDTOCIL,IDTOCOB,IDTSO4,IDTNA,IDTDUST,
     &     IDTECOB

      IMPLICIT NONE
C-----VARIABLE DECLARATIONS---------------------------------------------
      
      real mecil, mocil, mocob, mso4, mnacl,mdust, mtot
      real xocil, xso4, xnacl
      integer iso4, inacl, iocil,k
      integer getbinact
      real*8,dimension(nbins) ::  fract
      REAL*8, DIMENSION(ntm) :: TM
      LOGICAL TR_CONV
      
      do k=1, nbins
        mecil=TM(IDTECIL-1+k)
	mocil=TM(IDTOCIL-1+k)
	mocob=TM(IDTOCOB-1+k)
	mso4=TM(IDTSO4-1+k)*1.2 !account for ammonium sulfate
	mnacl=TM(IDTNA-1+k)
	mdust=TM(IDTDUST-1+k)
	mtot=mecil+mocil+mocob+mso4+mnacl+mdust+1.e-20
	xocil=mocil/mtot
	xso4=mso4/mtot
	xnacl=mnacl/mtot
	iso4=min(101,int(xso4*100)+1)
	inacl=min(101,int(xnacl*100)+1)
	iocil=min(101,int(xocil*100)+1)
        
        if(xso4.lt.0.or.xnacl.lt.0.or. xocil.lt.0)then
          print*,'wrong getfraction'
          print*,'mass',mso4,mnacl,mecil,mocob
     &         ,mocil,mdust
          call stop_model('wrong getfraction',255)
        endif
        
        if (tr_conv)then        !convective clouds
          getbinact=binact10(iso4,inacl,iocil)
          
          if(binact10(iso4,inacl,iocil).lt.0)then
            print*,'wrong binact',binact10(iso4,inacl,iocil)
            print*,'iso4',iso4,inacl,iocil
            print*,'mass',mso4,mnacl,mocil,mtot
            call stop_model('wrong binact',255)
          endif
          
          if (getbinact.gt.k) then
            fract(k)=0.         !not activated
          else if (getbinact.eq.k) then           
            fract(k)=fraction10(iso4,inacl,iocil) !partly activated
          else
            fract(k)=1.         !all sizebin activated
          endif
          if(getbinact.le.2)
     &         print*,'CONV CLD',getbinact,k,iso4,inacl,iocil
        else                    !large-scale
          getbinact=binact02(iso4,inacl,iocil)
          
          if(binact02(iso4,inacl,iocil).lt.0)then
            print*,'wrong binact',binact02(iso4,inacl,iocil)
            print*,'iso4',iso4,inacl,iocil
            print*,'mass',mso4,mnacl,mocil,mtot
          endif
          
          if (getbinact.gt.k) then
            fract(k)=0.         !not activated
          else if (getbinact.eq.k) then            
            fract(k)=fraction02(iso4,inacl,iocil) !partly activated
          else
            fract(k)=1.         !all sizebin activated            
            if(getbinact.le.4)
     &           print*,'STRAT CLD',getbinact,k,iso4,inacl,iocil
            
          endif
        endif
      enddo                     !k
      return
      end subroutine getfraction
      
C     **************************************************
C     *  aqoxid                                        *
C     **************************************************

C     WRITTEN BY Peter Adams, June 2000

C     This routine takes an amount of SO4 produced via in-cloud
C     oxidation and condenses it onto an existing aerosol size
C     distribution.  It assumes that only particles larger than the
C     critical activation diameter activate and that all of these have
C     grown to roughly the same size.  Therefore, the mass of SO4 
C     produced by oxidation is partitioned to the various size bins
C     according to the number of particles in that size bin.  
C     Values of tau are calculated for each size bin accordingly and
C     the cond subroutine is called to update Nk and Mk.

C-----INPUTS------------------------------------------------------------

C-----OUTPUTS-----------------------------------------------------------

      SUBROUTINE aqoxid(i,j,l,tr_conv)

C-----INCLUDE FILES-----------------------------------------------------

      USE TOMAS_AEROSOL
      USE TRACER_COM, only : ntm, IDTECIL,
     &       IDTOCIL,IDTOCOB,IDTSO4,IDTNA,IDTDUST,
     &       IDTECOB,IDTH2O,xk,nbins

      IMPLICIT NONE
C-----VARIABLE DECLARATIONS---------------------------------------------
      real*8 moxid !mass of new sulfate from in-cloud oxid.
      real*8, dimension(ibins) :: fraction    !fraction activated for every sizebin
      real*8 Nact, Mact  !#/mass of activated particles
      real*8 mpo   !initial particle mass (kg)
      real*8 mpw   !initial particle wet mass (kg)
      real*8 aqtau(ibins)
      integer k,mpnum,n,tracnum,i,j,l
      real*8 Nko(ibins), Mko(ibins, icomp) !input to cond routine
      real*8 Nkf(ibins), Mkf(ibins, icomp) !output from cond routine
      real*8 tdt      !the value 2/3
      real*8,parameter :: eps=1.d-40
      integer jc
      real*8 frac      
      real*8 WR                ! wet ratio = total mass/ dry mass (win, 5/15/06)
      real*8 mox(ibins) !mass of new sulfate per particle in each bin
      real*8 tot_aam ! total aerosol ammonia
      real*8 TM(ntm) ! total aerosol ammonia

      LOGICAL TR_CONV

C-----CODE--------------------------------------------------------------

              
      if (tr_conv) then
         moxid=AQSO4oxid_mc(i,j,l) 
      else
         moxid=AQSO4oxid_ls(i,j,l) 
      endif

      if (moxid.eq.0.d0) return

      tdt=2.d0/3.d0

      TM(:)=0.0
!only mass needed for getfraction
      do n=1,IBINS
         do jc=1,icomp-idiag
            tracnum=IDTSO4-1+n+ibins*(jc-1)
            TM(tracnum)=Mk(n,jc)               
         enddo
         tracnum=IDTH2O-1+n
         TM(tracnum)=Mk(n,srth2o)
      enddo

      if (tr_conv) then
         CALL getfraction (.true.,TM,FRACTION) !1% supersaturation assumption
      else
         CALL getfraction(.false.,TM,FRACTION) !0.2% supersaturation assumption
      endif

      Nact=0.0
      Mact=0.0
      do k=1,ibins
         Nact=Nact+Nk(k)*fraction(k)
         do jc=1,icomp-idiag
            Mact=Mact+Mk(k,jc)*fraction(k)
         enddo
      enddo

      if ((Mact+moxid)/(Nact+eps) .gt. xk(ibins)) then !YHL- I change xk(ibins-1) to xk(ibins)
!            if (TAU .gt. 8350.) then
cdebug               write(*,*) 'ERROR in aqoxidcc: Ave size out of bounds'
c$$$               write(*,*) 'Nact: ',Nact
c$$$               write(*,*) 'moxid/Mact: ',moxid,Mact
c$$$               do k=1,ibins
c$$$                  write(*,*) 'k, N, MSO4, MH2O: ',k,Nk(k),
c$$$     &                  Mk(k,srtso4),Mk(k,srth2o)
c$$$               enddo
               goto 20
c$$$            else
c$$$               !don't worry about the first two weeks
c$$$               goto 20
c$$$            endif
       endif

C Calculate tau for each size bin
      moxid=moxid/(Nact+eps) !now kg H2SO4 per activated particle
      do k=1,ibins
         mox(k)=fraction(k)*moxid
         if (fraction(k) .eq. 0.) then
            !too small to activate - no sulfate for this particle
            aqtau(k)=0.0
         else
            !activated particle - calculate appropriate tau
            mpo=0.0
            mpw=0.0
            !WIN'S CODE MODIFICATION 6/19/06
            !THIS MUST CHANGED WITH THE NEW dmdt_int.f
            do jc=1,icomp-idiag
               mpo = mpo+Mk(k,jc)  !accumulate dry mass
            enddo
            do jc=1,icomp
               mpw = mpw+Mk(k,jc)  ! have wet mass include amso4
            enddo
            WR = mpw/mpo  !WR = wet ratio = total mass/dry mass
            if (Nk(k) .gt. 0.d0) then
               mpw=mpw/Nk(k)
               aqtau(k)=1.5d0*((mpw+mox(k)*WR)**tdt-mpw**tdt)  !added WR to moxid term (win, 5/15/06)
            else
               !nothing in this bin - set tau to zero
               aqtau(k)=0.0
               mox(k)=0.d0
            endif
         endif
      enddo

      Nko(:)=Nk(:)
      Mko(:,:)=Mk(:,:)

      call tmcond(aqtau,xk,Mko,Nko,Mkf,Nkf,srtso4,mox)

      Nk(:)=Nkf(:)
      Mk(:,:)=Mkf(:,:)   


 20      continue   !go here if process is skipped  

      RETURN
      END SUBROUTINE aqoxid


C     **************************************************
C     *  aerodens                                      *
C     **************************************************

C     WRITTEN BY Peter Adams, May 1999
C     November, 2001 - extended to include NaCl and bug fixed
C                      the bug was that species densities (dan, ds0,
C                      etc...) are supposed to be calculated based on
C                      *total* solute concentration, not each species
C                      contribution as it had been.
Ckpc  Jan.,2002 - extended to include carbonaceous aerosols

C     This function calculates the density (kg/m3) of a sulfate-
C     nitrate-ammonium-nacl-water mixture that is assumed to be internally
C     mixed.  

C-----INPUTS------------------------------------------------------------

C     mso4, mno3, mnh4, mh2o, mnacl - These are the masses of each aerosol
C     component.  Since the density is an intensive property,
C     these may be input in a variety of units (ug/m3, mass/cell, etc.).

C-----OUTPUTS-----------------------------------------------------------

      real FUNCTION aerodens(mso4,mno3,mnh4,mnacl,mecil,
     & mecob,mocil,mocob,mdust,mh2o)

      IMPLICIT NONE

C-----ARGUMENT DECLARATIONS---------------------------------------------

      real mso4, mno3, mnh4, mnacl, mecil,mecob,mocil,mocob,mdust,mh2o

C-----VARIABLE DECLARATIONS---------------------------------------------

      real inodens, idensity, dec,doc,ddust
!      external inodens

C     In the lines above, "an" refers to ammonium nitrate, "s0" to 
C     sulfuric acid, "s1" to ammonium bisulfate, and "s2" to ammonium sulfate.
C     "nacl" or "ss" is sea salt.

C-----ADJUSTABLE PARAMETERS---------------------------------------------
     	parameter(dec=2200., doc=1400., ddust=2650.)
C-----CODE--------------------------------------------------------------=

      idensity=inodens(mso4, mno3,mnh4, mnacl, mh2o)

      aerodens=(idensity*(mso4+mno3+mnh4+mnacl+mh2o) !mno3 taken out! 
     &  +dec*(mecil+mecob)+doc*(mocil+mocob)+ddust*mdust)
     &  /(mso4+mno3+mnh4+mnacl+mh2o+mecil+mecob+mocil+mdust+mocob)
      RETURN
      END FUNCTION aerodens



C     **************************************************
C     *  inodens                                      *
C     **************************************************

C     WRITTEN BY Peter Adams, May 1999
C     November, 2001 - extended to include NaCl and bug fixed
C                      the bug was that species densities (dan, ds0,
C                      etc...) are supposed to be calculated based on
C                      *total* solute concentration, not each species
C                      contribution as it had been.

C     This function calculates the density (kg/m3) of a sulfate-
C     nitrate-ammonium-nacl-water mixture that is assumed to be internally
C     mixed.  

C-----Literature cited--------------------------------------------------
C     I. N. Tang and H. R. Munkelwitz, Water activities, densities, and
C       refractive indices of aqueous sulfates and sodium nitrate droplets
C       of atmospheric importance, JGR, 99, 18,801-18,808, 1994
C     Ignatius N. Tang, Chemical and size effects of hygroscopic aerosols
C       on light scattering coefficients, JGR, 101, 19,245-19,250, 1996
C     Ignatius N. Tang, Thermodynamic and optical properties of mixed-salt
C       aerosols of atmospheric importance, JGR, 102, 1883-1893, 1997

C-----INPUTS------------------------------------------------------------

C     mso4, mno3, mnh4, mh2o, mnacl - These are the masses of each aerosol
C     component.  Since the density is an intensive property,
C     these may be input in a variety of units (ug/m3, mass/cell, etc.).

C-----OUTPUTS-----------------------------------------------------------

      real FUNCTION inodens(mso4, mno3, mnh4, mnacl, mh2o)

      IMPLICIT NONE

C-----ARGUMENT DECLARATIONS---------------------------------------------

      real mso4, mno3, mnh4, mnacl, mh2o

C-----VARIABLE DECLARATIONS---------------------------------------------

      real so4temp, no3temp, nh4temp, nacltemp, h2otemp  !store initial values
      real mwso4, mwno3, mwnh4, mwnacl, mwh2o            !molecular weights
      real ntot, mtot, drytot                      !total number of moles, mass
      real nso4, nno3, nnh4, nnacl, nh2o       !moles of each species
      real xso4, xno3, xnh4, xnacl, xh2o       !mole fractions
      real rso4, rno3, rnh4, rnacl, rh2o       !partial molar refractions
      real ran, rs0, rs1, rs15, rs2       !same, but for solute species
      real asr                            !ammonium/sulfate molar ratio
      real nan, ns0, ns1, ns15, ns2, nss  !moles of dry solutes (nss = sea salt)
      real xan, xs0, xs1, xs15, xs2, xss  !mass % of dry solutes - Tang (1997) eq. 10
      real dan, ds0, ds1, ds15, ds2, dss  !binary solution densities - Tang (1997) eq. 10
      real mwan, mws0, mws1, mws15, mws2  !molecular weights
      real yan, ys0, ys1, ys15, ys2, yss  !mole fractions of dry solutes
      real yh2o
      real d                              !mixture density
      real xtot

C     In the lines above, "an" refers to ammonium nitrate, "s0" to 
C     sulfuric acid, "s1" to ammonium bisulfate, and "s2" to ammonium sulfate.
C     "nacl" or "ss" is sea salt.

C-----ADJUSTABLE PARAMETERS---------------------------------------------
      parameter(mwso4=96., mwno3=62., mwnh4=18., mwh2o=18., 
     &          mwnacl=58.45)
      parameter(mwan=mwnh4+mwno3, mws0=mwso4+2., mws1=mwso4+1.+mwnh4,
     &          mws2=2*mwnh4+mwso4)
C-----CODE--------------------------------------------------------------

C Save initial component masses to restore later 
!      mno3=0.d0
      so4temp=mso4
      no3temp=mno3
      nh4temp=mnh4
      h2otemp=mh2o
      nacltemp=mnacl

C Calculate mole fractions
      mtot = mso4+mno3+mnh4+mnacl+mh2o
      drytot = mso4+mno3+mnh4+mnacl
      if (drytot .lt. 1.e-15) then
      inodens=1000.
      return
      endif
      nso4 = mso4/mwso4
      nno3 = mno3/mwno3  !nno3 =zero 
      nnh4 = mnh4/mwnh4
      nnacl = mnacl/mwnacl
      nh2o = mh2o/mwh2o
      ntot = nso4+nno3+nnh4+nnacl+nh2o
      xso4 = nso4/ntot
      xno3 = nno3/ntot
      xnh4 = nnh4/ntot
      xnacl = nnacl/ntot
      xh2o = nh2o/ntot
!      call nanstop(mtot,92,0,0)
C If there are more moles of nitrate than ammonium, treat unneutralized
C HNO3 as H2SO4
      if (nno3 .gt. nnh4) then  !will never occur as no3 is always zero
         !make the switch
         nso4=nso4+(nno3-nnh4)
         nno3=nnh4
         mso4=nso4*mwso4
         mno3=nno3*mwno3

         !recalculate quantities
         mtot = mso4+mno3+mnh4+mnacl+mh2o
         nso4 = mso4/mwso4
         nno3 = mno3/mwno3
         nnh4 = mnh4/mwnh4
         nnacl = mnacl/mwnacl
         nh2o = mh2o/mwh2o
         ntot = nso4+nno3+nnh4+nnacl+nh2o
         xso4 = nso4/ntot
         xno3 = nno3/ntot
         xnh4 = nnh4/ntot
         xnacl = nnacl/ntot
         xh2o = nh2o/ntot

      endif

C Calculate the mixture density
C Assume that nitrate exists as ammonium nitrate and that other ammonium
C contributes to neutralizing sulfate
      nan=nno3
      if (nnh4 .gt. nno3) then 
         !extra ammonium
         asr=(nnh4-nno3)/nso4 
      else
         !less ammonium than nitrate - all sulfate is sulfuric acid
         asr=0.0                !if nnh4=0, then asr=0 
      endif
      if (asr .ge. 2.) asr=2.0
      if (asr .ge. 1.) then
         !assume NH4HSO4 and (NH4)2(SO4) mixture
         !NH4HSO4
         ns1=nso4*(2.-asr)
         !(NH4)2SO4
         ns2=nso4*(asr-1.)
         ns0=0.0
      else
         !assume H2SO4 and NH4HSO4 mixture
         !NH4HSO4
         ns1=nso4*asr
         !H2SO4
         ns0=nso4*(1.-asr)
         ns2=0.0
      endif

      !Calculate weight percent of solutes
      xan=nan*mwan/mtot*100.
      xs0=ns0*mws0/mtot*100.
      xs1=ns1*mws1/mtot*100.
      xs2=ns2*mws2/mtot*100.
      xnacl=nnacl*mwnacl/mtot*100.
      xtot=xan+xs0+xs1+xs2+xnacl
      !Calculate binary mixture densities (Tang, eqn 9)
      dan=0.9971 +4.05e-3*xtot +9.0e-6*xtot**2.
      ds0=0.9971 +7.367e-3*xtot -4.934e-5*xtot**2. +1.754e-6*xtot**3.
     &       -1.104e-8*xtot**4.
      ds1=0.9971 +5.87e-3*xtot -1.89e-6*xtot**2. +1.763e-7*xtot**3.
      ds2=0.9971 +5.92e-3*xtot -5.036e-6*xtot**2. +1.024e-8*xtot**3.
      dss=0.9971 +7.41e-3*xtot -3.741e-5*xtot**2. +2.252e-6*xtot**3.
     &       -2.06e-8*xtot**4.

      !Convert x's (weight percent of solutes) to fraction of dry solute (scale to 1)
      xtot=xan+xs0+xs1+xs2+xnacl
      xan=xan/xtot
      xs0=xs0/xtot
      xs1=xs1/xtot
      xs2=xs2/xtot
      xnacl=xnacl/xtot

      !Calculate mixture density
      d=1./(xan/dan+xs0/ds0+xs1/ds1+xs2/ds2+xnacl/dss)  !Tang, eq. 10
      if (abs(d).ge.0) then
      else
         write(*,*) d,xtot,xan,xs0,xs1,xs2,xnacl,dan,ds0,ds1,ds2,dss
         write(*,*) 'woo',asr,mtot,mso4,mno3,mnh4,mnacl,mh2o
         write(*,*) 'woo',so4temp,no3temp,nh4temp,nacltemp,h2otemp
         call stop_model('nan in inodens',255)
      endif
      if ((d .gt. 2.) .or. (d .lt. 0.997)) then
         write(*,*) 'ERROR in inodens'
         write(*,*) mso4,mno3,mnh4,mnacl,mh2o
         call stop_model('ERROR in inodens',255)
      endif

C Restore masses passed
      mso4=so4temp
      mno3=no3temp
      mnh4=nh4temp
      mnacl=nacltemp
      mh2o=h2otemp

C Return the density
      inodens=1000.*d    !Convert g/cm3 to kg/m3

      RETURN
      END FUNCTION inodens




C     **************************************************
C     *  ezwatereqm                                    *
C     **************************************************

C     WRITTEN BY Peter Adams, March 2000

C     This routine uses the current RH to calculate how much water is 
C     in equilibrium with the aerosol.  Aerosol water concentrations 
C     are assumed to be in equilibrium at all times and the array of 
C     concentrations is updated accordingly.

C-----INPUTS------------------------------------------------------------

C-----OUTPUTS-----------------------------------------------------------

      SUBROUTINE ezwatereqm(Mke)

      USE TOMAS_AEROSOL

      IMPLICIT NONE

C-----VARIABLE DECLARATIONS---------------------------------------------
      integer k
      real*8 Mke(ibins,icomp)
      real*8 so4mass, naclmass, ocilmass
      real*8 wrso4, wrnacl, wrocil
      real*8 rhe
      real*8 waterso4, waternacl, waterocil

C     VARIABLE COMMENTS...

C     This version of the routine works for sulfate and sea salt
C     particles.  They are assumed to be externally mixed and their
C     associated water is added up to get total aerosol water.
C     wr is the ratio of wet mass to dry mass of a particle.  Instead
C     of calling a thermodynamic equilibrium code, this routine uses a
C     simple curve fits to estimate wr based on the current humidity.
C     The curve fit is based on ISORROPIA results for ammonium bisulfate
C     at 273 K and sea salt at 273 K.


C-----CODE--------------------------------------------------------------

      rhe=100.d0*rh
      if (rhe .gt. 99.d0) rhe=99.d0
      if (rhe .lt. 1.d0) rhe=1.d0

      do k=1,ibins

         so4mass=Mke(k,srtso4)*1.1875  !1.2 converts kg so4 to kg nh4hso4
         naclmass=Mke(k,srtna)      !already as kg nacl - no conv necessary
         ocilmass=MKe(k,srtocil)    !already as kg ocil

         wrso4=waterso4(rhe)
         wrnacl=waternacl(rhe)
         wrocil=waterocil(rhe)

         Mke(k,srth2o)=so4mass*(wrso4-1.d0)+naclmass*(wrnacl-1.d0)
     &                 +ocilmass*(wrocil-1.d0)

      enddo

      RETURN
      END SUBROUTINE ezwatereqm


C     **************************************************
C     *  eznh3eqm                                      *
C     **************************************************

C     WRITTEN BY Jeff Pierce, April 2007

C     This subroutine puts ammonia to the particle phase until 
C     there is 2 moles of ammonium per mole of sulfate and the remainder
C     of ammonia is left in the gas phase.

C-----INPUTS------------------------------------------------------------

C-----OUTPUTS-----------------------------------------------------------

      SUBROUTINE eznh3eqm(Gce,Mke)

      USE TOMAS_AEROSOL
      IMPLICIT NONE

C-----VARIABLE DECLARATIONS---------------------------------------------

      integer k
      real*8 tot_nh3  !total kmoles of ammonia
      real*8 tot_so4  !total kmoles of so4
      real*8 sfrac    !fraction of sulfate that is in that bin
      real*8 Gce(icomp)
      real*8 Mke(ibins,icomp)

C-----ADJUSTABLE PARAMETERS---------------------------------------------

C-----CODE--------------------------------------------------------------

      ! get the total number of kmol nh3
      tot_nh3 = Gce(srtnh4)/17.d0
      do k=1,ibins
         tot_nh3 = tot_nh3 + Mke(k,srtnh4)/18.d0
      enddo

      ! get the total number of kmol so4
      tot_so4 = 0.d0
      do k=1,ibins
         tot_so4 = tot_so4 + Mke(k,srtso4)/96.d0
      enddo

      ! see if there is free ammonia
      if (tot_nh3/2.d0.lt.tot_so4)then  ! no free ammonia
         Gce(srtnh4) = 0.d0 ! no gas phase ammonia
         do k=1,ibins
            sfrac = Mke(k,srtso4)/96.d0/tot_so4
            Mke(k,srtnh4) = sfrac*tot_nh3*18.d0 ! put the ammonia where the sulfate is
         enddo
      else ! free ammonia
         do k=1,ibins
            Mke(k,srtnh4) = Mke(k,srtso4)/96.d0*2.d0*18.d0 ! fill the particle phase
         enddo
         Gce(srtnh4) = (tot_nh3 - tot_so4*2.d0)*17.d0 ! put whats left over in the gas phase
      endif

      RETURN
      END 



C     **************************************************
C     *  waterso4                                      *
C     **************************************************

C     Adaptation of ezwatereqm used in size-resolved sulfate only sim
C     November, 2001
C     ezwatereqm WRITTEN BY Peter Adams, March 2000

C     This function uses the current RH to calculate how much water is 
C     in equilibrium with the sulfate.  Aerosol water concentrations 
C     are assumed to be in equilibrium at all times and the array of 
C     concentrations is updated accordingly.

      real*8 FUNCTION waterso4(rhe)

      IMPLICIT NONE

C-----VARIABLE DECLARATIONS---------------------------------------------
      real*8 rhe   !relative humidity (0-100 scale)


C     VARIABLE COMMENTS...

C     waterso4 is the ratio of wet mass to dry mass of a particle.  Instead
C     of calling a thermodynamic equilibrium code, this routine uses a
C     simple curve fit to estimate wr based on the current humidity.
C     The curve fit is based on ISORROPIA results for ammonium bisulfate
C     at 273 K.

C-----CODE--------------------------------------------------------------

      if (rhe .gt. 99.) rhe=99.
      if (rhe .lt. 1.) rhe=1.

         if (rhe .gt. 96.) then
            waterso4=
     &      0.7540688*rhe**3-218.5647*rhe**2+21118.19*rhe-6.801999e5
         else
         if (rhe .gt. 91.) then
            waterso4=8.517e-2*rhe**2 -15.388*rhe +698.25
         else
         if (rhe .gt. 81.) then
            waterso4=8.2696e-3*rhe**2 -1.3076*rhe +53.697
         else
         if (rhe .gt. 61.) then
            waterso4=9.3562e-4*rhe**2 -0.10427*rhe +4.3155
         else
         if (rhe .gt. 41.) then
            waterso4=1.9149e-4*rhe**2 -8.8619e-3*rhe +1.2535
         else
            waterso4=5.1337e-5*rhe**2 +2.6266e-3*rhe +1.0149
         endif
         endif
         endif
         endif
         endif

         !check for error
         if (waterso4 .gt. 30.) then
            write(*,*) 'ERROR in waterso4'
            write(*,*) rhe,waterso4
            call stop_model('ERROR in waterso4',255)
         endif

      RETURN
      END  FUNCTION waterso4


C     **************************************************
C     *  waternacl                                     *
C     **************************************************

C     WRITTEN BY Peter Adams, November 2001

C     This function uses the current RH to calculate how much water is 
C     in equilibrium with the seasalt.  Aerosol water concentrations 
C     are assumed to be in equilibrium at all times and the array of 
C     concentrations is updated accordingly.

      real*8 FUNCTION waternacl(rhe)

      IMPLICIT NONE

C-----VARIABLE DECLARATIONS---------------------------------------------
      real*8 rhe   !relative humidity (0-100 scale)

C     VARIABLE COMMENTS...

C     waternacl is the ratio of wet mass to dry mass of a particle.  Instead
C     of calling a thermodynamic equilibrium code, this routine uses a
C     simple curve fit to estimate waternacl based on the current humidity.
C     The curve fit is based on ISORROPIA results for sodium sulfate
C     at 273 K.

C-----CODE--------------------------------------------------------------

      if (rhe .gt. 99.) rhe=99.
      if (rhe .lt. 1.) rhe=1.

         if (rhe .gt. 90.) then
            waternacl=5.1667642e-2*rhe**3-14.153121*rhe**2
     &               +1292.8377*rhe-3.9373536e4
         else
         if (rhe .gt. 80.) then
            waternacl=
     &      1.0629e-3*rhe**3-0.25281*rhe**2+20.171*rhe-5.3558e2
         else
         if (rhe .gt. 50.) then
            waternacl=
     &      4.2967e-5*rhe**3-7.3654e-3*rhe**2+.46312*rhe-7.5731
         else
         if (rhe .gt. 20.) then
            waternacl=
     &      2.9443e-5*rhe**3-2.4739e-3*rhe**2+7.3430e-2*rhe+1.3727
         else
            waternacl=1.17
         endif
         endif
         endif
         endif

         !check for error
         if (waternacl .gt. 45.) then
            write(*,*) 'ERROR in waternacl'
            write(*,*) rhe,waternacl
            call stop_model('ERROR in waternacl',255)
         endif

      RETURN
      END  FUNCTION waternacl


C     **************************************************
C     *  waterocil                                     *
C     **************************************************

C     MODIFIED BY YUNHA LEE, AUG, 2006

C     This function uses the current RH to calculate how much water is 
C     in equilibrium with the hydrophillic OA.  Aerosol water concentrations 
C     are assumed to be in equilibrium at all times and the array of 
C     concentrations is updated accordingly.

      real*8 FUNCTION waterocil(rhe)

      IMPLICIT NONE

C-----VARIABLE DECLARATIONS---------------------------------------------
      real*8 rhe   !relative humidity (0-100 scale)
      real*8 a, b, c, d, e, f, prefactor, activcoef
      parameter(a=1.0034, b=0.1614, c=1.1693,d=-3.1,
     & e=6.0)

C     VARIABLE COMMENTS...

C     waterocil is the ratio of wet mass to dry mass of a particle.  Instead
C     of calling a thermodynamic equilibrium code, this routine uses a
C     simple curve fit to estimate waterocil based on the current humidity.
C     The curve fit is based on observations of Dick et al. JGR D1 1471-1479

C-----CODE--------------------------------------------------------------
      
      if (rhe .gt. 99.) rhe=99.
      if (rhe .lt. 1.) rhe=1.

      if (rhe .gt. 85.) then
         waterocil=d+e*(rhe/100) 
cyhl Growth factor above RH 85% is not available, so it assumes linear growth 
cyhl at above 85%.  
      else
         waterocil=a+b*(rhe/100)+c*(rhe/100)**2. 
cyhl This eq is based on the extrapolation curve obtained from  
cyhl Dick et al 2000 figure 5.(High organic,density=1400g/cm3)
      endif
      
         !check for error
      if (waterocil .gt. 10.) then
         write(*,*) 'ERROR in waterocil'
         write(*,*) rhe,waterocil
         call stop_model('ERROR in waterocil',255)
      endif

      RETURN
      END  FUNCTION waterocil


C     **************************************************
C     *  GCM_ezwatereqm                                    *
C     **************************************************

C     WRITTEN BY Peter Adams, March 2000
C     MODIFIED BY Yunha Lee, March 2011 - to use this in GCM subroutine 


C     This routine uses the current RH to calculate how much water is 
C     in equilibrium with the aerosol.  Aerosol water concentrations 
C     are assumed to be in equilibrium at all times and the array of 
C     concentrations is updated accordingly.

C     This version of the routine works for sulfate and sea salt
C     particles.  They are assumed to be externally mixed and their
C     associated water is added up to get total aerosol water.
C     wr is the ratio of wet mass to dry mass of a particle.  Instead
C     of calling a thermodynamic equilibrium code, this routine uses a
C     simple curve fits to estimate wr based on the current humidity.
C     The curve fit is based on ISORROPIA results for ammonium bisulfate
C     at 273 K and sea salt at 273 K.


      SUBROUTINE aeroupdate

C-----INCLUDE FILES-----------------------------------------------------
      USE DOMAIN_DECOMP_ATM, only : GRID, write_parallel,
     &     am_i_root, getDomainBounds
      USE TOMAS_AEROSOL 
      USE GEOM, only: imaxj
      USE TRACER_COM, only : IDTSO4, IDTNA, IDTOCIL,IDTH2O,NBINS
     &     ,trm,IDTECOB,IDTECIL,IDTOCOB,IDTDUST,IDTNUMD,TRNAME
     *     ,ntm,ntm_TOMAS

      USE TRDIAG_COM, only : taijs=>taijs_loc !,taijls=>taijls_loc
      USE RESOLUTION, only : im,jm,lm     ! dimensions
      USE ATM_COM, only :   t            ! potential temperature (C)
     $                     ,q            ! saturated pressure
      USE MODEL_COM, only : dtsrc
      USE CONSTANT,   only:  lhe
      USE ATM_COM,   only: pmid,pk ! midpoint pressure in hPa (mb)
      IMPLICIT NONE

C-----VARIABLE DECLARATIONS---------------------------------------------

      INTEGER J_0, J_1, I_0, I_1
      integer i,j,l,k,n,jc,mpnum  !counters
      real*8 qsat         !used in RH calculation
      integer tracnum
      real*8 frac, Nkout(iBINS),Mkout(iBINS,icomp),Gcout(icomp-1)
      real*8 rhe
      real*8 waterso4, waternacl, waterocil


C-----CODE-----------------------------------------------------------    
C****
C**** Extract useful local domain parameters from "grid"
C****
      call getDomainBounds(grid, J_STRT=J_0,       J_STOP=J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP


C     Loop over all grid cells
      DO L=1,LM                            
         DO J=J_0,J_1                          
            DO I=I_0,IMAXJ(J)

               temp = pk(l,i,j)*t(i,j,l) !should be in [K]
               rh = MIN(1.,q(i,j,l)/QSAT(temp,lhe,pmid(l,i,j))) ! rH [0-100%]
C     Swap GCM variables into aerosol algorithm variables
               do n=1,NBINS
                  Nk(n)=trm(i,j,l,IDTNUMD-1+n)
                  Mk(n,srtso4)=trm(i,j,l,IDTSO4-1+n)
                  Mk(n,srtna )=trm(i,j,l,IDTNA -1+n)
                  Mk(n,srtnh4)=0.1875*Mk(n,srtso4) ! artificial for now.. 0.0!t0m(i,j,l,IDTNH4-1+n)
                  MK(n,srtecob)=trm(i,j,l,IDTECOB -1+n)
                  MK(n,srtecil)=trm(i,j,l,IDTECIL -1+n)
                  MK(n,srtocob)=trm(i,j,l,IDTOCOB -1+n)
                  MK(n,srtocil)=trm(i,j,l,IDTOCIL -1+n) 
                  MK(n,srtdust)=trm(i,j,l,IDTDUST -1+n) 
                  Mk(n,srth2o)= trm(i,j,l,IDTH2O-1+n) !I don't think this is necessary!
               enddo

C ****************
C Aerosol dynamics
C ****************

      !Do water eqm at appropriate times

      call storenm()
! I won't update Nk and Mk changes to aerodiag for now. 
      call mnfix(Nk,Mk) 
      mpnum=7
      call aerodiag(mpnum,i,j,l)

      call ezwatereqm(Mk)

C ***********************
C End of aerosol dynamics
C ***********************

C Swap Nk, Mk, and Gc arrays back to T0M
      do n=1,NBINS
         tracnum=IDTNUMD-1+n
         if (Nk(n) .ge. TRM(i,j,l,tracnum)) then
            TRM(i,j,l,tracnum)=Nk(n)
         else
            frac=Nk(n)/TRM(i,j,l,tracnum)
            call scalemom(i,j,l,tracnum,frac)
         endif
         do jc=1,icomp-idiag
            tracnum=IDTSO4-1+n+ibins*(jc-1)
            if (Mk(n,jc) .ge. TRM(i,j,l,tracnum)) then
               TRM(i,j,l,tracnum)=Mk(n,jc)
            else
               frac=Mk(n,jc)/TRM(i,j,l,tracnum)
               call scalemom(i,j,l,tracnum,frac)
            endif
         enddo
         tracnum=IDTH2O-1+n
         if (Mk(n,srth2o) .ge. TRM(i,j,l,tracnum)) then
            TRM(i,j,l,tracnum)=Mk(n,srth2o)
         else
            frac=Mk(n,srth2o)/TRM(i,j,l,tracnum)
            call scalemom(i,j,l,tracnum,frac)
         endif    
      enddo

C Check for negative tracer problems
      do n=1,ntm_TOMAS
         if (TRM(i,j,l,IDTSO4+n-1) .lt. 0.0) then
            if (abs(TRM(i,j,l,IDTSO4+n-1)) .gt. 1.e-10) then
               !serious problem - report error
               write(*,*) 'ERROR: Tracer ',trname(IDTSO4+n-1),
     &              trm(i,j,l,IDTSO4+n-1)
               write(*,*) ' < 0 in box ', i,j,l
               call stop_model('TRM<0 in aeroupdate',255)
!            else
!               !numerical problem - set to zero
!               TRM(i,j,l,IDTSO4+n-1)=0.0!1.d-42 !5??
            endif
         endif
      enddo

            enddo
         enddo
      enddo
      
      RETURN
      END SUBROUTINE aeroupdate


C     **************************************************
C     *  scalemom                                      *
C     **************************************************

C     WRITTEN BY Peter Adams, September 2000

C     When a tracer concentration decreases, call this routine to
C     decrease T0M and all higher order moments in proportion.

C-----INPUTS------------------------------------------------------------

C-----OUTPUTS-----------------------------------------------------------

      SUBROUTINE scalemom(i,j,l,tn,f)

      USE QUSDEF, only : nmom
      USE TRACER_COM, only : trm, trmom

C-----ARGUMENT DECLARATIONS---------------------------------------------

      integer i,j,l       !grid box
      integer tn          !tracer id number
      real*8 f  !factor (0-1) by which to decrease all moments

C-----VARIABLE DECLARATIONS---------------------------------------------

C     VARIABLE COMMENTS...

C-----ADJUSTABLE PARAMETERS---------------------------------------------

C-----CODE--------------------------------------------------------------

      TRM(i,j,l,tn) = TRM(i,j,l,tn)*f

      do n=1,nmom
         trmom(n,i,j,l,tn)=f*trmom(n,i,j,l,tn)
      enddo
      
      RETURN
      END SUBROUTINE scalemom


C     **************************************************
C     *  storenm                                       *
C     **************************************************

C     WRITTEN BY Peter Adams, June 2000

C     Stores values of Nk and Mk into Nkd and Mkd for diagnostic
C     purposes.  Also do gas phase concentrations.

C-----INPUTS------------------------------------------------------------

C-----OUTPUTS-----------------------------------------------------------

      SUBROUTINE storenm()

      USE TOMAS_AEROSOL

      IMPLICIT NONE

      integer j,k

C     VARIABLE COMMENTS...   

C-----ADJUSTABLE PARAMETERS---------------------------------------------

C-----CODE--------------------------------------------------------------

      do j=1,icomp-1
         Gcd(j)=Gc(j)
      enddo
      do k=1,ibins
         Nkd(k)=Nk(k)
         do j=1,icomp
            Mkd(k,j)=Mk(k,j)
         enddo
      enddo

      RETURN
      END SUBROUTINE storenm


C     **************************************************
C     *  aerodiag                                      *
C     **************************************************

C     WRITTEN BY Peter Adams, June 2000

C     Accumulates diagnostics on aerosol microphysical processes.

C-----INPUTS------------------------------------------------------------

C-----OUTPUTS-----------------------------------------------------------

      SUBROUTINE aerodiag(pt,i,j,l)

      USE TOMAS_AEROSOL
      USE TRACER_COM, only : IDTSO4,IDTNUMD,n_H2SO4
      IMPLICIT NONE

C-----VARIABLE DECLARATIONS---------------------------------------------

      integer pt, i, j, l, jc, n
      integer k,kk,tracnum

C-----CODE--------------------------------------------------------------

      if(pt.eq.4.or.pt.eq.5.or.pt.eq.7)then ! Aqoxid_mc,Aqoxid_ls,Aeroupdate

C Bulk species
        AEROD(i,j,l,n_H2SO4,pt)=
     &     (Gc(srtso4)-Gcd(srtso4))

         do n=1,ibins       
                       
!     Aerosol number
            tracnum=IDTNUMD-1+n
            AEROD(i,j,l,tracnum,pt)= 
     &           (Nk(n)-Nkd(n))
!     Aerosol mass
            do jc=1,icomp-idiag
               tracnum=IDTSO4-1+n+ibins*(jc-1)
               AEROD(i,j,l,tracnum,pt)= 
     &              (Mk(n,jc)-Mkd(n,jc))
               
            enddo   
         enddo
         
      else

C Bulk species
        AEROD(i,j,l,n_H2SO4,pt)=AEROD(i,j,l,n_H2SO4,pt)
     &     +(Gc(srtso4)-Gcd(srtso4))

         do n=1,ibins  
!     Aerosol number
            tracnum=IDTNUMD-1+n
            
            AEROD(i,j,l,tracnum,pt)= 
     &           AEROD(i,j,l,tracnum,pt) 
     &           +(Nk(n)-Nkd(n))  
            
!     Aerosol mass
            do jc=1,icomp-idiag
               tracnum=IDTSO4-1+n+ibins*(jc-1)
               AEROD(i,j,l,tracnum,pt)= 
     &           AEROD(i,j,l,tracnum,pt) 
     &              + (Mk(n,jc)-Mkd(n,jc))
               
            enddo 
         enddo
      endif     
 
      RETURN
      END SUBROUTINE aerodiag


C     **************************************************
C     *  subgridcoag_drv                                   *
C     **************************************************

C     Written by Yunha Lee, July 2011 
C     In order to accommodate 3-D emission in GISS ModelE, the original 
C     subgridcoag subroutine, written by Jeff Pierce, is modified into two. 
C 
C     No moments updated here because aerosol emission are positive! 

      SUBROUTINE subgridcoag_drv(dtstep)

C-----INCLUDE FILES--------------------------------------------------

      USE DOMAIN_DECOMP_ATM, ONLY : GRID,getDomainBounds,write_parallel
      USE TOMAS_AEROSOL
      USE TRDIAG_COM, only : taijs=>taijs_loc,ijts_subcoag,itcon_subcoag
      USE FLUXES, only : tr3dsource

      use resolution, only     : lm
      USE ATM_COM, only :   t   ! potential temperature (C)
     $     ,q                   ! saturated pressure
     $     ,pmid                ! midpoint pressure in hPa (mb)
     $     ,pk                  ! pk is t mess up factor
     $     ,am                 ! BYAM  1/Air mass (m^2/kg) 
      
      USE MODEL_COM, only : dtsrc
      USE GEOM, only : imaxj,axyp,BYAXYP
      USE CONSTANT, ONLY : pi,gasc,mair 

      USE TRACER_COM, only : nbins,xk,ntm,trm,trmom,ntsurfsrc,
     &     IDTSO4,IDTNA,IDTECOB,IDTECIL,IDTOCOB,
     &     IDTOCIL,IDTDUST,IDTNUMD,n_SO2,IDTH2O
 
      IMPLICIT NONE

C-----PASSED VARIABLE DECLARATIONS-----------------------------------
      integer :: J_1, J_0, I_1, I_0
      INTEGER :: L,I,J

      REAL*8, INTENT(IN) :: dtstep
      INTEGER n,ns,c,k,tot_src,tracnum
      INTEGER tomas_ntsurf !same as ntsurfsrc
      real*8 ndistinit(nbins) !the number of particles being added to the gridbox before subgrid coag
      real*8,dimension(nbins) ::  ndist, ndist2, ndist0 !the number of particles in the box
      real*8,dimension(nbins,icomp) :: mdist,mdist2,mdist0 ! the mass of each component in the box. (kg)
      real tscale ! the scale time for mixing (s)
      real*8 ndistfinal(nbins),tot_ndistinit(nbins) !the number of particles being added to the gridbox after subgrid coag
      real*8 maddfinal(nbins) !the mass that should be added to each bin due to coagulation (kg)

C-----CODE-----------------------------------------------------------

      call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

      DO L=1,LM; DO J=J_0,J_1; DO I=I_0,imaxj(j)

c$$$        IF(I.EQ.25.AND.J.EQ.62)THEN
c$$$          open (1044,file='debug_coag.dat',access='append',
c$$$     &         status='unknown')
c$$$          do k=1,nbins
c$$$            write(1044,*)'begin',l,trm_emis(i,j,l,idtecob+k-1)
c$$$     $           ,trm_emis(i,j,l,idtnumd+k-1),
c$$$     $   trm(i,j,l,idtecob+k-1),trm(i,j,l,idtnumd+k-1),
c$$$     $           taijs(i,j,ijts_subcoag(IDTECOB+k-1))
c$$$            enddo
c$$$      ENDIF

      tscale=5.*3600.
      temp = pk(l,i,j)*t(i,j,l) !should be in [K]
      pres= pmid(l,i,j)*100.    ! pmid in [hPa]
      boxvol=am(l,i,j)*axyp(i,j)/mair*1000.d0
     &     *gasc*temp/pres*1e6  !cm3

      do k=1,nbins
        ndist0(k)=TRM_EMIS(I,J,L,IDTNUMD+k-1)
        do c=1,icomp-idiag
          mdist0(k,c)=TRM_EMIS(I,J,L,IDTSO4+(c-1)*nbins+k-1)
        enddo
        mdist0(k,srtnh4)=0.0
        mdist0(k,srth2o)=TRM_EMIS(I,J,L,IDTH2O+k-1)
        ndistfinal(k)=0
        maddfinal(k)=0
      enddo
      
      ndist(:)=0.
      mdist(:,:)=0.
      ndistinit(:)=0.
        
      DO NS=1,3
!        print*,'subcoag 3d',ns,dtstep,nbins,idtso4

        do k=1,nbins
          if(ns.lt.3) 
     &         ndistinit(k)=tr3Dsource(i,j,l,ns,IDTNUMD+K-1)*dtstep
          if(ns.eq.3) 
     &         ndistinit(k)=tr3Dsource(i,j,l,ns+1,IDTNUMD+K-1)*dtstep 
        enddo
        
c$$$        IF(I.EQ.25.AND.J.EQ.62)THEN
c$$$          
c$$$          do k=1, nbins
c$$$          write(1044,*)'middle',l,ns,ndist0(k),mdist0(k,srtecob),
c$$$     $         ndist(k),mdist(k,srtecob),ndist2(k),mdist2(k,srtecob)   
c$$$          enddo
c$$$        ENDIF
c$$$
        if(sum(ndistinit(1:nbins)).gt.0.)then
!     only when there is emission! 
          call subgridcoag(ndistinit,ndist0,mdist0,boxvol,
     &         tscale,ndistfinal,maddfinal) ! account for subgrid coagulation
          
          do k=1,nbins
            ndist(k)=ndist(k)+ndistfinal(k)
            
            IF(NS.EQ.1)THEN
                                !SO4
              mdist(k,srtso4)=
     &             ndistfinal(k)*(sqrt(xk(k)*xk(k+1)))+
     &             maddfinal(k)
              
            ELSEIF(NS.EQ.2)THEN
                                !EC
              mdist(k,srtecil)=
     &             ndistfinal(k)*0.2*(sqrt(xk(k)*xk(k+1)))+
     &             maddfinal(k)*0.2
              mdist(k,srtecob)=
     &             ndistfinal(k)*0.8*(sqrt(xk(k)*xk(k+1)))+
     &             maddfinal(k)*0.8   
              
            ELSEIF(NS.EQ.3)THEN
                                !OC
                mdist(k,srtocil)=
     &               ndistfinal(k)*0.5*(sqrt(xk(k)*xk(k+1)))+
     &               maddfinal(k)*0.5
                mdist(k,srtocob)=
     &               ndistfinal(k)*0.5*(sqrt(xk(k)*xk(k+1)))+
     &               maddfinal(k)*0.5   
              ENDIF
            enddo               !k
          endif
        enddo                   !ns

!     fix the inconsistancies in the distribution
        do k=1,nbins
          ndist2(k)=ndist(k)+ndist0(k)
          do c=1,icomp !-idiag
            mdist2(k,c)=mdist(k,c)+mdist0(k,c)
          enddo
        enddo
 
c$$$        IF(I.EQ.25.AND.J.EQ.62)THEN
c$$$          do k=1,nbins
c$$$            write(1044,*)'middle 2',l,ndist(k),ndist2(k)
c$$$     $           ,mdist(k,srtecob),mdist2(k,srtecob) 
c$$$          enddo
c$$$
c$$$      ENDIF
        if(sum(ndist(1:nbins)).gt.0.) call mnfix(ndist2,mdist2)
        
c$$$        IF(I.EQ.25.AND.J.EQ.62)THEN
c$$$          do k=1,nbins
c$$$            write(1044,*)'after mnfix',l,ndist(k),ndist2(k)
c$$$     $           ,mdist(k,srtecob),mdist2(k,srtecob),sum(ndist(1:nbins)) 
c$$$          enddo   
c$$$
c$$$      ENDIF
!        call mnfix(ndist2,mdist2)
!     DIAGNOSTICS! 
        
        do k=1,nbins  
          
          tracnum=IDTNUMD-1+k  
          N_subgridcg(i,j,l,k,2)=(ndist2(k)- !this is emission after subgrid
     &         trm(i,j,l,tracnum)) 


!          if(ndist2(k)/trm(i,j,l,tracnum).GT.1000.)THEN
!            PRINT*,'too large subcoag',i,j,l,k,ndist2(k),ndist(k),
!     $           ndist0(k)
!          endif

          if(l.eq.1)then
            N_subgridcg(i,j,l,k,2)=N_subgridcg(i,j,l,k,2)+ 
     &           N_subgridcg(i,j,l,k,1) !from 2-d emission subgrid coagulation
          endif
          
          trm(i,j,l,tracnum)=ndist2(k)
          
          taijs(i,j,ijts_subcoag(tracnum)) 
     &         =taijs(i,j,ijts_subcoag(tracnum))
     &         +N_subgridcg(i,j,l,k,2) ! /adt
          
          if (itcon_subcoag(tracnum).gt.0) 
     &         call inc_diagtcb(i,j,N_subgridcg(i,j,l,k,2) ,
     &         itcon_subcoag(tracnum),tracnum)
          
          do c=1,icomp-idiag
            
            tracnum=IDTSO4-1+k+nbins*(c-1) 
            M_subgridcg(i,j,l,k,c,2)=mdist2(k,c)- !trm + emission after subgrid 
     &           trm(i,j,l,tracnum) !trm + emission before subgrid (which is computed in apply_tracer3d)


          if(l.eq.1)then
            M_subgridcg(i,j,l,k,c,2)=M_subgridcg(i,j,l,k,c,2)+
     &           M_subgridcg(i,j,l,k,c,1)
          endif

            trm(i,j,l,tracnum)=mdist2(k,c)
            
            taijs(i,j,ijts_subcoag(tracnum)) 
     &           =taijs(i,j,ijts_subcoag(tracnum))
     &           +M_subgridcg(i,j,l,k,c,2) ! /adt

            if (itcon_subcoag(tracnum).gt.0) 
     &           call inc_diagtcb(i,j,M_subgridcg(i,j,l,k,c,2),
     &           itcon_subcoag(tracnum),tracnum)
            
            
          enddo
        enddo

c$$$        IF(I.EQ.25.AND.J.EQ.62)THEN
c$$$          do k=1,nbins
c$$$            
c$$$        write(1044,*) 'subcoag result',l,trm_emis(i,j,l,idtecob+k-1),
c$$$     $       trm(i,j,l,idtecob+3),M_subgridcg(i,j,l,k,3,2),
c$$$     $           N_subgridcg(i,j,l,k,2),
c$$$     $       taijs(i,j,ijts_subcoag(IDTECOB+k-1))
c$$$      enddo 
c$$$
c$$$      ENDIF

        M_subgridcg(i,j,l,:,:,:)=0.0
        N_subgridcg(i,j,l,:,:)=0.0
        
       enddo; enddo; enddo
        
       return
       end subroutine subgridcoag_drv
    
C     **************************************************
C     *  subgridcoag_drv                                   *
C     **************************************************

C     Written by Yunha Lee, July 2011 
C     In order to accommodate 2-D and 3-D emission in GISS ModelE, the original 
C     subgridcoag subroutine, written by Jeff Pierce, is modified into two. 
C 
C     No moments updated here because aerosol emission are positive! 

      SUBROUTINE subgridcoag_drv_2D(dtstep)

C-----INCLUDE FILES--------------------------------------------------

      USE DOMAIN_DECOMP_ATM, ONLY : GRID,getDomainBounds,write_parallel
      USE TOMAS_AEROSOL
      USE FLUXES, only : trsource,atmsrf

      use resolution, only     : lm
      USE ATM_COM, only :   t   ! potential temperature (C)
     $     ,q                   ! saturated pressure
     $     ,pmid                ! midpoint pressure in hPa (mb)
     $     ,pk                  ! pk is t mess up factor
     $     ,am                 ! BYAM  1/Air mass (m^2/kg) 
      
      USE MODEL_COM, only : dtsrc
      USE GEOM, only : imaxj,axyp,BYAXYP
      USE CONSTANT, ONLY : pi,gasc,mair 

      USE TRACER_COM, only : nbins,xk,ntm,trm,trmom,ntsurfsrc,
     &     IDTSO4,IDTNA,IDTECOB,IDTECIL,IDTOCOB,
     &     IDTOCIL,IDTDUST,IDTNUMD,n_SO2,IDTH2O
 
      IMPLICIT NONE

C-----VARIABLE DECLARATIONS-----------------------------------

      integer :: J_1, J_0, I_1, I_0,L,I,J

      REAL*8, INTENT(IN) :: dtstep
      INTEGER n,ns,c,k,tot_src,tracnum
      INTEGER tomas_ntsurf !same as ntsurfsrc
      real*8 ndistinit(nbins) !the number of particles being added to the gridbox before subgrid coag
      real*8,dimension(nbins) ::  ndist, ndist0 !the number of particles in the box
      real*8,dimension(nbins,icomp) :: mdist,mdist0 ! the mass of each component in the box. (kg)
      real*8 tscale ! the scale time for mixing (s)
      real*8 ndistfinal(nbins),tot_ndistinit(nbins) !the number of particles being added to the gridbox after subgrid coag
      real*8 maddfinal(nbins) !the mass that should be added to each bin due to coagulation (kg)
      real*8 trflux_tom ! temporal variable

C-----CODE-----------------------------------------------------------

      call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP


      l=1 !2-D emission is at 1st layer

      DO J=J_0,J_1
        DO I=I_0,imaxj(j)
          
!     subgrid timescale and met conditions
 
        tscale=5.*3600.
        temp = pk(l,i,j)*t(i,j,l) !should be in [K]
        pres= pmid(l,i,j)*100.  ! pmid in [hPa]
        boxvol=am(l,i,j)*axyp(i,j)/mair*1000.d0
     &       *gasc*temp/pres*1e6 !cm3
        
!     Amount of tracer before emission is applied. 
        
        do k=1,nbins
          ndist0(k)=TRM(I,J,L,IDTNUMD+k-1)
          do c=1,icomp-idiag
            mdist0(k,c)=TRM(I,J,L,IDTSO4+(c-1)*nbins+k-1)
          enddo
          mdist0(k,srtnh4)=0.0
          mdist0(k,srth2o)=TRM(I,J,L,IDTH2O+k-1)
          ndistfinal(k)=0
          maddfinal(k)=0
        enddo
              
!     Only initialize when 2-D emission starts! 
        tot_ndistinit(:)=0.
        ndist(:)=0.
        mdist(:,:)=0.
        ndistinit(:)=0.
                
        DO ns=1,ntsurfsrc(idtnumd)
!          print*,'subcoag 2d',ns,ntsurfsrc(idtnumd),dtstep
!     ns=1 for so4; ns=2 for ec; ns=3 for oc
          
          do k=1,nbins
            
            ndistinit(k)=trsource(i,j,NS,IDTNUMD+K-1)*dtstep
            
            tot_ndistinit(k)=tot_ndistinit(k)+ndistinit(k) !sum of number emission for SO4, EC, and OC
          enddo
          
          if(sum(ndistinit(1:nbins)).gt.0.)then
!     only when there is emission! 
            
            call subgridcoag(ndistinit,ndist0,mdist0,boxvol,
     &           tscale,ndistfinal,maddfinal) ! account for subgrid coagulation
            
            do k=1,nbins

              ndist(k)=ndistfinal(k)+ ndist(k) !sum of all emission type after subgrid coag
              
              IF(NS.EQ.1)       !SO4
     &             mdist(k,srtso4)= 
     &             ndistfinal(k)*(sqrt(xk(k)*xk(k+1)))+
     &             maddfinal(k)
              
              IF(NS.EQ.2)       !EC
     &             mdist(k,srtecil)= 
     &             ndistfinal(k)*0.2*(sqrt(xk(k)*xk(k+1)))+
     &             maddfinal(k)*0.2                
              IF(NS.EQ.2)       !EC
     &             mdist(k,srtecob)= 
     &             ndistfinal(k)*0.8*(sqrt(xk(k)*xk(k+1)))+
     &             maddfinal(k)*0.8   
              
              IF(NS.EQ.3)       !OC
     &             mdist(k,srtocil)=
     &             ndistfinal(k)*0.5*(sqrt(xk(k)*xk(k+1)))+
     &             maddfinal(k)*0.5
              IF(NS.EQ.3)       !OC
     &             mdist(k,srtocob)=
     &             ndistfinal(k)*0.5*(sqrt(xk(k)*xk(k+1)))+
     &             maddfinal(k)*0.5 

            enddo

          endif !positive emission
          
        enddo !ntsurfsrc    
        
c$$$!     fix the inconsistancies in the distribution
c$$$        do k=1,nbins
c$$$          ndist2(k)=ndist(k)+ndist0(k)
c$$$          do c=1,icomp-idiag
c$$$            mdist2(k,c)=mdist(k,c)+mdist0(k,c)
c$$$          enddo
c$$$        enddo
        
        if(sum(ndist(1:nbins)).gt.0.) call mnfix(ndist,mdist)
        
!     DIAGNOSTICS! 
      
        do k=1,nbins  
          
          tracnum=IDTNUMD-1+k  
          
          N_subgridcg(i,j,l,k,1)=N_subgridcg(i,j,l,k,1)+ndist(k)- ! emission after subgrid
     &         tot_ndistinit(k) ! emission before subgrid
          
!          trflux_tom=atmsrf%trflux_prescr(tracnum,i,j)
!     &         /byaxyp(i,j)+
!     &         (ndist(k)- tot_ndistinit(k))/dtstep !kg/sec

          atmsrf%trflux_prescr(tracnum,i,j)= ndist(k)/dtstep*byaxyp(i,j) ! kg/s to kg/m2/s
          
          do c=1,icomp-idiag
            
            tracnum=IDTSO4-1+k+nbins*(c-1)
            
            if(c.eq.2.or.c.eq.7)then
                                !no subgrid coagulation
              M_subgridcg(i,j,l,k,c,1) =0.

            else
              
              if(c.eq.1) tomas_ntsurf=ntsurfsrc(n_SO2)
              if(c.eq.3.or.c.eq.4) tomas_ntsurf=ntsurfsrc(IDTECOB) !ecob
              if(c.eq.5.or.c.eq.6)  tomas_ntsurf=ntsurfsrc(IDTOCOB) !ecob
              
              M_subgridcg(i,j,l,k,c,1)=M_subgridcg(i,j,l,k,c,1)
     &             + mdist(k,c)-
     &             atmsrf%trflux_prescr(tracnum,i,j)/byaxyp(i,j)*dtstep
!     &             (sum(trsource(i,j,1:tomas_ntsurf,tracnum))*dtstep) ! emission before subgrid
              
            endif
!            trflux_tom=atmsrf%trflux_prescr(tracnum,i,j)
!     &           /byaxyp(i,j)+ (mdist(k,c)/dtstep-
!     &             atmsrf%trflux_prescr(tracnum,i,j)/byaxyp(i,j))

!     &           (sum(trsource(i,j,1:tomas_ntsurf,tracnum))))


            atmsrf%trflux_prescr(tracnum,i,j)=
     &           mdist(k,c)/dtstep*byaxyp(i,j)
            
          enddo !c
        enddo ! k

      enddo

!      if(M_subgridcg(i,j,l,4,3,1).gt.0.) 
!     $     print*,'subcoag 2D result',trsource(i,j,2,IDTNUMD+3),
!     $       trm(i,j,1,idtecob+3),M_subgridcg(i,j,l,4,3,1)
      
      enddo ! i,j
      
      return
      end subroutine subgridcoag_drv_2D
     


C     **************************************************
C     *  subgridcoag                                   *
C     **************************************************
C     WRITTEN BY Jeff Pierce, December, 2006

C     This subroutine will determine how much of each size of freshly emitted aerosol will 
C     be scavenged by coagulation prior to being completely mixed in the gridbox and will
C     give the new emissions size distribution along with where the mass of coagulated
C     particles should be added.

      SUBROUTINE subgridcoag(ndistinit,ndist2,mdist2,boxvolume,
     & tscale,ndistfinal,maddfinal)

C-----INCLUDE FILES--------------------------------------------------
      USE TRACER_COM, only : nbins,xk

      USE TOMAS_AEROSOL
      USE CONSTANT, ONLY : pi,gasc,mair
      IMPLICIT NONE

C-----PASSED VARIABLE DECLARATIONS-----------------------------------
 
      INTEGER n,k,c,kk
      real*8, intent(in) :: ndistinit(nbins) !the number of particles being added to the gridbox before subgrid coag
      real*8, intent(in) :: ndist2(nbins) !the number of particles in the box
      real*8, intent(in) :: mdist2(nbins,icomp) ! the mass of each component in the box. (kg)
      real*8, intent(in) :: boxvolume  ! volume of box in cm3
      real, intent(in) :: tscale ! the scale time for mixing (s)
      real*8, intent(out) :: ndistfinal(nbins) !the number of particles being added to the gridbox after subgrid coag
      real*8, intent(out) :: maddfinal(nbins) !the mass that should be added to each bin due to coagulation (kg)

C-----VARIABLE DECLARATIONS------------------------------------------

      real*8 mp ! mass of the particle (kg)
      real density                !density (kg/m3) of particles
      real*8 diameter(nbins) ! diamter of the particle (m)
      real*8 diaml(nbins) ! total diamter of particles larger (m/cm3)
      real*8 fracdiaml(nbins,nbins) ! fraction of coagulation that occurs with each bin larger
      real*8 kcoag(nbins) ! the coagulation rate for the particles in each bin (s^-1)
      real aerodens

      real mso4, mh2o, mno3, mnh4  !mass of each component (kg/grid box)
      real mecil,mecob,mocil,mocob
      real mdust,mnacl   

!      external aerodens

C-----VARIABLE COMMENTS----------------------------------------------

C-----ADJUSTABLE PARAMETERS------------------------------------------

      real*8 v1,v2,v3  !for coag rate calculation
      parameter(v1=8.5708E-13,
     &          v2=-1.4174,
     &           v3=4.3047E-4)

C-----CODE-----------------------------------------------------------

C     get the wet diameter of particles in each size bin
      do k=1,nbins
         mp=0.1875*mdist2(k,srtso4)
         do c=1,icomp
            mp = mp + mdist2(k,c)
         enddo
         if (ndist2(k).eq.0.)then
            mp=sqrt(xk(k)*xk(k+1))
         else
            mp = mp / ndist2(k)
         endif
         if((mdist2(k,srtso4)+mdist2(k,srtna)+mdist2(k,srtocil)+
     &        mdist2(k,srtdust)).eq.0)then
            density=1400.
         else
         mso4=mdist2(k,srtso4) 
         mnacl=mdist2(k,srtna)
         mno3=0.e0
         if ((mso4+mno3) .lt. 1.e-8) mso4=1.e-8
         mnh4=0.1875*mso4  !assume ammonium bisulfate
         mecob=mdist2(k,srtecob)
         mecil=mdist2(k,srtecil)
         mocil=mdist2(k,srtocil)
         mocob=mdist2(k,srtocob)
         mdust=mdist2(k,srtdust)          
         mh2o=mdist2(k,srth2o)   

         density=aerodens(mso4,mno3,mnh4 !mno3 taken off!
     *        ,mnacl,mecil,mecob,mocil,mocob,mdust,mh2o) !assume bisulfate 
         endif
         diameter(k)=2.*(3./4./pi*mp/density)**(1./3.) ! m
      enddo

C     get the total diameter of particles larger than each size bin
      diaml(nbins)=0. !no diameter larger than largest bin
      do kk=1,nbins-1
         k=nbins-kk
         diaml(k) = diaml(k+1) + diameter(k+1)*ndist2(k+1)/boxvolume ! m/cm3
      enddo
      
C     get the fraction of the diameter larger that comes from each bin larger
      do k=1,nbins
         do kk=1,nbins
            fracdiaml(k,kk)=0.
         enddo
      enddo
      do k=1,nbins-1
         do kk=k+1,nbins
            if (diaml(k).ne.0.0)then
              fracdiaml(k,kk)=diameter(kk)*ndist2(kk)/boxvolume/diaml(k)
            else
               fracdiaml(k,kk)=0.0
            endif
         enddo
      enddo

C     determine the coagulation rate for each size bin
      do k=1,nbins
         if (diameter(k).gt.0.d0)then
            kcoag(k) = (v1*diameter(k)**(v2)+v3)*diaml(k)
         else
            kcoag(k) = 0.d0
         endif
      enddo

C     determine the number of new particles left after coagulation
      do k=1,nbins
         ndistfinal(k)=ndistinit(k)*exp(-kcoag(k)*tscale)
      enddo

C     determine the mass added to each bin coagulation
      do k=1,nbins
         maddfinal(k)=0.
      enddo
      do k=1,nbins-1
         do kk=k+1,nbins
            maddfinal(kk)=maddfinal(kk) + (ndistinit(k)-ndistfinal(k))*
     &           fracdiaml(k,kk)*sqrt(xk(k)*xk(k+1))
         enddo
      enddo

      return
      end SUBROUTINE subgridcoag

      subroutine alloc_tracer_TOMAS_com(grid)
!@SUM  To alllocate arrays whose sizes now need to be determined
!@+    at run-time
!@auth Yunha Lee
      use domain_decomp_atm, only : dist_grid, getDomainBounds 
      use resolution, only     : lm
      use TOMAS_aerosol

      IMPLICIT NONE

      type (dist_grid), intent(in) :: grid
      integer :: ier, J_1H, J_0H, I_1H, I_0H
      logical :: init = .false.

      if(init)return
      init=.true.
    
      call getDomainBounds( grid , J_STRT_HALO=J_0H, J_STOP_HALO=J_1H )
      I_0H=GRID%I_STRT_HALO
      I_1H=GRID%I_STOP_HALO 

      allocate(  AQSO4oxid_mc(I_0H:I_1H,J_0H:J_1H,LM)   )
      allocate(  AQSO4oxid_ls(I_0H:I_1H,J_0H:J_1H,LM)   )
      allocate(  H2SO4_chem(I_0H:I_1H,J_0H:J_1H,LM)  )
      allocate(  AEROD(I_0H:I_1H,J_0H:J_1H,LM,NTM,ptype) )
      allocate(  N_subgridcg(I_0H:I_1H,J_0H:J_1H,LM,IBINS,2) )
      allocate(  M_subgridcg(I_0H:I_1H,J_0H:J_1H,LM,IBINS,
     *     ICOMP-IDIAG,2))
     
      allocate(  TRM_EMIS(I_0H:I_1H,J_0H:J_1H,LM,NTM) )

      return
      end subroutine alloc_tracer_TOMAS_com
