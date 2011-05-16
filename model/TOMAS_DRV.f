#include "rundeck_opts.h"

      MODULE TOMAS_AEROSOL

      USE MODEL_COM, only : im,jm,lm     ! dimensions
      USE TRACER_COM, only : ntm
      IMPLICIT NONE 

C-----INCLUDE FILES--------------------------------------------------
 
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
      integer bin_nuc, tern_nuc, ion_nuc   !flags for binary and ternary nuc
      parameter(bin_nuc=1, tern_nuc=0, ion_nuc=0) ! 1 = on

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
      integer, parameter :: ptype=6 ! number of microphysics process 
      real*8, ALLOCATABLE,dimension(:,:,:,:,:) :: AEROD
      real*8, ALLOCATABLE,DIMENSION(:,:,:) :: AQSO4oxid_mc,AQSO4oxid_ls  !1 for Convective and 2 for large-scale
      real*8, ALLOCATABLE,DIMENSION(:,:,:)  ::  h2so4_chem  !h2so4 formation rate from so2+oh [kg of H2SO4/sec\
      integer ncomp
      parameter(ncomp=8)


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

      USE DOMAIN_DECOMP_ATM, only : GRID, GET, write_parallel
     &     ,am_i_root
      USE TOMAS_AEROSOL 
      USE TRACER_COM

      USE TRDIAG_COM, only : taijs=>taijs_loc,taijls=>taijls_loc
     *     ,ijts_TOMAS,itcon_TOMAS !ijlt_AMPm,ijlt_AMPext,ijts_AMPpdf
c$$$     *     ,itcon_AMP,itcon_AMPm
!      USE AEROSOL_SOURCES, only: off_HNO3
      USE FLUXES, only: tr3Dsource

      USE MODEL_COM, only : im,jm,lm     ! dimensions
     $                     ,t            ! potential temperature (C)
     $                     ,q            ! saturated pressure
     $                     ,dtsrc
      USE GEOM, only: axyp,imaxj,BYAXYP
      USE CONSTANT,   only:  lhe,mair,gasc   
      USE DYNAMICS,   only: pmid,pk,byam,gz, am   ! midpoint pressure in hPa (mb)
!                                           and pk is t mess up factor
!                                           BYAM  1/Air mass (m^2/kg)

c$$$      USE PBLCOM,     only: EGCM !(LM,IM,JM) 3-D turbulent kinetic energy [m^2/s^2]

c$$$#ifndef NO_HDIURN
c$$$c for the hourly diagnostic
c$$$#ifdef CLD_AER_CDNC 
c$$$      USE CLOUDS_COM, only: CDN3D  ! CDNC
c$$$#endif
c$$$      USE DIAG_COM, only: adiurn=>adiurn_loc,ndiuvar,iwrite,
c$$$     *     jwrite,itwrite,ndiupt,idd_diam
c$$$     *                    ,ijdd, idd_ccn, idd_cdnc, 
c$$$     *     idd_lwp, idd_numb, idd_mass, idd_so2,
c$$$     *                     idd_lwc, idd_ncL
c$$$     *     ,hdiurn=>hdiurn_loc
c$$$#endif
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
!      real*8, dimension(ntm) :: tomas_init
!      real*4 Mocob
             
      real*8 H2SO4rate_o ! H2SO4rate for the specific gridcell
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


c$$$      ! variables for recording timesteps in cond_nuc
c$$$      real tbinlimits(101) ! the time limits for each bin
c$$$      integer stepsinbin(100) ! the number of time steps in each bin
c$$$      real avgstep ! the average timestep taken during cond_nuc
c$$$      integer intTAU


C-----CODE-----------------------------------------------------------    
C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT=J_0,       J_STOP=J_1)
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

               adt=dtsrc        ! 30 mins

               H2SO4rate_o = H2SO4_chem(i,j,l) !kg of h2so4/sec  (from SO2+OH)

                                !calculate SOA to condense
c$$$               if(l.eq.1)then
c$$$                  SOArate =  SOA_chem(i,j)*(1.d0-
c$$$     &                 exp(-adt/(tau_soa*3600.*24.))) ! SOA_chem is kg/s 
c$$$               else
c$$$                  SOArate=0.0
c$$$               endif

               SOArate = TRM(i,j,l,n_SOAgas)*(1.d0-
     &              exp(-adt/(tau_soa*3600.*24.)))/adt

               if(H2SO4rate_o.le.0.) THEN
                  if(am_i_root())
     &                 print*,'problem in soa_amp',i,j,l,
     &                 h2so4rate_o,SOArate

                  soa_amp=0
               else

               soa_amp = SOArate/H2SO4rate_o  ! TOMAS- floating invalid due to zero division??
               endif
              
               
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
                  Mk(n,srtnh4)=0.! 1875*Mk(n,srtso4) !TRM(i,j,l,IDTNH4A-1+n) - TOMAS -WILL FIX IT SOON. 
               enddo

               INIT_NK(:) = NK(:)
               INIT_MK(:,:)=MK(:,:)
               INIT_H2SO4 = H2SO4_chem(I,J,L)*adt
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
               

C     ****************
C     Aerosol dynamics
C     ****************               
!               do mt=1,10        ! 6 * 10 min inside
                  
                                !Do water eqm at appropriate times
               call ezwatereqm(Mk)
               

               call storenm()
               call mnfix(Nk,Mk)
               mpnum=6 
               call aerodiag(mpnum,i,j,l)

c$$$               if(am_i_root())then
c$$$                  print*,'i',i,'j',j,'l',l
c$$$                  print*,'mnfix NK',Nk(1),Nkd(1)
c$$$                  print*,'mnfix MK',Mk(8,srtso4),Mkd(8,srtso4),
c$$$     *                 INIT_Mk(8,srtso4)
c$$$               endif

!! in-cloud oxidation! move from clouds2.f to here            
               call storenm()
               if(AQSO4oxid_mc(I,J,L).gt.0.)then
                  call aqoxid(AQSO4oxid_mc(I,J,L),.TRUE.,Nk,Mk,
     &                 Nkout,Mkout) ! Moist Convective clouds
                  Nk(:)=Nkout(:)
                  Mk(:,:)=Mkout(:,:)
! Nk/Mk changes by mnfix will go to aqoxid. 
                  call mnfix(Nk,Mk)
               endif

               mpnum=4 
               call aerodiag(mpnum,i,j,l)

c$$$               if(am_i_root())then
c$$$                  print*,'i',i,'j',j,'l',l
c$$$                  print*,'aq mc NK',Nk(1),Nkd(1)
c$$$                  print*,'aq mc MK',Mk(8,srtso4),Mkd(8,srtso4)
c$$$               endif


!Should I call aerodiag?  Can I call more than once?

               call storenm()
               if(AQSO4oxid_ls(I,J,L).gt.0.)then
                  call aqoxid(AQSO4oxid_ls(I,J,L),.false.,Nk,Mk,
     &                 Nkout,Mkout) ! Large scale clouds
                  Nk(:)=Nkout(:)
                  Mk(:,:)=Mkout(:,:)
                  call mnfix(Nk,Mk)
               endif

               mpnum=5 
               call aerodiag(mpnum,i,j,l)
c$$$               if(am_i_root())then
c$$$                  print*,'i',i,'j',j,'l',l
c$$$                  print*,'aq ls NK',Nk(1),Nkd(1)
c$$$                  print*,'aq ls MK',Mk(8,srtso4),Mkd(8,srtso4)
c$$$               endif

! get the total mass of S
               tot_s_1 = H2SO4rate_o*adt*32.d0/98.d0
               do k=1,ibins
                  tot_s_1 = tot_s_1 + Mk(k,srtso4)*32.d0/96.d0
               enddo

               Gcavg = 0.0
               call storenm()

C If any Nk are zero, then set them to a small value to avoid division by zero
               call cond_nuc(Nk,Mk,Gc,Nkout,Mkout,Gcout,fn,fn1,
     &             H2SO4rate_o,adt,num_iter,Nknuc,Mknuc,Nkcond,Mkcond)
                                !get nucleation diagnostic
 
               Mk(:,:)=Mknuc(:,:)
               Nk(:)=Nknuc(:)
            
               mpnum=3 
               call aerodiag(mpnum,i,j,l)

c$$$               if(am_i_root())then
c$$$                  print*,'i',i,'j',j,'l',l
c$$$                  print*,'nuc NK',Nk(1),Nkd(1)
c$$$                  print*,'nuc MK',Mk(8,srtso4),Mkd(8,srtso4)
c$$$               endif
c$$$               

c$$$               if(am_i_root())then
c$$$!                  if(j.eq.50)then
c$$$               print*,'nuc NK',Nk(1), Nkd(1),Nk(2)
c$$$               print*,'nuc MK',Mk(1,1),Mk(1,2),Gc(srtso4)
c$$$!                   endif
c$$$               endif

               Mk(:,:)=Mkcond(:,:)
               Nk(:)=Nkcond(:)

               
               Gc(srtnh4)=Gcout(srtnh4)
               Gc(srtso4)=Gcout(srtso4)

               mpnum=1
               call aerodiag(mpnum,i,j,l)

c$$$               if(am_i_root())then
c$$$                  print*,'i',i,'j',j,'l',l
c$$$                  print*,'cond NK',Nk(1),Nkd(1)
c$$$                  print*,'cond MK',Mk(8,srtso4),Mkd(8,srtso4)
c$$$               endif
               

c$$$               if(am_i_root())then
c$$$!                  if(j.eq.50)then
c$$$               print*,'cond_nuc NK',Nk(1), Nkd(1),Nk(2)
c$$$               print*,'cond_nuc MK',Mk(1,1),Mk(1,2),Gc(srtso4)
c$$$!                   endif
c$$$               endif

               Mk(:,:)=Mkout(:,:)
               Nk(:)=Nkout(:)

c$$$               if(am_i_root())then
c$$$                  print*,'i',i,'j',j,'l',l
c$$$                  print*,'cond/nuc out NK',Nk(1),Nkd(1)
c$$$                  print*,'cond/nuc out MK',Mk(8,srtso4),Mkd(8,srtso4)
c$$$               endif
               

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
               
               
                                !Coagulation
               call storenm()
               call multicoag(adt)

c$$$               if(am_i_root())then
c$$$!                  if(j.eq.50)then
c$$$               print*,'multicoag NK',Nk(1), Nkd(1),INIT_Nk(1)
c$$$               print*,'multicoag MK',Mk(1,1),Mk(1,2),Gc(srtso4)
c$$$!                   endif
c$$$               endif

               mpnum=2
               call aerodiag(mpnum,i,j,l)

c$$$               if(am_i_root())then
c$$$                  print*,'i',i,'j',j,'l',l
c$$$                  print*,'coag NK',Nk(1),Nkd(1)
c$$$                  print*,'coag MK',Mk(8,srtso4),Mkd(8,srtso4),
c$$$     &                 aerod(i,j,l,idtso4+7,mpnum)
c$$$               endif
               
C     Do water eqm at appropriate times
               call eznh3eqm(Gc,Mk)
               call ezwatereqm(Mk)
                  
!            enddo
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
               

               if(am_i_root())then

                  if (abs(tot_n_2-tot_n_1)/tot_n_1.gt.1.0D-4)then
                     print*,'Nitrogen not conserved in aerophys'
                     print*,'i',i,'j',j,'l',l
                     print*,'Init,Init1,Intm,Final',tot_n_i,tot_n_1,
     *                    tot_n_1a,tot_n_2
                  endif
                  
                  if (abs(tot_s_2-tot_s_1)/tot_s_1.gt.1.0D-4)then !TOMAS - increase from 1.0D-4 
                     print*,'Sulfur not conserved in aerophys'
                     print*,'i',i,'j',j,'l',l
                     print*,'Init,Init1,Intm,Final',tot_s_1,
     *                    tot_s_1b,tot_s_2
c$$$  print*,'Intermediate',tot_s_1a
c$$$  print*,'Final',tot_s_2
c$$$  print*,'Gen',H2SO4rate_o*adt*32.d0/98.d0
c     STOP
                  endif
               endif
              
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
                                ! regular bins

!TOMAS - TRM should not be updated here.. 
c$$$!!!!!  diagnostics !!!!!!
c$$$
c$$$               do n=1,NBINS
c$$$                  tracnum=IDTNUMD-1+n
c$$$                  if (Nk(n) .ge. TRM(i,j,l,tracnum)) then
c$$$                     TRM(i,j,l,tracnum)=Nk(n)
c$$$                  else
c$$$                     frac=Nk(n)/TRM(i,j,l,tracnum)
c$$$                     call scalemom(i,j,l,tracnum,frac)
c$$$                  endif
c$$$                  do jc=1,icomp-idiag
c$$$                     tracnum=IDTSO4-1+n+ibins*(jc-1)
c$$$                     if (Mk(n,jc) .ge. TRM(i,j,l,tracnum)) then
c$$$                        TRM(i,j,l,tracnum)=Mk(n,jc)
c$$$                     else
c$$$                        frac=Mk(n,jc)/TRM(i,j,l,tracnum)
c$$$                        call scalemom(i,j,l,tracnum,frac)
c$$$                     endif
c$$$                  enddo
c$$$                  tracnum=IDTH2O-1+n
c$$$                  if (Mk(n,srth2o) .ge. TRM(i,j,l,tracnum)) then
c$$$                     TRM(i,j,l,tracnum)=Mk(n,srth2o)
c$$$                  else
c$$$                     frac=Mk(n,srth2o)/TRM(i,j,l,tracnum)
c$$$                     call scalemom(i,j,l,tracnum,frac)
c$$$                  endif         
c$$$               enddo
c$$$               
c$$$                                ! gas phase sulfuric acid
c$$$               if (Gc(srtso4) .ge. TRM(i,j,l,n_H2SO4)) then
c$$$                  TRM(i,j,l,n_H2SO4)=Gc(srtso4)
c$$$               else
c$$$                  frac=Gc(srtso4)/TRM(i,j,l,n_H2SO4)
c$$$                  call scalemom(i,j,l,n_H2SO4,frac)
c$$$               endif
c$$$                                               ! gas phase ammonia
c$$$               if (Gc(srtnh4) .ge. TRM(i,j,l,n_NH3)) then
c$$$                  TRM(i,j,l,n_NH3)=Gc(srtnh4)
c$$$               else
c$$$                  frac=Gc(srtnh4)/TRM(i,j,l,n_NH3)
c$$$                  call scalemom(i,j,l,n_NH3,frac)
c$$$               endif
c$$$                                ! aerosol ammonia
c$$$               tot_aam = 0.d0
c$$$               do n=1,NBINS
c$$$                  tot_aam = tot_aam + Mk(n,srtnh4)
c$$$               enddo
c$$$               if (tot_aam .ge. TRM(i,j,l,n_NH4)) then
c$$$                  TRM(i,j,l,n_NH4)=tot_aam
c$$$               else
c$$$                  frac=tot_aam/TRM(i,j,l,n_NH4)
c$$$                  call scalemom(i,j,l,n_NH4,frac)
c$$$               endif
c$$$                                ! SOAgas
c$$$               TRM(i,j,l,n_SOAgas)=INIT_SOA-SOArate*adt
c$$$              
c$$$               frac=TRM(i,j,l,n_SOAgas)/INIT_SOA
c$$$               call scalemom(i,j,l,n_SOAgas,frac)

c$$$C     Check for negative tracer problems
c$$$               do n=1,ibins
c$$$
c$$$                  if (Nk(n) .lt. 0.0) then
c$$$                     if (abs(Nk(n)) .gt. 1.e-5) then
c$$$                                !serious problem - report error
c$$$                    write(*,*) 'ERROR: Tracer ',n,' < 0 in box ', i,j,l
c$$$                    call stop_model ('ERROR: Tracer in TOMAS',255)
c$$$                     else
c$$$                                !numerical problem - set to zero
c$$$                        NK(n)=0.0
c$$$                     endif
c$$$                  endif
c$$$
c$$$                  do jc=1,icomp-idiag
c$$$                  if (Mk(n,jc) .lt. 0.0) then
c$$$                     if (abs(Mk(n,jc)) .gt. 1.e-5) then
c$$$                                !serious problem - report error
c$$$                    write(*,*) 'ERROR: Tracer ',n,' < 0 in box ', i,j,l
c$$$                    call stop_model ('ERROR: Tracer in TOMAS',255)
c$$$                     else
c$$$                                !numerical problem - set to zero
c$$$                        Mk(n,jc)=0.0
c$$$                     endif
c$$$                  endif
c$$$                  enddo
c$$$               enddo

               do n=1,ibins       
!     Aerosol number             
                  tracnum=IDTNUMD-1+n 
                  tr3Dsource(i,j,l,nOther,tracnum)=
     &                 (NK(N)-INIT_NK(N))/adt
                  
                  do np=1,ptype
                     taijs(i,j,ijts_TOMAS(np,tracnum)) 
     &                    =taijs(i,j,ijts_TOMAS(np,tracnum))
     &                    +AEROD(i,j,l,tracnum,np) ! /adt
                     if (itcon_TOMAS(np,tracnum).gt.0) 
     &                    call inc_diagtcb(i,j,AEROD(i,j,l,tracnum,np),
     &                    itcon_TOMAS(np,tracnum),tracnum)
                  enddo

                  do jc=1,icomp-idiag
                     tracnum=IDTSO4-1+n+ibins*(jc-1)
                     tr3Dsource(i,j,l,nOther,tracnum)=
     &                    (MK(n,jc)-INIT_Mk(n,jc))/adt

                  do np=1,ptype
                     taijs(i,j,ijts_TOMAS(np,tracnum)) 
     &                    =taijs(i,j,ijts_TOMAS(np,tracnum))
     &                    +AEROD(i,j,l,tracnum,np) ! /adt

                     if (itcon_TOMAS(np,tracnum).gt.0) 
     &                    call inc_diagtcb(i,j,AEROD(i,j,l,tracnum,np),
     &                    itcon_TOMAS(np,tracnum),tracnum)
                  enddo

                  enddo  
c$$$                  
                  tracnum=IDTH2O-1+n 
                  tr3Dsource(i,j,l,nOther,tracnum)=
     &                 (MK(N,SRTH2O)-INIT_MK(N,SRTH2O))/adt
               enddo
               

               tr3Dsource(i,j,l,nChemistry,n_H2SO4) =
     *              (Gc(srtSO4)-INIT_H2SO4)/adt
               
               tr3Dsource(i,j,l,nChemistry,n_NH3)=
     *              (Gc(srtNH4)-INIT_NH3)/adt

                                ! aerosol ammonia
               tot_aam = 0.d0
               do n=1,NBINS
                  tot_aam = tot_aam + Mk(n,srtnh4)
               enddo
               
               tr3Dsource(i,j,l,nChemistry,n_NH4)=
     *              (tot_aam-INIT_NH4)/adt

               tr3Dsource(i,j,l,nChemistry,n_SOAgas)=
     *              -SOArate
               
               AEROD(i,j,l,:,7)=0.0 !for aeroupdate

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


c$$$C     **************************************************
c$$$C     *  getdp                                         *
c$$$C     **************************************************
c$$$
c$$$C     WRITTEN BY Peter Adams
c$$$
c$$$C     This function calculates the average diameter of aerosol
c$$$C     particles in a given GCM grid cell and size bin.
c$$$
c$$$      subroutine dep_getdp(tm,getdp,size_density)                                            
c$$$!      USE TOMAS_AEROSOL
c$$$      USE TRACER_COM, only : nbins,IDTSO4,IDTNA,IDTECIL,
c$$$     &     IDTECOB,IDTOCIL,IDTOCOB,IDTDUST,IDTH2O,
c$$$     &     IDTNUMD,ntm,xk
c$$$      USE CONSTANT,   only : pi 
c$$$
c$$$      IMPLICIT NONE
c$$$
c$$$C-----ARGUMENT DECLARATIONS------------------------------------------
c$$$
c$$$      integer i,j,l  !coordinate of GCM grid cell
c$$$      integer n      !tracer index
c$$$
c$$$C-----VARIABLE DECLARATIONS------------------------------------------
c$$$
c$$$      integer k     !size bin index
c$$$      real density                 !density (kg/m3) of current size bin
c$$$      real mso4, mh2o, mno3, mnh4  !mass of each component (kg/grid box)
c$$$      real mecil,mecob,mocil,mocob
c$$$      real mdust,mtot,mnacl             
c$$$      real*8 mp          !particle mass (kg)
c$$$      real*8 mu          !air viscosity (kg/m s)
c$$$      real aerodens
c$$$!      external aerodens
c$$$!      REAL*8, INTENT(IN) :: TEMPK
c$$$      REAL*8, intent(in), DIMENSION(ntm) :: TM
c$$$      real,intent(out),DIMENSION(nbins) :: getdp,size_density
c$$$C-----VARIABLE COMMENTS----------------------------------------------
c$$$
c$$$C-----ADJUSTABLE PARAMETERS------------------------------------------
c$$$      real*8 Neps  !a small number of particles (#/box)
c$$$      parameter (Neps=1.d-20)
c$$$
c$$$C-----CODE-----------------------------------------------------------
c$$$
c$$$!     Compute particle diameter for each bins - YUNHA LEE 
c$$$      do k=1,nbins
c$$$              
c$$$         if (TM(IDTNUMD-1+k) .eq. 0.0) then
c$$$            if (TM(IDTSO4-1+k) .gt. 1.) then
c$$$               print*, 'ERROR in getdp - # = but mass > 0'
c$$$               print*, 'bin=',k
c$$$               print*, 'TM(#)=',TM(IDTNUMD-1+k)
c$$$               print*, 'TM(SO4)=',TM(IDTSO4-1+k)
c$$$               print*, 'TM(NACL)=',TM(IDTNA-1+k)
c$$$               print*, 'TM(OCIL)=',TM(IDTOCIL-1+k)
c$$$               call stop_model('ERROR IN getdp',255)
c$$$            else
c$$$!Not allow due to intent(in)?
c$$$!               print*, 'zero for mso4 and mh2o in getdp'
c$$$            endif
c$$$         endif
c$$$
c$$$         mso4=TM(IDTSO4-1+k) 
c$$$         mnacl=TM(IDTNA-1+k)
c$$$         mno3=0.e0
c$$$         if ((mso4+mno3) .lt. 1.e-8) mso4=1.e-8
c$$$         mnh4=0.1875*mso4  !assume ammonium bisulfate
c$$$         mecob=TM(IDTECOB-1+k)
c$$$         mecil=TM(IDTECIL-1+k)
c$$$         mocil=TM(IDTOCIL-1+k)
c$$$         mocob=TM(IDTOCOB-1+k)
c$$$         mdust=TM(IDTDUST-1+k)          
c$$$         mh2o=TM(IDTH2O-1+k)   
c$$$
c$$$!in CLOUDS2.f - some tracers goes negative..
c$$$!So, set to zero to prevent a problem in density calculation
c$$$         if(mnacl.lt.0) mnacl=0.
c$$$         if(mecob.lt.0) mecob=0.
c$$$         if(mecil.lt.0) mecil=0.
c$$$         if(mocob.lt.0) mocob=0.
c$$$         if(mocil.lt.0) mocil=0.
c$$$         if(mdust.lt.0) mdust=0.
c$$$         if(mh2o.lt.0) mh2o=0.
c$$$
c$$$         mtot= 1.1875*TM(IDTSO4-1+k)+mnacl+mecil+mecob+
c$$$     *        mocil+mocob+mdust+mh2o
c$$$
c$$$
c$$$
c$$$         density=aerodens(mso4, mno3,mnh4 !mno3 taken off!
c$$$     *        ,mnacl,mecil,mecob,mocil,mocob,mdust,mh2o) !assume bisulfate     
c$$$         size_density(k)=density      
c$$$         
c$$$!        print*,'density', K,density, size_density(k)
c$$$         
c$$$         if (TM(IDTNUMD-1+k) .gt. Neps.and.mtot.gt.0.) then
c$$$            mp=mtot/(TM(IDTNUMD-1+k))
c$$$         else
c$$$            mp=sqrt(xk(k+1)*xk(k))
c$$$            if(TM(IDTNUMD-1+k) .gt. Neps) 
c$$$     &           print*,'Warning in getdp:#>Neps but mtot=0',
c$$$     &           k,mtot,TM(IDTNUMD-1+k)
c$$$         endif
c$$$         
c$$$!     fix unrealistically large mp for low aerosol conc.
c$$$         if (mp .gt. 1.d3*xk(NBINS+1)) then
c$$$            
c$$$            if ((TM(IDTNUMD-1+k) .lt. 1.d5) .and. !negligible amount of aerosol - fudge mp
c$$$     &           (TM(IDTSO4-1+k) .lt. 3.)) then
c$$$               mp=sqrt(xk(k+1)*xk(k))
c$$$            else
c$$$               if (TM(IDTNUMD-1+k) .gt. 1.d12) then
c$$$!MODELE-TOMAS: during CONDSE, TM(H2O) is so large that causes too big mp. 
c$$$!MODELE-TOMAS: So, if dry mass is less than the max size boundary, just take the max mp. 
c$$$                  if((mtot-TM(IDTH2O-1+k)).lt.
c$$$     &                 1.d1*xk(nbins+1)*TM(IDTNUMD-1+k))then
c$$$                     print*,'Fudge mp in getdp: large mp by AH2O'
c$$$                     mp=sqrt(xk(nbins+1)*xk(nbins))                     
c$$$                  else
c$$$                  print*,'ERROR in getdp: mp too large'
c$$$                  print*, 'bin=',k
c$$$                  print*, 'TM(#)=', TM(IDTNUMD-1+k)
c$$$                  print*, 'TM(SO4)=', mso4, mh2o
c$$$                  print*, 'TM(NACL)=', mnacl, mdust
c$$$                  print*, 'TM(OC)=',mocob,mocil
c$$$                  print*, 'TM(EC)=',mecob,mecil
c$$$                  call stop_model('mp too large getdp',255)
c$$$                  endif
c$$$               endif
c$$$            endif
c$$$         endif
c$$$!         print*,'getdp',k, mp,size_density(k)
c$$$         getdp(k)=(6.d0*mp/(pi*size_density(k)))**(1.d0/3.d0)
c$$$c$$$         if(getdp(k).le.0.or.size_density(k).le.0.) then
c$$$c$$$            print*,'Dp=0 in getdp',k,getdp(k),size_density(k),mtot
c$$$c$$$            getdp(k)=
c$$$c$$$            CALL stop_model('Dp=0 in getdp',255)
c$$$c$$$         endif
c$$$      enddo
c$$$
c$$$! Finish computing particle diameter \
c$$$
c$$$      RETURN
c$$$      END subroutine dep_getdp
c$$$


C     **************************************************
C     *  drydep_getdp                                  *
C     **************************************************

C     WRITTEN BY Yunha Lee

C     This function calculates the average diameter of aerosol
C     particles in a given GCM grid cell and size bin.

      subroutine dep_getdp(i,j,l,getdp,size_density)                                            
!      USE TOMAS_AEROSOL
      USE TRACER_COM, only : nbins,IDTSO4,IDTNA,IDTECIL,
     &     IDTECOB,IDTOCIL,IDTOCOB,IDTDUST,IDTH2O,
     &     IDTNUMD,ntm,xk,trm
      USE CONSTANT,   only : pi 

      IMPLICIT NONE

C-----ARGUMENT DECLARATIONS------------------------------------------

      integer i,j,l  !coordinate of GCM grid cell
      integer n      !tracer index

C-----VARIABLE DECLARATIONS------------------------------------------

      integer k     !size bin index
      real density                 !density (kg/m3) of current size bin
      real mso4, mh2o, mno3, mnh4  !mass of each component (kg/grid box)
      real mecil,mecob,mocil,mocob
      real mdust,mtot,mnacl             
      real*8 mp          !particle mass (kg)
      real*8 mu          !air viscosity (kg/m s)
      real aerodens
      real,intent(out),DIMENSION(nbins) :: getdp,size_density
C-----VARIABLE COMMENTS----------------------------------------------

C-----ADJUSTABLE PARAMETERS------------------------------------------
      real*8 Neps  !a small number of particles (#/box)
      parameter (Neps=1.d-20)

C-----CODE-----------------------------------------------------------

!     Compute particle diameter for each bins - YUNHA LEE 
      do k=1,nbins
              
         if (TRM(I,J,L,IDTNUMD-1+k) .eq. 0.0) then
            if (TRM(I,J,L,IDTSO4-1+k) .gt. 1.) then
               print*, 'ERROR in getdp - # = but mass > 0',i,j
               print*, 'bin=',k
               print*, 'TRM(#)=',TRM(I,J,L,IDTNUMD-1+k)
               print*, 'TRM(SO4)=',TRM(I,J,L,IDTSO4-1+k)
               print*, 'TRM(NACL)=',TRM(I,J,L,IDTNA-1+k)
               print*, 'TRM(OCIL)=',TRM(I,J,L,IDTOCIL-1+k)
               call stop_model('ERROR IN getdp',255)
            endif
         endif

         mso4=TRM(I,J,L,IDTSO4-1+k) 
         mnacl=TRM(I,J,L,IDTNA-1+k)
         mno3=0.e0
         if ((mso4+mno3) .lt. 1.e-8) mso4=1.e-8
         mnh4=0.1875*mso4  !assume ammonium bisulfate
         mecob=TRM(I,J,L,IDTECOB-1+k)
         mecil=TRM(I,J,L,IDTECIL-1+k)
         mocil=TRM(I,J,L,IDTOCIL-1+k)
         mocob=TRM(I,J,L,IDTOCOB-1+k)
         mdust=TRM(I,J,L,IDTDUST-1+k)          
         mh2o=TRM(I,J,L,IDTH2O-1+k)   

!in CLOUDS2.f - some tracers goes negative..
!So, set to zero to prevent a problem in density calculation
         if(mnacl.lt.0) mnacl=0.
         if(mecob.lt.0) mecob=0.
         if(mecil.lt.0) mecil=0.
         if(mocob.lt.0) mocob=0.
         if(mocil.lt.0) mocil=0.
         if(mdust.lt.0) mdust=0.
         if(mh2o.lt.0) mh2o=0.

         mtot= 1.1875*TRM(I,J,L,IDTSO4-1+k)+mnacl+mecil+mecob+
     *        mocil+mocob+mdust+mh2o

         density=aerodens(mso4, mno3,mnh4 !mno3 taken off!
     *        ,mnacl,mecil,mecob,mocil,mocob,mdust,mh2o) !assume bisulfate     
         size_density(k)=density      
         
!        print*,'density', K,density, size_density(k)
         
         if (TRM(I,J,L,IDTNUMD-1+k) .gt. Neps.and.mtot.gt.0.) then
            mp=mtot/(TRM(I,J,L,IDTNUMD-1+k))
         else
            mp=sqrt(xk(k+1)*xk(k))
            if(TRM(I,J,L,IDTNUMD-1+k) .gt. Neps) 
     &           print*,'Warning in getdp:#>Neps but mtot=0',
     &           k,mtot,TRM(I,J,L,IDTNUMD-1+k)
         endif
         
!     fix unrealistically large mp for low aerosol conc.
         if (mp .gt. 1.d3*xk(NBINS+1)) then
            
            if ((TRM(I,J,L,IDTNUMD-1+k) .lt. 1.d5) .and. !negligible amount of aerosol - fudge mp
     &           (TRM(I,J,L,IDTSO4-1+k) .lt. 3.)) then
               mp=sqrt(xk(k+1)*xk(k))
            else
               if (TRM(I,J,L,IDTNUMD-1+k) .gt. 1.d12) then
!MODELE-TOMAS: during CONDSE, TM(H2O) is so large that causes too big mp. 
!MODELE-TOMAS: So, if dry mass is less than the max size boundary, just take the max mp. 
                  if((mtot-TRM(I,J,L,IDTH2O-1+k)).lt.
     &                 1.d1*xk(nbins+1)*TRM(I,J,L,IDTNUMD-1+k))then
                     print*,'Fudge mp in getdp: large mp by AH2O'
                     mp=sqrt(xk(nbins+1)*xk(nbins))                     
                  else
                  print*,'ERROR in getdp: mp too large'
                  print*, 'bin=',k
                  print*, 'TM(#)=', TRM(I,J,L,IDTNUMD-1+k)
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

! Finish computing particle diameter \

      RETURN
      END subroutine dep_getdp


C	readfraction.f
	subroutine readfraction(infile,fraction)

        implicit none

	character*80 infile
	integer innum, ii, jj, kk

	real*8,intent(out),dimension(101,101,101):: fraction
	parameter (innum=580)
 1	format(f6.5)
	open(unit=innum,file=infile,FORM='FORMATTED',STATUS='OLD')
	do ii=1,101
		do jj=1,101
			do kk=1,101
			read(innum,1) fraction(kk,jj,ii)
			if (fraction(kk,jj,ii).gt.1.) fraction(kk,jj,ii)=0.
			enddo
		enddo
	enddo
!        print*,'fraction last',fraction(101,101,101)
	close(innum)
	return 
	end subroutine readfraction
			 

C	readbinact.f
	subroutine readbinact(infile,binact)


        USE TOMAS_AEROSOL, ONLY : ibins

        implicit none

	character*80 infile
	integer innum, ii, jj, kk
	integer,intent(out),dimension(101,101,101):: binact
	parameter (innum=590)
 1	format(I2)
	open(unit=innum,file=infile,FORM='FORMATTED',STATUS='OLD')
	do ii=1,101
		do jj=1,101
			do kk=1,101
			read(innum,1) binact(kk,jj,ii)
			if (binact(kk,jj,ii).eq.0) binact(kk,jj,ii)=ibins+1
			enddo
		enddo
	enddo
!        print*,'binact last',binact(101,101,101)
	close(innum)
	return
	end subroutine readbinact
			 

	subroutine getfraction(tr_conv,tm,fraction)
Ckpc  change binact values
        USE TOMAS_AEROSOL, ONLY : binact02,binact10,
     &       fraction02,fraction10 !,GCM_fraction
        USE TRACER_COM, only : nbins,ntm,IDTECIL,
     &       IDTOCIL,IDTOCOB,IDTSO4,IDTNA,IDTDUST,
     &       IDTECOB

        implicit none

	real mecil, mocil, mocob, mso4, mnacl,mdust, mtot
	real xocil, xso4, xnacl
	integer iso4, inacl, iocil,k
	integer getbinact
	real*8,intent(out),dimension(nbins) ::  fraction
        REAL*8,  INTENT(IN), DIMENSION(ntm) :: TM
        LOGICAL TR_CONV

        do k=1, nbins
        mecil=TM(IDTECIL-1+k)
	mocil=TM(IDTOCIL-1+k)
	mocob=TM(IDTOCOB-1+k)
	mso4=TM(IDTSO4-1+k)*1.2  !account for ammonium sulfate
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
     &          ,mocil,mdust
           call stop_model('wrong getfraction',255)
        endif

        if (tr_conv)then !convective clouds
           getbinact=binact10(iso4,inacl,iocil)
           
           if(binact10(iso4,inacl,iocil).lt.0)then
              print*,'wrong binact',binact10(iso4,inacl,iocil)
              print*,'iso4',iso4,inacl,iocil
              print*,'mass',mso4,mnacl,mocil,mtot
           call stop_model('wrong binact',255)
           endif
           
           if (getbinact.gt.k) then
              fraction(k)=0.    !not activated
!              GCM_fraction(i,j,l,k)=0.0
           else if (getbinact.eq.k) then
              
              fraction(k)=fraction10(iso4,inacl,iocil) !partly activated
!              GCM_fraction(i,j,l,k)=fraction(k)
           else
              fraction(k)=1.    !all sizebin activated
!              GCM_fraction(i,j,l,k)=1.
           endif
           if(getbinact.le.2)
     &          print*,'CONV CLD',getbinact,k,iso4,inacl,iocil
        else                    !large-scale
           getbinact=binact02(iso4,inacl,iocil)

           if(binact02(iso4,inacl,iocil).lt.0)then
              print*,'wrong binact',binact02(iso4,inacl,iocil)
              print*,'iso4',iso4,inacl,iocil
              print*,'mass',mso4,mnacl,mocil,mtot
           endif

           if (getbinact.gt.k) then
              fraction(k)=0.    !not activated
!              GCM_fraction(i,j,l,k)=0.0
           else if (getbinact.eq.k) then
              
              fraction(k)=fraction02(iso4,inacl,iocil) !partly activated
!              GCM_fraction(i,j,l,k)=fraction(k)
           else
              fraction(k)=1.    !all sizebin activated
!              GCM_fraction(i,j,l,k)=1.
         
           if(getbinact.le.4)
     &          print*,'STRAT CLD',getbinact,k,iso4,inacl,iocil

           endif
        endif
        enddo                   !k
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

      SUBROUTINE aqoxid(moxid,tr_conv,Nko,Mko,Nkf,Mkf)

C-----INCLUDE FILES-----------------------------------------------------

      USE TOMAS_AEROSOL
      USE TRACER_COM, only : ntm, IDTECIL,
     &       IDTOCIL,IDTOCOB,IDTSO4,IDTNA,IDTDUST,
     &       IDTECOB,IDTH2O,xk,nbins

      IMPLICIT NONE
C-----ARGUMENT DECLARATIONS---------------------------------------------

      real*8 moxid !mass of new sulfate from in-cloud oxid.
      real*8, dimension(ibins) :: fraction    !fraction activated for every sizebin

c      real*8 tracc(NTM)

C-----VARIABLE DECLARATIONS---------------------------------------------

      real*8 Nact, Mact  !#/mass of activated particles
      real*8 mpo   !initial particle mass (kg)
      real*8 mpw   !initial particle wet mass (kg)
      real*8 aqtau(ibins)
      integer k,mpnum,n,tracnum
      real*8 Nko(ibins), Mko(ibins, icomp) !input to cond routine
      real*8 Nkf(ibins), Mkf(ibins, icomp) !output from cond routine
      real*8 tdt      !the value 2/3
      real*8,parameter :: eps=1.d-40
      integer jc,j
      real*8 frac      
      real*8 WR                ! wet ratio = total mass/ dry mass (win, 5/15/06)
      real*8 mox(ibins) !mass of new sulfate per particle in each bin
      real*8 tot_aam ! total aerosol ammonia
      real*8 TM(ntm) ! total aerosol ammonia

      LOGICAL TR_CONV

C-----CODE--------------------------------------------------------------

              
      if (moxid.eq.0.d0) return

      tdt=2.d0/3.d0

      TM(:)=0.0
!only mass needed for getfraction
      do n=1,NBINS
         do jc=1,icomp-idiag
            tracnum=IDTSO4-1+n+ibins*(jc-1)
            TM(tracnum)=Mko(n,jc)               
         enddo
         tracnum=IDTH2O-1+n
         TM(tracnum)=Mko(n,srth2o)
      enddo

      if (tr_conv) then
         CALL getfraction (.true.,TM,FRACTION) !1% supersaturation assumption
      else
         CALL getfraction(.false.,TM,FRACTION) !0.2% supersaturation assumption
      endif

      Nact=0.0
      Mact=0.0
      do k=1,ibins
         Nact=Nact+Nko(k)*fraction(k)
         do jc=1,icomp-idiag
            Mact=Mact+Mko(k,jc)*fraction(k)
         enddo
      enddo

      if ((Mact+moxid)/(Nact+eps) .gt. xk(ibins)) then !YHL- I change xk(ibins-1) to xk(ibins)
!            if (TAU .gt. 8350.) then
               write(*,*) 'ERROR in aqoxidcc: Ave size out of bounds'
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
               mpo = mpo+Mko(k,jc)  !accumulate dry mass
            enddo
            do jc=1,icomp
               mpw = mpw+Mko(k,jc)  ! have wet mass include amso4
            enddo
            WR = mpw/mpo  !WR = wet ratio = total mass/dry mass
            if (Nko(k) .gt. 0.d0) then
               mpw=mpw/Nko(k)
               aqtau(k)=1.5d0*((mpw+mox(k)*WR)**tdt-mpw**tdt)  !added WR to moxid term (win, 5/15/06)
            else
               !nothing in this bin - set tau to zero
               aqtau(k)=0.0
               mox(k)=0.d0
            endif
         endif
      enddo

      call tmcond(aqtau,xk,Mko,Nko,Mkf,Nkf,srtso4,mox)

C Do water eqm at appropriate times
!      call eznh3eqm(Gc,Mkf)
      call ezwatereqm(Mkf)      

      do k=1,ibins  
         if(Nkf(k).lt.0)then
            print*,'ERROR:Nk < 0',k, Nkf(k)
         endif
         do j=1,icomp
            if(Mkf(k,j).lt.0)then
               print*,'ERROR:Mk < 0',k,j,Mkf(k,j)
            endif
         enddo
      enddo

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


!      print*,'inodens init',so4temp,no3temp,nh4temp,h2otemp,nacltemp
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
!      call nanstop(xtot,153,0,0)
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
c      call nanstop(d,173,0,0)
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
!      print*,'inodens end',mso4,mno3,so4temp,no3temp
!      call stop_model('density',13)
      mso4=so4temp
      mno3=no3temp
      mnh4=nh4temp
      mnacl=nacltemp
      mh2o=h2otemp

!      print*,'inodens end2',mso4,mno3,mnh4,mh2o,mnacl

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

C-----ARGUMENT DECLARATIONS---------------------------------------------


      real*8 Mke(ibins,icomp)

C-----VARIABLE DECLARATIONS---------------------------------------------

      integer k
      real*8 so4mass, naclmass, ocilmass
      real*8 wrso4, wrnacl, wrocil
      real*8 rhe

      real*8 waterso4, waternacl, waterocil
!      external waterso4, waternacl, waterocil

C     VARIABLE COMMENTS...

C     This version of the routine works for sulfate and sea salt
C     particles.  They are assumed to be externally mixed and their
C     associated water is added up to get total aerosol water.
C     wr is the ratio of wet mass to dry mass of a particle.  Instead
C     of calling a thermodynamic equilibrium code, this routine uses a
C     simple curve fits to estimate wr based on the current humidity.
C     The curve fit is based on ISORROPIA results for ammonium bisulfate
C     at 273 K and sea salt at 273 K.

C-----ADJUSTABLE PARAMETERS---------------------------------------------

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

C-----ARGUMENT DECLARATIONS---------------------------------------------

      real*8 Gce(icomp)
      real*8 Mke(ibins,icomp)

C-----VARIABLE DECLARATIONS---------------------------------------------

      integer k
      real*8 tot_nh3  !total kmoles of ammonia
      real*8 tot_so4  !total kmoles of so4
      real*8 sfrac    !fraction of sulfate that is in that bin

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

C-----INPUTS------------------------------------------------------------

C-----OUTPUTS-----------------------------------------------------------

      real*8 FUNCTION waterso4(rhe)

      IMPLICIT NONE

C-----ARGUMENT DECLARATIONS---------------------------------------------

      real*8 rhe   !relative humidity (0-100 scale)

C-----VARIABLE DECLARATIONS---------------------------------------------

C     VARIABLE COMMENTS...

C     waterso4 is the ratio of wet mass to dry mass of a particle.  Instead
C     of calling a thermodynamic equilibrium code, this routine uses a
C     simple curve fit to estimate wr based on the current humidity.
C     The curve fit is based on ISORROPIA results for ammonium bisulfate
C     at 273 K.

C-----ADJUSTABLE PARAMETERS---------------------------------------------

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

C-----INPUTS------------------------------------------------------------

C-----OUTPUTS-----------------------------------------------------------

      real*8 FUNCTION waternacl(rhe)

      IMPLICIT NONE


C-----ARGUMENT DECLARATIONS---------------------------------------------

      real*8 rhe   !relative humidity (0-100 scale)

C-----VARIABLE DECLARATIONS---------------------------------------------

C     VARIABLE COMMENTS...

C     waternacl is the ratio of wet mass to dry mass of a particle.  Instead
C     of calling a thermodynamic equilibrium code, this routine uses a
C     simple curve fit to estimate waternacl based on the current humidity.
C     The curve fit is based on ISORROPIA results for sodium sulfate
C     at 273 K.

C-----ADJUSTABLE PARAMETERS---------------------------------------------

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

C-----INPUTS------------------------------------------------------------

C-----OUTPUTS-----------------------------------------------------------

      real*8 FUNCTION waterocil(rhe)

      IMPLICIT NONE

C-----ARGUMENT DECLARATIONS---------------------------------------------

      real*8 rhe   !relative humidity (0-100 scale)
      real*8 a, b, c, d, e, f, prefactor, activcoef
      parameter(a=1.0034, b=0.1614, c=1.1693,d=-3.1,
     & e=6.0)

C-----VARIABLE DECLARATIONS---------------------------------------------

C     VARIABLE COMMENTS...

C     waterocil is the ratio of wet mass to dry mass of a particle.  Instead
C     of calling a thermodynamic equilibrium code, this routine uses a
C     simple curve fit to estimate waterocil based on the current humidity.
C     The curve fit is based on observations of Dick et al. JGR D1 1471-1479

C-----ADJUSTABLE PARAMETERS---------------------------------------------

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


      USE DOMAIN_DECOMP_ATM, only : GRID, GET, write_parallel,
     &     am_i_root
      USE TOMAS_AEROSOL 
      USE GEOM, only: imaxj
      USE TRACER_COM, only : IDTSO4, IDTNA, IDTOCIL,IDTH2O,NBINS
     &     ,trm,IDTECOB,IDTECIL,IDTOCOB,IDTDUST,IDTNUMD,TRNAME
     *     ,ntm,ntm_TOMAS

      USE TRDIAG_COM, only : taijs=>taijs_loc !,taijls=>taijls_loc

      USE MODEL_COM, only : im,jm,lm     ! dimensions
     $                     ,q            ! saturated pressure
     $                     ,t
     $                     ,dtsrc
      USE CONSTANT,   only:  lhe
      USE DYNAMICS,   only: pmid,pk ! midpoint pressure in hPa (mb)

      IMPLICIT NONE

      INTEGER J_0, J_1, I_0, I_1

      integer i,j,l,k,n,jc,mpnum  !counters
      real*8 qsat         !used in RH calculation
      integer tracnum
      
      real*8 frac, Nkout(iBINS),Mkout(iBINS,icomp),Gcout(icomp-1)
      
      real*8 rhe


      real*8 waterso4, waternacl, waterocil
!      REAL(8):: TK
!      external waterso4, waternacl, waterocil

C     VARIABLE COMMENTS...


C-----CODE-----------------------------------------------------------    
C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT=J_0,       J_STOP=J_1)
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

c$$$      !Do water eqm at appropriate times

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
c$$$            write(*,*) 'TRM<=0 in aeroupdate ',
c$$$     &           trname(IDTSO4+n-1),trm(i,j,l,IDTSO4+n-1)
            if (abs(TRM(i,j,l,IDTSO4+n-1)) .gt. 1.e-10) then
               !serious problem - report error
               write(*,*) 'ERROR: Tracer ',trname(IDTSO4+n-1),
     &              trm(i,j,l,IDTSO4+n-1)
               write(*,*) ' < 0 in box ', i,j,l
               call stop_model('TRM<0 in aeroupdate',255)
            else
               !numerical problem - set to zero
               TRM(i,j,l,IDTSO4+n-1)=0.0!1.d-42 !5??
            endif
         endif
      enddo

            enddo
         enddo
      enddo
      
      RETURN
      END SUBROUTINE aeroupdate



c$$$C     **************************************************
c$$$C     *  GCM_ezwatereqm                                    *
c$$$C     **************************************************
c$$$
c$$$C     WRITTEN BY Peter Adams, March 2000
c$$$C     MODIFIED BY Yunha Lee, March 2011 - to use this in GCM subroutine 
c$$$
c$$$
c$$$C     This routine uses the current RH to calculate how much water is 
c$$$C     in equilibrium with the aerosol.  Aerosol water concentrations 
c$$$C     are assumed to be in equilibrium at all times and the array of 
c$$$C     concentrations is updated accordingly.
c$$$
c$$$C     This version of the routine works for sulfate and sea salt
c$$$C     particles.  They are assumed to be externally mixed and their
c$$$C     associated water is added up to get total aerosol water.
c$$$C     wr is the ratio of wet mass to dry mass of a particle.  Instead
c$$$C     of calling a thermodynamic equilibrium code, this routine uses a
c$$$C     simple curve fits to estimate wr based on the current humidity.
c$$$C     The curve fit is based on ISORROPIA results for ammonium bisulfate
c$$$C     at 273 K and sea salt at 273 K.
c$$$
c$$$
c$$$      SUBROUTINE aeroupdate
c$$$
c$$$
c$$$      USE DOMAIN_DECOMP_ATM, only : GRID, GET, write_parallel,
c$$$     &     am_i_root
c$$$      USE TOMAS_AEROSOL 
c$$$      USE GEOM, only: imaxj
c$$$      USE TRACER_COM, only : IDTSO4, IDTNA, IDTOCIL,IDTH2O,NBINS
c$$$     &     ,trm,IDTECOB,IDTECIL,IDTOCOB,IDTDUST,IDTNUMD,TRNAME
c$$$     *     ,ntm,xk
c$$$
c$$$      USE TRDIAG_COM, only : taijs=>taijs_loc,taijls=>taijls_loc
c$$$
c$$$      USE MODEL_COM, only : im,jm,lm     ! dimensions
c$$$     $                     ,q            ! saturated pressure
c$$$
c$$$      USE CONSTANT,   only:  lhe
c$$$      USE DYNAMICS,   only: pmid ! midpoint pressure in hPa (mb)
c$$$
c$$$      IMPLICIT NONE
c$$$
c$$$      INTEGER J_0, J_1, I_0, I_1
c$$$
c$$$      integer i,j,l,k,n,jc  !counters
c$$$
c$$$      real Dp(nbins)            !particle diameter (m)
c$$$      real density                 !density (kg/m3) of current size bin
c$$$      real mso4, mh2o, mno3, mnh4  !mass of each component (kg/grid box)
c$$$      real mecil,mecob,mocil,mocob
c$$$      real mdust,mtot,mnacl             
c$$$      double precision mp          !particle mass (kg)
c$$$
c$$$      real aerodens
c$$$      external aerodens
c$$$
c$$$
c$$$C-----CODE-----------------------------------------------------------    
c$$$C****
c$$$C**** Extract useful local domain parameters from "grid"
c$$$C****
c$$$      CALL GET(grid, J_STRT=J_0,       J_STOP=J_1)
c$$$      I_0 = grid%I_STRT
c$$$      I_1 = grid%I_STOP
c$$$
c$$$
c$$$C     Loop over all grid cells
c$$$      DO L=1,LM                            
c$$$         DO J=J_0,J_1                          
c$$$            DO I=I_0,IMAXJ(J)
c$$$               DO K=1,NBINS
c$$$                  mso4=TRM(i,j,l,IDTSO4-1+k) 
c$$$                  mnacl=TRM(i,j,l,IDTNA-1+k)
c$$$                  mno3=0.e0
c$$$                  if ((mso4+mno3) .lt. 1.e-8) mso4=1.e-8
c$$$                  mnh4=0.1875*mso4 !assume ammonium bisulfate
c$$$                  mecob=TRM(i,j,l,IDTECOB-1+k)
c$$$                  mecil=TRM(i,j,l,IDTECIL-1+k)
c$$$                  mocil=TRM(i,j,l,IDTOCIL-1+k)
c$$$                  mocob=TRM(i,j,l,IDTOCOB-1+k)
c$$$                  mdust=TRM(i,j,l,IDTDUST-1+k)          
c$$$                  mh2o=TRM(i,j,l,IDTH2O-1+k)   
c$$$                  
c$$$                  if ((mnacl) .lt. 0) mnacl=0.
c$$$                  if ((mecob) .lt. 0) mecob=0.
c$$$                  if ((mecil) .lt. 0) mecil=0.
c$$$                  if ((mocob) .lt. 0) mocob=0.
c$$$                  if ((mocil) .lt. 0) mocil=0.
c$$$                  if ((mdust) .lt. 0) mdust=0.
c$$$                  
c$$$                  
c$$$                  mtot= 1.1875*mso4+mnacl+mecil+mecob+
c$$$     *                 mocil+mocob+mdust+mh2o
c$$$                  
c$$$                  density=aerodens(mso4, mno3,mnh4 !mno3 taken off!
c$$$     *                 ,mnacl,mecil,mecob,mocil,mocob,mdust,mh2o) !assume bisulfate     
c$$$                  
c$$$                  
c$$$                  if (TRM(i,j,l,IDTNUMD-1+k) .gt.1.d-20) then ! arbituary number !
c$$$                     mp=mtot/(TRM(i,j,l,IDTNUMD-1+k))
c$$$                  else
c$$$                     mp=sqrt(xk(k+1)*xk(k))
c$$$                  endif
c$$$                  
c$$$!     fix unrealistically large mp for low aerosol conc.
c$$$             if (mp .gt. 1.d3*xk(NBINS+1)) then            
c$$$                if ((TRM(i,j,l,IDTNUMD-1+k) .lt. 1.d5) .and. !negligible amount of aerosol - fudge mp
c$$$     &               (TRM(i,j,l,IDTSO4-1+k) .lt. 3.)) then
c$$$                   mp=sqrt(xk(k+1)*xk(k))
c$$$                else
c$$$                   if (TRM(i,j,l,IDTNUMD-1+k) .gt. 1.d12) then
c$$$                      print*,'ERROR in aeroupdate: mp too large'
c$$$                      print*, 'bin=',k,'i,j,',i,j,l
c$$$                      print*, 'TRM(#)=', TRM(i,j,l,IDTNUMD-1+k)
c$$$                      print*, 'TRM(SO4)=', mso4, mh2o
c$$$                      print*, 'TRM(NACL)=', mnacl, mdust
c$$$                      print*, 'TRM(OC)=',mocob,mocil
c$$$                      print*, 'TRM(EC)=',mecob,mecil
c$$$                      call stop_model('mp too large aeroupdate',13)
c$$$                   endif
c$$$                endif
c$$$             endif   
c$$$          enddo
c$$$
c$$$C Check for negative tracer problems
c$$$      do n=1,ntm
c$$$         if (TRM(i,j,l,n) .lt. 0.0) then
c$$$            if (abs(TRM(i,j,l,n)) .gt. 1.e-10) then
c$$$               !serious problem - report error
c$$$               write(*,*) 'ERROR: Tracer ',trname(n),trm(i,j,l,n)
c$$$               write(*,*) ' < 0 in box ', i,j,l
c$$$               call stop_model('TRM<0 in aeroupdate',13)
c$$$            else
c$$$               !numerical problem - set to zero
c$$$               TRM(i,j,l,n)=1.e-20 !5??
c$$$            endif
c$$$         endif
c$$$      enddo
c$$$            enddo
c$$$         enddo
c$$$      enddo
c$$$      
c$$$      RETURN
c$$$      END SUBROUTINE aeroupdate


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
      USE TRACER_COM, only : IDTSO4,IDTNUMD

      implicit none 
      integer pt, i, j, l, jc, n
      real adt            !aerosol microphysics time step (seconds)

      integer k,kk,tracnum
C-----CODE--------------------------------------------------------------

      if(pt.ne.7)then 

         do n=1,ibins       
                       
!     Aerosol number
            tracnum=IDTNUMD-1+n
            AEROD(i,j,l,tracnum,pt)= 
     &           Nk(n)-Nkd(n)
!     Aerosol mass
            do jc=1,icomp-idiag
               tracnum=IDTSO4-1+n+ibins*(jc-1)
               AEROD(i,j,l,tracnum,pt)= 
!     &           AEROD(i,j,l,tracnum,pt) 
     &              Mk(n,jc)-Mkd(n,jc)
               
            enddo   
         enddo
         
      else

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


      subroutine alloc_tracer_TOMAS_com(grid)
!@SUM  To alllocate arrays whose sizes now need to be determined
!@+    at run-time
!@auth Susanne Bauer
!@ver  1.0
      use domain_decomp_atm, only : dist_grid, get
      use model_com, only     : lm
!      use tracer_com, only    : ntmAMP
      use TOMAS_aerosol

      IMPLICIT NONE

      type (dist_grid), intent(in) :: grid
      integer :: ier, J_1H, J_0H, I_1H, I_0H
      logical :: init = .false.

      if(init)return
      init=.true.
    
      call get( grid , J_STRT_HALO=J_0H, J_STOP_HALO=J_1H )
      I_0H=GRID%I_STRT_HALO
      I_1H=GRID%I_STOP_HALO 

! I,J,L
      allocate(  AQSO4oxid_mc(I_0H:I_1H,J_0H:J_1H,LM)   )
      allocate(  AQSO4oxid_ls(I_0H:I_1H,J_0H:J_1H,LM)   )
! other dimensions
!      allocate(  SOA_chem(I_0H:I_1H,J_0H:J_1H)  )
      allocate(  H2SO4_chem(I_0H:I_1H,J_0H:J_1H,LM)  )

      allocate( AEROD(I_0H:I_1H,J_0H:J_1H,LM,NTM,ptype) )
     
      return
      end subroutine alloc_tracer_TOMAS_com
