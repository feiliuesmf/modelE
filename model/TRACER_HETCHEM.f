#include "rundeck_opts.h"

#ifdef TRACERS_SPECIAL_Shindell
      SUBROUTINE HETCDUST
!
! Version 1.   (version 2 needs to be written... without integration over ndr)
!
!-----------------------------------------------------------------------
!   Computation of heterogeneous reaction rates on dust aerosol surfaces
!   Susanne E Bauer, 2003
!-----------------------------------------------------------------------

      USE MODEL_COM, only : im,jm,lm     ! dimensions (72, 46, 12)
     $                     ,jday         ! time in ITU
     $                     ,t            ! potential temperature (C)
     $                     ,q            ! saturatered pressure

      USE TRACER_COM                   !, only : rxts,trm
      USE CONSTANT,   only:  lhe       ! latent heat of evaporation at 0 C
      USE GEOM,       only:  bydxyp
      USE DYNAMICS,   only:  byam ,pmid,pk   ! midpoint pressure in hPa (mb)
c                                          and pk is t mess up factor
      USE CONSTANT,   only:  pi, avog, gasc

      IMPLICIT NONE
!-----------------------------------------------------------------------
!       ... Dummy arguments
!-----------------------------------------------------------------------
      real                  :: erf
      integer, parameter    :: ndtr = 8  ! # dust bins
      real*8                 :: rxtnox(im,jm,lm,ndtr,rhet) !rate for single dust tracers
      real*8                :: dusttx(im,jm,lm,ndtr) ! dust (kg/kg) 
      real                  :: dustnc(im,jm,lm,ndtr) ! dust number concentr.
!-----------------------------------------------------------------------
!       ... Look up variables
!----------------------------------------------------------------------- 
      integer, parameter :: klo = 1000
      integer ip,imd,np1,np2,nh1,nh2
      real klook,phelp,nmd
      real look_p, look_t,hx,px,hp1,hp2
!-----------------------------------------------------------------------
!       ... Local variables
!-----------------------------------------------------------------------
      integer :: i, j, k, nd, l ,ll, il
      integer, parameter :: ktoa = 300
! 1-SO2
c      real, parameter :: alph1  = 0.0001 !uptake coeff of Rossi EPFL (independent of humidity)
!      real, parameter :: alph(rhet)=(/0.001,0.001,0.003/) !uptake coeff for HNO3,N2O5,NO3
      real, parameter :: alph(rhet)=(/0.0001,0.001,0.003/) !uptake coeff for HNO3,N2O5,NO3
      real, parameter :: mQ(rhet)=(/0.063,0.108,0.062/)
      real, parameter :: xx    = 0.  !correction factor anisotropic movement
      real, parameter :: Bolz  = 1.3807e-23 !Boltzmann kg m2/s2 K molec.
      real, parameter :: Mgas  = 28.97 /1000. ! Molekular Gewicht Luft
      real, parameter :: Diaq  = 4.5e-10      ! m Molecul Diameter
C**** functions
      real*8 :: QSAT,RH,te,temp

      real :: Kn(rhet), Mdc(rhet), Kdj(rhet)
      real :: lamb(rhet), wrk(rhet),VSP(rhet)
      real :: lsig0,drada,dn,Roh

      logical, save  :: enteredb = .false.
      real, save, dimension(ktoa) :: rada
      real, save     :: lookS(11,klo,rhet),Rrange,md_look(klo)

!-----------------------------------------------------------------
!     Dust variables: Radius of dust particles in the 8 size bins
!-----------------------------------------------------------------

      real, parameter  :: Dradi(8)=(/0.15e-6,0.25e-6,0.4e-6,0.8e-6,
     $                               1.5e-6,2.5e-6,4.e-6,8.e-6/)
      real, parameter  :: rop(8) = (/ 2500.,  2500.,  2500.,  2500.
     $                              , 2600.,  2600.,  2600.,  2600./)
!-----------------------------------------------------------------
!    1000 Intervals for Radius = 0.01ym ->10ym
!-----------------------------------------------------------------
!      Integration of radius:  0.01 ym to 10 ym
       rada(1) = 0.01E-6    ! smallest radius
       drada   = 0.1E-6     ! delta radius


      if (.not. enteredb) then
      enteredb = .true.
      PRINT*, 'CALCULATING LOOK UP TABLE FOR HETEROGENEOUS CHEMISTY'
      DO i   = 2, ktoa
      rada(i) = rada(i-1) + drada
      END DO

c      enteredb = .true.

      lsig0 = LOG(2.)

!-----------------------------------------------------------------
!     HNO3 + DUSTM =>    Dust Aerosol Reaction
!-----------------------------------------------------------------

  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!                         LOOK UP TABLE
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       md_look(1) = 1.E-10       ! smallest md
       Rrange=5.E-8

       DO i   = 2, klo
       md_look(i) = md_look(i-1) +  Rrange 
       END DO

      DO il  = 1, rhet  ! no loop over rhet, only HNO3 uptake  
      wrk(il)=  (mQ(il) + Mgas) / mQ(il)

      DO ip  = 1, 11  !pressure from 1000 to 0 hPa 

      look_p=max(0.001,1.1-ip*0.1)         !atmosphere (minimum is 1 hPa)
      look_t=max(210.,288.*(look_p/1.)**((1.40-1)/1.40))
      look_p=look_p * 100000.              ! pressure in Pa

       Roh    = look_p / look_t / 287.
C Molecular diffusion coefficient for a trace gas in air [ m/s ]
       Mdc(il)  = 3. / (8.* Avog * Roh * (Diaq**2.))
       Mdc(il)  = Mdc(il) * SQRT( ((gasc*look_t*Mgas)/(2.*pi))*wrk(il))
C thermal velocity of a trace gas molecule [m/s2]
       VSP(il)  = SQRT((8. * Bolz * look_t)/(Pi * mQ(il)/ Avog))
C lamb  mean free pathway  [m]
       lamb(il)   = 3. * Mdc(il)/VSP(il)
C Loop over radius

      DO imd  = 1,klo
       lookS(ip,imd,il) = 0.
      DO k = 1,ktoa-1
C Knudsen Number 
       Kn(il)= lamb(il) / rada(k) ! Radius in [m]
C Mass Transfer Coefficient
       Kdj(il) =(4. * pi * rada(k) * Mdc(il)) 
     .      / (1. + Kn(il) * (xx + 4 *(1.- alph(il))/(3. *alph(il))))
C Number distribution
       dn=abs(erf(log( rada(k)/md_look(imd)) / lsig0
     .       /sqrt(2.))-erf(log(rada(k+1)/ 
     .       md_look(imd)) / lsig0 /sqrt(2.)))/2.
C Net removal rate [s-1]
       lookS(ip,imd,il)= lookS(ip,imd,il) + Kdj(il)  * dn    
      END DO              !radius loop
      END DO              !median diameter loop
      END DO              !pressure loop
      END DO              !reaction loop

      ENDIF

  
c--------------------------------------------------------------

#ifdef TRACERS_DUST      

      dusttx(:,:,:,:)= 0.
      
      DO l  = 1,lm   
      DO j  = 1,jm   
      
      dusttx(:,j,l,5)= trm(:,j,l,n_clay) * byam(l,:,j)* bydxyp (j)
      dusttx(:,j,l,6)= trm(:,j,l,n_silt1)* byam(l,:,j)* bydxyp (j)
      dusttx(:,j,l,7)= trm(:,j,l,n_silt2)* byam(l,:,j)* bydxyp (j)
      dusttx(:,j,l,8)= trm(:,j,l,n_silt3)* byam(l,:,j)* bydxyp (j)
      
      enddo
      enddo

#else
c Read dust fields from files
      CALL READDUST(dusttx)
      print*, 'OFF-LINE DUST TRACER'
c      print*,'MAX OFF LINE:  ',maxval(dusttx(:,:,:,:))
c      print*,dusttx(30,23,3,:)
#endif
c--------------------------------------------------------------
c--------------------------------------------------------------

c INTERPOLATION FROM LOOK UP TABLES
 
C Net removal rates [s-1]  
        krate(:,:,:,:,:) = 0.
        rxtnox(:,:,:,:,:)=0
      DO il = 1,rhet    ! Loop over het reactions
#ifdef TRACERS_DUST
      DO nd = 5,ndtr    ! Loop over dust tracers
#else
      DO nd = 1,ndtr    ! Loop over dust tracers
#endif
      DO l  = 1,lm   
      DO j  = 1,jm
      DO i  = 1,im

       if(dusttx(i,j,l,nd).GT.0.) then
c number concentration
        dustnc(i,j,l,nd) = dusttx(i,j,l,nd)/pi*6/rop(nd)/
     .                     Dradi(nd)**3*exp(4.5*log(2.)**2)
       if(dustnc(i,j,l,nd).GT.0.) then
c median number diameter
         nmd = (dusttx(i,j,l,nd) /dustnc(i,j,l,nd)*6./pi/rop(nd))
     .         **0.33333* exp(1.5*log(2.)**2)
c pressure 
        phelp = pmid(l,i,j)*100.
          if (pmid(l,i,j).gt.100000.) then 
        phelp=99999.
          endif
c potential temperature, temperature
        te=pk(l,i,j)*t(i,j,l)
c pressure interpolation 
        np1=min(11,1+nint((10.-phelp/10000.)-0.499))  !pressure
        np1=max(1,np1)
        np2=min(11,np1+1)
c radii interpolation
        nh1=max(1,nint((nmd/Rrange)+0.499))      !median diameter 
        nh1=min(klo,nh1)              
        nh2=min(klo,nh1+1)
        px=((11-np1)*10000.- phelp)/10000.
        hx=((nh1*Rrange+md_look(1)) - nmd)/Rrange
        hp1=px*lookS(np2,nh1,il)+(1.-px)*lookS(np1,nh1,il)
        hp2=px*lookS(np2,nh2,il)+(1.-px)*lookS(np1,nh2,il)
        klook=hx*hp1+(1.-hx)*hp2

        if  (dustnc(i,j,l,nd).gt.1000..and.dustnc(i,j,l,nd).lt.(1.e30)) 
     *       then
        rxtnox(i,j,l,nd,il) = klook* dustnc(i,j,l,nd)
     .              / (287.054 * te / (pmid(l,i,j)*100.))
        else
        rxtnox(i,j,l,nd,il) = 0.    
        endif

        else
        rxtnox(i,j,l,nd,il) = 0.
        endif
        ENDIF
      ENDDO ! i
      ENDDO ! j
      ENDDO ! l
c      print*,' KRATE NR:,',il,'D TRA',nd,'SAHARA',rxtnox(36,28,1,nd,il)
c      print*,' STAUB', il,dustnc(36,28,1,nd)
      ENDDO ! nd
      
#ifdef TRACERS_DUST
      DO nd = 5,ndtr  !1,ndtr
        krate(:,:,:,1,il) = krate(:,:,:,1,il) + rxtnox(:,:,:,nd,il)
      ENDDO
        krate(:,:,:,2,il) = rxtnox(:,:,:,5,il)
#else
      DO nd = 1,ndtr  !1,ndtr
        krate(:,:,:,1,il) = krate(:,:,:,1,il) + rxtnox(:,:,:,nd,il)
      ENDDO
      DO nd = 1,5  !1,ndtr
        krate(:,:,:,2,il) = krate(:,:,:,2,il) + rxtnox(:,:,:,nd,il)
      ENDDO
#endif
        krate(:,:,:,3,il) = rxtnox(:,:,:,6,il)
        krate(:,:,:,4,il) = rxtnox(:,:,:,7,il)
        krate(:,:,:,5,il) = rxtnox(:,:,:,8,il)
      ENDDO ! il

c      print*,' KRATE NR:,', krate(36,28,1,:,:)

      return
      end subroutine 
#endif   ! TRACERS_SPECIAL_Shindell


#ifdef TRACERS_AEROSOLS_Koch
      SUBROUTINE SULFDUST
!
! Version 1.   (version 2 needs to be written... without integration over ndr)
!
!-----------------------------------------------------------------------
!   Computation of heterogeneous reaction rates on dust aerosol surfaces
!   Susanne E Bauer, 2003
!-----------------------------------------------------------------------

      USE MODEL_COM, only : im,jm,lm     ! dimensions (72, 46, 12)
     $                     ,jday         ! time in ITU
     $                     ,t            ! potential temperature (C)
     $                     ,q            ! saturatered pressure

      USE TRACER_COM                   !, only : rxts,trm
      USE CONSTANT,   only:  lhe       ! latent heat of evaporation at 0 C
      USE GEOM,       only:  bydxyp
      USE DYNAMICS,   only:  byam ,pmid,pk   ! midpoint pressure in hPa (mb)
      USE CONSTANT,   only:  pi, avog, gasc

      IMPLICIT NONE
!-----------------------------------------------------------------------
!       ... Dummy arguments
!-----------------------------------------------------------------------
      real                   :: erf
      integer, parameter     :: ndtr = 8  ! # dust bins
      real*8                 :: rxt(im,jm,lm,ndtr) !rate for single dust tracers
      real*8                   :: dusttx(im,jm,lm,ndtr) ! dust (kg/kg) 
      real                  :: dustnc(im,jm,lm,ndtr) ! dust number concentr.
!-----------------------------------------------------------------------
!       ... Look up variables
!----------------------------------------------------------------------- 
      integer, parameter :: klo = 1000
      integer ip,imd,np1,np2,nh1,nh2
      real klook,phelp,nmd
      real look_p, look_t,hx,px,hp1,hp2
!-----------------------------------------------------------------------
!       ... Local variables
!-----------------------------------------------------------------------
      integer :: i, j, k, nd, l ,ll
      integer, parameter :: ktoa = 300
! 1-SO2
c      real, parameter :: alph1  = 0.0001 !uptake coeff of Rossi EPFL (independent of humidity)
      real, parameter :: alph1  = 0.000001 !uptake coeff for SO2: RH < 60 %
      real, parameter :: alph2  = 0.0001    !uptake coeff for SO2: RH > 60 %
c      real, parameter :: alph1  = 0.1 !reactive uptake coefficient for HNO3
      real, parameter :: mQ1    = 64./1000.    ! kg/mol SO2
c      real, parameter :: mQ1    = 63./1000. ! kg/mol HNO3
      real, parameter :: xx    = 0.  !correction factor anisotropic movement
      real, parameter :: Bolz  = 1.3807e-23 !Boltzmann kg m2/s2 K molec.
      real, parameter :: Mgas  = 28.97 /1000. ! Molekular Gewicht Luft
      real, parameter :: Diaq  = 4.5e-10      ! m Molecul Diameter
C**** functions
      real*8 :: QSAT,RH,te,temp

      real :: Kn(rhet), Mdc(rhet), Kdj(2)
      real :: lamb(rhet), wrk(rhet),VSP(rhet)
      real :: lsig0,drada,dn,Roh!,te,temp

      logical, save             :: entereda = .false.
      real, save, dimension(ktoa) :: rada
      real, save                :: look(11,klo,2),Rrange,md_look(klo)

!-----------------------------------------------------------------
!     Dust variables: Radius of dust particles in the 8 size bins
!-----------------------------------------------------------------

      real, parameter  :: Dradi(8)=(/0.15e-6,0.25e-6,0.4e-6,0.8e-6,
     $                               1.5e-6,2.5e-6,4.e-6,8.e-6/)
      real, parameter  :: rop(8) = (/ 2500.,  2500.,  2500.,  2500.
     $                              , 2600.,  2600.,  2600.,  2600./)
!-----------------------------------------------------------------
!    1000 Intervals for Radius = 0.01ym ->10ym
!-----------------------------------------------------------------
!      Integration of radius:  0.01 ym to 10 ym
       rada(1) = 0.01E-6    ! smallest radius
       drada   = 0.1E-6     ! delta radius


      if (.not. entereda) then
      entereda = .true.
      PRINT*, 'CALCULATING LOOK UP TABLE FOR HETEROGENEOUS CHEMISTY'
      DO i   = 2, ktoa
      rada(i) = rada(i-1) + drada
      END DO

c      entereda = .true.

      lsig0 = LOG(2.)

!-----------------------------------------------------------------
!     SO2 + DUSTM =>    Dust Aerosol Reaction
!-----------------------------------------------------------------

      wrk(1)=  (mQ1 + Mgas) / mQ1 
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!                         LOOK UP TABLE
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       md_look(1) = 1.E-10       ! smallest md
       Rrange=5.E-8

       DO i   = 2, klo
       md_look(i) = md_look(i-1) +  Rrange 
c       write(*,*) 'MD_LOOK ', i, md_look(i)
       END DO

      DO ip  = 1, 11  !pressure from 1000 to 0 hPa 

      look_p=max(0.001,1.1-ip*0.1)         !atmosphere (minimum is 1 hPa)
      look_t=max(210.,288.*(look_p/1.)**((1.40-1)/1.40))
      look_p=look_p * 100000.              ! pressure in Pa

       Roh    = look_p / look_t / 287.
C Molecular diffusion coefficient for a trace gas in air [ m/s ]
       Mdc(1)  = 3. / (8.* Avog * Roh * (Diaq**2.))
       Mdc(1)  = Mdc(1) * SQRT( ((gasc*look_t*Mgas)/(2.*pi))*wrk(1))

C thermal velocity of a trace gas molecule [m/s2]
       VSP(1)  = SQRT((8. * Bolz * look_t)/(Pi * mQ1/ Avog))

C lamb  mean free pathway  [m]
       lamb(1)   = 3. * Mdc(1)/VSP(1)

C Loop over radius

      DO imd  = 1,klo
       look(ip,imd,:) = 0.

      DO k = 1,ktoa-1
C Knudsen Number 
       Kn(1)= lamb(1) / rada(k) ! Radius in [m]

C Mass Transfer Coefficient
c RH < 60 %
       Kdj(1) =(4. * pi * rada(k) * Mdc(1)) 
     .      / (1. + Kn(1) * (xx + 4 *(1.- alph1)/(3. *alph1)))
c RH > 60 %
       Kdj(2) =(4. * pi * rada(k) * Mdc(1)) 
     .      / (1. + Kn(1) * (xx + 4 *(1.- alph2)/(3. *alph2)))

C Number distribution
       dn=abs(erf(log( rada(k)/md_look(imd)) / lsig0
     .       /sqrt(2.))-erf(log(rada(k+1)/ 
     .       md_look(imd)) / lsig0 /sqrt(2.)))/2.

C Net removal rate [s-1]
       
       look(ip,imd,1)= look(ip,imd,1) + Kdj(1)  * dn    
       look(ip,imd,2)= look(ip,imd,2) + Kdj(2)  * dn    

      END DO              !radius loop
      END DO              !median diameter loop
      END DO              !pressure loop

      ENDIF

  
c--------------------------------------------------------------

c  Or use online dust 
#ifdef TRACERS_DUST      
c      CALL READDUST(dusttx)
c      print*, 'OFF-LINE DUST TRACER'
c      print*,'MAX OFF LINE:  ',maxval(dusttx(:,:,:,:))
      
c
      dusttx(:,:,:,:)= 0.
      
      DO l  = 1,lm   
      DO j  = 1,jm   
      dusttx(:,j,l,5)= trm(:,j,l,n_clay) * byam(l,:,j)* bydxyp (J)
      dusttx(:,j,l,6)= trm(:,j,l,n_silt1)* byam(l,:,j)* bydxyp (j)
      dusttx(:,j,l,7)= trm(:,j,l,n_silt2)* byam(l,:,j)* bydxyp (j)
      dusttx(:,j,l,8)= trm(:,j,l,n_silt3)* byam(l,:,j)* bydxyp (j)
      enddo
      enddo
#else
c Read dust fields from files
      CALL READDUST(dusttx)
      print*, 'OFF-LINE DUST TRACER'
#endif
c--------------------------------------------------------------
c--------------------------------------------------------------

c INTERPOLATION FROM LOOK UP TABLES
 
C Net removal rate for SO2 [s-1]  

#ifdef TRACERS_DUST
      DO nd = 5,ndtr    ! Loop over dust tracers
#else
      DO nd = 1,ndtr    ! Loop over dust tracers
#endif
      DO l  = 1,lm   
      DO j  = 1,jm
      DO i  = 1,im

c number concentration
        dustnc(i,j,l,nd) = dusttx(i,j,l,nd)/pi*6/rop(nd)/
     .                     Dradi(nd)**3*exp(4.5*log(2.)**2)
c median number diameter
        if(dustnc(i,j,l,nd).GT.0.) then
        nmd = (dusttx(i,j,l,nd) /dustnc(i,j,l,nd)*6./pi/rop(nd))
     .         **0.33333* exp(1.5*log(2.)**2)
c pressure 
        phelp = pmid(l,i,j)*100.
          if (pmid(l,i,j).gt.100000.) then 
        phelp=99999.
          endif
c potential temperature, temperature
        te=pk(l,i,j)*t(i,j,l)
c compute relative humidity
        RH=Q(i,j,l)/QSAT(te,lhe,pmid(l,i,j))    !temp in K, pres in mb
        IF(RH.LT.0.6) ll = 1
        IF(RH.GE.0.6) ll = 2
c pressure interpolation 
        np1=min(11,1+nint((10.-phelp/10000.)-0.499))  !pressure
        np1=max(1,np1)
        np2=min(11,np1+1)
c radii interpolation
        nh1=max(1,nint((nmd/Rrange)+0.499))      !median diameter 
        nh1=min(klo,nh1)              
        nh2=min(klo,nh1+1)
        px=((11-np1)*10000.- phelp)/10000.
        hx=((nh1*Rrange+md_look(1)) - nmd)/Rrange
        hp1=px*look(np2,nh1,ll)+(1.-px)*look(np1,nh1,ll)
        hp2=px*look(np2,nh2,ll)+(1.-px)*look(np1,nh2,ll)
        klook=hx*hp1+(1.-hx)*hp2

        if  (dustnc(i,j,l,nd).gt.1000..and.dustnc(i,j,l,nd).lt.(1.e30)) 
c        if  (dustnc(i,j,l,nd).gt.1000.) 
     *       then
        rxt(i,j,l,nd) = klook* dustnc(i,j,l,nd)
     .              / (287.054 * te / (pmid(l,i,j)*100.))
        else
        rxt(i,j,l,nd) = 0.    
        endif

        else
        rxt(i,j,l,nd) = 0.
        endif
      ENDDO ! i
      ENDDO ! j
      ENDDO ! l
      ENDDO ! nd

         rxts(:,:,:) = 0.

#ifdef TRACERS_DUST
      DO nd = 5,ndtr  !1,ndtr
        rxts(:,:,:) = rxts(:,:,:) + rxt(:,:,:,nd)
      ENDDO
        rxts1(:,:,:) = rxt(:,:,:,5)
#else
      DO nd = 1,ndtr  !1,ndtr
        rxts(:,:,:) = rxts(:,:,:) + rxt(:,:,:,nd)
      ENDDO
      DO nd = 1,5  !1,ndtr
        rxts1(:,:,:) = rxts1(:,:,:) + rxt(:,:,:,nd)
      ENDDO
#endif
        rxts2(:,:,:) = rxt(:,:,:,6)
        rxts3(:,:,:) = rxt(:,:,:,7)
        rxts4(:,:,:) = rxt(:,:,:,8)


      end subroutine sulfdust

      SUBROUTINE READDUST(dusttx)  
                                  
!@sum     READ IN MONTHLY MEAN DESERT DUST (CLAY,SILT) CONCENTRATION (KG/KG) 
!@auth    Original by Ina Tegen, Modified by Jan     Modified by Susana 
!@ver     05/15/03
!@calls   -
!@cont    List of program units contained within a file or module + other info.
!@fun     FUNCNAME Denotes a function
!@param   PARAMNAME Denotes a parameter (in the FORTRAN sense)
!@var     VARNAME Denotes a variable (in the FORTRAN sense)
!@dbparam VARNAME Denotes a database parameter (in the modelE sense)
!@nlparam VARNAME Denotes a NAMELIST parameter (in the modelE sense)
!@+       Continuation line for @sum/@calls/@con
C     ------------------------------------------------------------------ 

      USE MODEL_COM, only    : im,jm,lm,jday
      USE FILEMANAGER, only  : openunit, closeunit
      IMPLICIT NONE

      integer, parameter     :: ndtr = 8  ! # dust tracers
      integer nn,n,j,i,l,ma,mb
      REAL*8    wmb,wma,xmo        ! alles integer ???
      integer odust_trcl
      real*8, intent(inout)       :: dusttx(im,jm,lm,ndtr) ! dust (kg/kg) 

      real*4    ::  OFFDUST(im,jm,9,ndtr,12), DUST(im,jm,lm,ndtr,12)
                
      call openunit ("OFFDUST",odust_trcl,.TRUE.,.TRUE.)
      READ(odust_trcl) OFFDUST
      call closeunit(odust_trcl)

      do nn=1,12
         do n=1,ndtr
            do j=1,jm
               do i=1,im                                                 
                  do l=1,6                                               
                     DUST(i,j,l,n,nn)=OFFDUST(i,j,l,n,nn)                    
                  end do                                                 
      DUST(i,j,7,n,nn)=(30.*OFFDUST(i,j,6,n,nn)+45.
     .                 *OFFDUST(i,j,7,n,nn))/75.  
      DUST(i,j,8,n,nn)=OFFDUST(i,j,7,n,nn)  
      DUST(i,j,9,n,nn)=OFFDUST(i,j,8,n,nn)             
      DUST(i,j,10,n,nn)=(30.*OFFDUST(i,j,8,n,nn)+10.
     .                  *OFFDUST(i,j,9,n,nn))/40.  
      DUST(i,j,11,n,nn)=OFFDUST(i,j,9,n,nn)         
      DUST(i,j,12,n,nn)=OFFDUST(i,j,9,n,nn)                                
                  do l=13,lm                                               
                     DUST(i,j,l,n,nn)=OFFDUST(i,j,9,n,nn)                    
                  end do    
               end do                                                    
            end do                                                       
         end do                
      end do                                                        
C
C------------------------                                          
      ENTRY DAYDST                                                 
C------------------------
c     XMO=(JDAY+14.75)/30.5                                        
      XMO=(JDAY+15.25)/30.5                                        
      MA=XMO                                                       
      MB=MA+1                                                      
      WMB=XMO-MA                                                   
      WMA=1.0-WMB                                                  
      IF(MA.LT.1)  MA=12                                           
      IF(MB.GT.12) MB=1                                            
C                                                                  

      DUSTTX(:,:,:,:)=WMA*DUST(:,:,:,:,MA)+WMB*DUST(:,:,:,:,MB) 

      RETURN
      END
#endif   ! TRACERS_AEROSOLS_Koch


