#include "rundeck_opts.h"
      SUBROUTINE masterchem
!@sum masterchem main chemistry routine
!@auth Drew Shindell (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on masterchem0C9.2_M23p)
!@calls photoj,checktracer,Crates,Oxinit,HOxfam,NOxfam,chemstep
C
C PLEASE SEE THE WARNINGS IN THE STRATOSPHERIC OVERWRITE SECTION.
c
C**** GLOBAL parameters and variables:
c
      USE MODEL_COM, only: q,JDAY,IM,JM
      USE DYNAMICS, only: pedn
      USE RADNCB, only : ALB,COSZ1
      USE GEOM, only : BYDXYP, DXYP, LAT_DG
      USE FLUXES, only : tr3Dsource
      USE TRACER_COM, only: n_Ox,n_NOx,n_N2O5,n_HNO3,n_H2O2,n_CH3OOH,
     &                  n_HCHO,n_HO2NO2,n_CO,n_CH4,n_PAN,n_Isoprene,
     &                  n_AlkylNit,n_Alkenes,n_Paraffin
      USE CONSTANT, only: radian
      USE TRACER_DIAG_COM, only : jls_N2O5sulf,tajls
      USE TRCHEM_Shindell_COM
c
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
!@param by35 1/35 used for spherical geometry constant
      REAL*8, PARAMETER :: by35=1.d0/35.d0
!@var FASTJ_PFACT temp factor for vertical pressure-weighting
!@var FACT1 temp variable for start overwrite
!@var bydtsrc reciprocal of the timestep dtsrc
!@var imonth month index for Ox strat correction factor
!@var m dummy loop variable for Ox strat correction factor
!@var local logical for error checking 
      REAL*8 FASTJ_PFACT, FACT1, bydtsrc
      INTEGER imonth, m
      LOGICAL error
!@var I,J,L,N,igas,inss,LL,Lqq,JJ,J3 dummy loop variables
      INTEGER igas,LL,I,J,L,N,inss,Lqq,JJ,J3
c
C Note: mass2vol is no longer hardcoded here, since it does not have
C to be redefined each hour (It is now done in subroutine init_tracer
C in TRACERS_DRV.f.)
C
C++++ First, some INITIALIZATIONS :
      bydtsrc = 1./dtsrc
      BYFJM=1./float(JM)
C reset change due to chemistry to zero:
      change(:,:,:,:) = 0.
C
C set surface albedo variable used in fastj, based on ALB(..1)
C from radiation code:
      SALBFJ(:,:)= ALB(:,:,1)
C
c Set "chemical time step". Really this is a method of applying only
c a fraction of the chemistry change to the tracer mass for the first
c 30 hours.  That fraction is: dt2/dtscr.  E.g. in the first hour it
c is (dtsrc/24)/dtsrc = 1/24th of the chemistry change is applied.
c This is to work around initial instabilities.
c
      if(Itime-ItimeI.le.3)then
        dt2=dtsrc/24.            ! 150.
      elseif(Itime-ItimeI.gt.3.and.Itime-ItimeI.le.6)then
        dt2=dtsrc/12.            ! 300.
      elseif(Itime-ItimeI.gt.6.and.Itime-ItimeI.le.11)then
        dt2=dtsrc/6.             ! 600.
      elseif(Itime-ItimeI.gt.11.and.Itime-ItimeI.le.30)then
        dt2=dtsrc/2.4d0          ! 1500.
      elseif(Itime-ItimeI.gt.30)then
        dt2=dtsrc                ! 3600
      endif
c
c     Calculate new photolysis rates every n_phot main timesteps
      MODPHOT=MOD(Itime-ItimeI,n_phot)
C****
C**** CALCULATE TX, THE REAL TEMPERATURE
C**** (note this section is already done in DIAG.f)
      DO L=1,LM
        TX(1,1,L)=T(1,1,L)*PK(L,1,1)
        TX(1,JM,L)=T(1,JM,L)*PK(L,1,JM)
        DO I=2,IM
          TX(I,1,L)=TX(1,1,L)
          TX(I,JM,L)=TX(1,JM,L)
        END DO
        DO J=2,JM-1
        DO I=1,IM
          TX(I,J,L)=T(I,J,L)*PK(L,I,J)
        END DO
        END DO
      END DO
C
!$OMP  PARALLEL DO PRIVATE (FASTJ_PFACT, LL, 
!$OMP*    I,igas,inss, J, L,Lqq, N,error )

      DO J=1,JM                          ! >>>> MAIN J LOOP BEGINS <<<<
c     ILIMT=IM
c     if(J.eq.1.or.J.eq.JM)ILIMT=1
C
      DO I=1,IM                          ! >>>> MAIN I LOOP BEGINS <<<<
      DO L=1,LM
c      save presure and temperature in local arrays:
       pres(L)=PMID(L,I,J)
       ta(L)=TX(I,J,L)
c      calculate M and set fixed ratios for O2 & H2:
       y(nM,L)=pres(L)/(ta(L)*1.38d-19)
       y(nO2,L)=y(nM,L)*pfix_O2
       y(nH2,L)=y(nM,L)*pfix_H2
c
c      Tracers (converted from mass mixing ratio to number density)
       do igas=1,ntm
         y(igas,L)=trm(I,J,L,igas)*y(nM,L)*mass2vol(igas)*
     *   BYDXYP(J)*BYAM(L,I,J)
       enddo
c      If desired, fix the methane concentration used in chemistry
       if(fix_CH4_chemistry) THEN
         if(J.lt.JEQ)then ! SH
           y(n_CH4,L)=y(nM,L)*pfix_CH4_S
         else             ! NH
           y(n_CH4,L)=y(nM,L)*pfix_CH4_N
         endif
       end if
c
c Limit N2O5 number density:
       if(y(n_N2O5,L).lt.1.)y(n_N2O5,L)=1.0
c Set H2O, based on Q:
       y(nH2O,L)=Q(I,J,L)*MWabyMWw*y(nM,L)
c
       if(Itime.eq.ItimeI)then !initial startup
        y(nAldehyde,L)=y(nM,L)*pfix_Aldehyde
       else
        y(nAldehyde,L)=yAldehyde(I,J,L)
       endif
c      set [NO]=0 for first HOx calc, NO2 = NOx
c      set reactive species for use in Oxinit & nighttime NO2
       y(nNO2,L)     =y(n_NOx,L)*pNOx(I,J,L)
       y(nNO,L)      =y(n_NOx,L)*(1.-pNOx(I,J,L))
       y(nO3,L)      =pOx(I,J,L)*y(n_Ox,L)
       y(nCH3O2,L)   =yCH3O2(I,J,L)
       y(nC2O3,L)    =yC2O3(I,J,L)
       y(nXO2,L)     =yXO2(I,J,L)
       y(nXO2N,L)    =yXO2N(I,J,L)
       y(nRXPAR,L)   =yRXPAR(I,J,L)
       y(nROR,L)     =yROR(I,J,L)
      END DO ! L

C For solar zenith angle, we now use the arccosine of the COSZ1
C from the radiation code, which is the cosine of the solar zenith 
C angle averaged over the physics time step.
C If the solar zenith angle (sza) from the radiation code is > 90 deg,
C (and hence COSZ1 is set to 0), recalculate it with get_sza routine:
      IF(COSZ1(I,J).eq.0.) THEN
        call get_sza(I,J,sza)
      ELSE
        sza = acos(COSZ1(I,J))*byradian
      END IF
c
c     update radiation and temperatures, call PHOTOLYSIS every
c     [desired number] of hours:
      if(MODPHOT.eq.0)then             !  >>>> PHOTOLYSIS IF BEGIN <<<<
c      additional SUNLIGHT criterion (see also fam chem criterion):
       if((SALBFJ(I,J).ne.0.).AND.(sza.lt.szamax))then!>>SUNLIGHT<<<
c
c       define column temperatures to be sent to FASTJ:
        TFASTJ = ta
c
c       define pressures to be sent to FASTJ (centers):
        DO LL=2,2*LM,2
          PFASTJ(LL) = PMID(LL/2,I,J)
        END DO
c       define pressures to be sent to FASTJ (edges):
        PFASTJ(1)=PEDN(1,I,J)
        DO LL=3,(2*LM)+1,2
          PFASTJ(LL) = PEDN((LL+1)/2,I,J)
        END DO
        PFASTJ((2*LM)+2) = 0.1d0*PFASTJ((2*LM)+1)
#ifdef TRACERS_SPECIAL_Shindell
C       This is a fudge, so that we don't have to get mesosphere data:
        PFASTJ((2*LM)+2) = 0.00058d0 ! for 23 layer model...
#endif
c
c Interpolate O3 (in ppm) from bottom LS1-1 model sigma levels
c (ie those levels normally in troposphere) onto bottom 2*(LS1-1) FASTJ
c levels. Above these levels fastj uses climatological (Nagatani) O3
c read in by chem_init.
c
c       define O3 to be sent to FASTJ (centers):
        DO LL=2,2*(LS1-1),2
          O3_FASTJ(LL)=y(nO3,LL/2)
        ENDDO
c       define O3 to be sent to FASTJ (edges):
        O3_FASTJ(1)=y(nO3,1)*O3_1_fact ! see parameter declaration...
        DO LL=3,2*(LS1-1)-1,2
c         interpolation factor, based on pressure:
          FASTJ_PFACT=(PFASTJ(LL)-PFASTJ(LL-1))/
     &    (PFASTJ(LL+1)-PFASTJ(LL-1))
          O3_FASTJ(LL)=y(nO3,(LL-1)/2)+
     &    (y(nO3,((LL-1)/2)+1)-y(nO3,(LL-1)/2))*FASTJ_PFACT
c         lower limit on O3:
          IF(O3_FASTJ(LL).LT.0.) O3_FASTJ(LL)=-O3_FASTJ(LL)
        ENDDO
c       Convert to ppm (units used in fastj)
        DO LL=1,2*(LS1-1)
          O3_FASTJ(LL)=O3_FASTJ(LL)/(PFASTJ(LL)*2.55d10)
        ENDDO
C
C CALL THE PHOTOLYSIS SCHEME:
C
       call photoj(I,J)
c
C And fill in the photolysis coefficients: ZJ --> ss:

        DO L=1,JPNL
          do inss=1,JPPJ
           ss(inss,L,I,J)=zj(L,inss)*
     &     by35*SQRT(1.224d3*(cos(ABS(LAT_DG(J,1))*radian))**2.+1.d0)
          enddo
        END DO

       endif                               ! >>>> SUNLIGHT IF END <<<<
      endif                             !  >>>> PHOTOLYSIS IF END <<<<
c
c     calculate the chemical reaction rates:
      call Crates (I,J)      
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
cc    Main chemistry calculations    cc
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
c     Family partitioning:
c
      if((SALBFJ(I,J).ne.0.).AND.(sza.lt.szamax))then ! S U N L I G H T
       call Oxinit(LS1-1,I,J)
       call HOxfam(LS1-1,I,J)
       call NOxfam(LS1-1,I,J)
c
C Some Chemistry Diagnostics:
      if(prnchg.and.J.eq.jprn.and.I.eq.iprn)then
       l=lprn
       write(6,*) ' '
       write(6,*) 'Family ratios at I,J,L: ',i,j,l
       write(6,*) 'OH/HO2 = ',y(nOH,l)/y(nHO2,l)
       write(6,*) 'O/O3 = ',y(nO,l)/y(nO3,l)
       write(6,*) 'O1D/O3 = ',y(nO1D,l)/y(nO3,l),
     *  '  J(O1D) = ',ss(2,l,I,J)
c     *  ,' J(O3) = ',ss(3,l,I,J)
       write(6,*) 'NO/NO2 = ',y(nNO,l)/y(nNO2,l),
     *  '   J(NO2) = ',ss(1,l,I,J)
c       write(6,*) 'pHOx, pNOx, pOx = ',pHOx(I,J,l),
c     * pNOx(I,J,l),pOx(I,J,l)
       write(6,*) 'conc OH = ',y(nOH,l)
       write(6,*)'sun, SALBFJ,sza,I,J,Itime= ',SALBFJ(I,J),sza,I,J,Itime
       do inss=1,JPPJ
        write(6,195) ' J',inss,ay(ks(inss)),' = ',
     *   (ss(inss,Lqq,I,J),Lqq=1,11)
       enddo
       write(6,196) ' RCloud',(RCLOUDFJ(Lqq,I,J),Lqq=1,11)
       write(6,196) ' Ozone ',(y(nO3,Lqq),Lqq=1,11)
       write(6,*) ' '
      endif
 195  format(a2,i2,1x,a8,a3,11(1x,e9.2))
 196  format(a7,9x,11(1x,e9.2))
c
cc    Non-family chemistry:
c
      call chemstep(I,J,change)
c
      else                           !             D A R K N E S S
c
C*****************************************************************
c      >> N2O5 sink on sulfate aerosol :: <<
C      Raw (ppt) sulfate data converted to surface area (cm2/cm3)
C      assuming a monodispersed aerosol with radius 0.078 microns
C      (Dentener and Crutzen, 1993). Also assuming density = 1.1 g/cm3,
C      Mr(sulfate) = 98.0g; 1ppt=2.55E7 molecules/cm3
C
C       Conversion factor used (ppt->surface area)
C       =2.55E7*98.*3./(6.022E23*1.1*1.E-5)
C
C      To evaluate 'k' for N2O5 + H2O -> 2HNO3,
C       assume k = GAMMA.SA.v / 4 (Dentener and Crutzen, 1993)
C
C      GAMMA = 0.1, SA = given, v = molecular velocity (cm/s)
C      where v  = SQRT (8.Kb.T / PI*M); Kb = 1.38062E-23;
C      T = Temp (K); M = mass N2O5 (Kg)
C*****************************************************************
c
       do L=1,LS1-1 ! (troposphere)
         pfactor=dxyp(J)*AM(L,I,J)/y(nM,L)
         bypfactor=1.D0/pfactor
         RVELN2O5=SQRT(TX(I,J,L)*RKBYPIM)*100.
C        Calculate sulfate sink, and cap it at 90% of N2O5:
         wprod_sulf=DT2*sulfate(I,J,L)*y(n_N2O5,L)*RGAMMASULF*RVELN2O5
         if(wprod_sulf.gt.0.9*y(n_N2O5,L))wprod_sulf=0.9d0*y(n_N2O5,L)
         prod_sulf=wprod_sulf*pfactor
         TAJLS(J,L,jls_N2O5sulf)=TAJLS(J,L,jls_N2O5sulf)
     &     -(prod_sulf*bymass2vol(n_N2O5))
C
C*****************************************************************
c        g signifies gas phase
c        while prod_sulf and wprod_sulf are sulfate rxn
c        wprods are in molecules/cm3/s
c        prods/mass2vol are in mass units to add to tracers
c
c        NO3 amounts are a function of reaction 7 (NO2 + O3 -> NO3),
c        24, 25 (leave out 28, 0.9*32, 46, 52 outside NOx family)
c        NO2, similarly leave out 29, 45, and 46.
c        Keep NOx unchanged as this is only intrafamily
C*****************************************************************
C
c        define O3 from Ox:
         y(nO3,L)=y(n_Ox,L)*pOx(I,J,L)
c
c        calculate change in NO3:
         dNO3=rr(7,L)*y(nNO2,L)*y(n_Ox,L)-(rr(24,L)*y(nNO2,L)+
     &        2*rr(25,L)*yNO3(I,J,L))*yNO3(I,J,L) -(rr(36,L)
     &        *y(n_Alkenes,L)+rr(32,L)*y(n_Isoprene,L))*yNO3(I,J,L)
         dNO3=dNO3-(rr(28,L)*y(n_HCHO,L)+rr(52,L)*y(nNO2,L))
     &        *yNO3(I,J,L)+rr(46,L)*y(n_N2O5,L)
         dNO3=dNO3*dt2!portions of dNO3 reflect fixes by dts 12/19/01
c
c        limit the change in NO3:
         if(-dNO3.gt.0.66*yNO3(I,J,L))dNO3=-0.66d0*yNO3(I,J,L)
c
c        apply the NO3 change limit the value to positive and 1/2 NOx:
         yNO3(I,J,L)=yNO3(I,J,L)+dNO3
         if(yNO3(I,J,L).lt.0.)yNO3(I,J,L)=0.
         if(yNO3(I,J,L).gt.y(n_NOx,L)*0.5)yNO3(I,J,L)=y(n_NOx,L)*0.5
c
c        calculate and limit NO2
         y(nNO2,L)=y(n_NOx,L)-yNO3(I,J,L)
         pNOx(I,J,L)=y(nNO2,L)/(y(n_NOx,L)+1.E-10)
         if(pNOx(I,J,L).gt.1)pNOx(I,J,L)=1.
         if(pNOx(I,J,L).lt.0.5)pNOx(I,J,L)=0.5
C
C        LOWER LIMIT ON N2O5:
         if(y(n_N2O5,L).le.1.)y(n_N2O5,L)=1.
C
C Calculate and limit gaseous changes to HNO3, HCHO, N2O5, Aldehyde,
C Alkenes, Isoprene, and AlkylNit:
C
         gwprodHNO3=(y(n_HCHO,L)*rr(28,L)+2.5d-15*
     &              yAldehyde(I,J,L))*yNO3(I,J,L)*dt2
         if(gwprodHNO3.gt.0.5*y(n_NOx,L))gwprodHNO3=0.5*y(n_NOx,L)
         if(gwprodHNO3.gt.y(n_HCHO,L))gwprodHNO3=y(n_HCHO,L)
         gprodHNO3=gwprodHNO3*pfactor
C
         gwprodN2O5=(yNO3(I,J,L)*y(nNO2,L)*rr(52,L)-y(n_N2O5,L)
     &              *rr(46,L))*dt2
         if(gwprodN2O5.gt.0.25*y(n_NOx,L))gwprodN2O5=0.25*y(n_NOx,L)
         if(-gwprodN2O5.gt.0.5*y(n_N2O5,L))
     *        gwprodN2O5=-0.49d0*y(n_N2O5,L)
C
         changeAldehyde=(rr(36,L)*y(n_Alkenes,L)+rr(32,L)*
     &                  y(n_Isoprene,L)*0.12-2.5E-15*yAldehyde(I,J,L))*
     &                  yNO3(I,J,L)*dt2
         if(-changeAldehyde.gt.0.75*yAldehyde(I,J,L))changeAldehyde=
     &   0.75*yAldehyde(I,J,L)
C
         changeAlkenes=(rr(32,L)*y(n_Isoprene,L)*0.45d0-rr(36,L)*
     &                y(n_Alkenes,L))*yNO3(I,J,L)*dt2+(rr(31,L)*
     &                y(n_Isoprene,L)-rr(35,L)*y(n_Alkenes,L))*
     &                y(nO3,L)*dt2
         if(-changeAlkenes.gt.0.75*y(n_Alkenes,L))changeAlkenes=
     &   0.75*y(n_Alkenes,L)
C
         changeIsoprene=-(rr(32,L)*yNO3(I,J,L)
     &                  +rr(31,L)*y(nO3,L))*y(n_Isoprene,L)*dt2
         if(-changeIsoprene.gt.0.75*y(n_Isoprene,L))changeIsoprene=
     &   0.75*y(n_Isoprene,L)
C
         changeHCHO=(rr(36,L)*y(n_Alkenes,L)+
     &              rr(32,L)*y(n_Isoprene,L)*0.03)*yNO3(I,J,L)*dt2
     &              -gwprodHNO3+(rr(31,L)*y(n_Isoprene,L)*0.9d0
     &              +rr(35,L)*y(n_Alkenes,L))*y(nO3,L)*0.64*dt2
C
         changeAlkylNit=rr(32,L)*y(n_Isoprene,L)*
     &   yNO3(I,J,L)*dt2*0.9d0
         if(-changeAlkylNit.gt.0.75*y(n_AlkylNit,L))changeAlkylNit=
     &   0.75*y(n_AlkylNit,L)
C
c        Convert some changes to molecules/cm3/s:
         changeHNO3=gwprodHNO3+2*wprod_sulf  !always positive
         changeNOx=-gwprodHNO3-2*gwprodN2O5-(0.9d0*rr(32,L)*
     &    y(n_Isoprene,L)+2.5E-15*yAldehyde(I,J,L))*yNO3(I,J,L)*dt2
         if(-changeNOx.gt.y(n_NOx,L))changeNOx=-.95d0*y(n_NOx,L)!dts9/01
         changeN2O5=gwprodN2O5-wprod_sulf
c
c Ensure nitrogen conservation (presumably dNOx<0, others >0):
c
         rlossN=0.
         rprodN=0.
         if(changeNOx.lt.0)then
          rlossN=rlossN+changeNOx
         else
          rprodN=rprodN+changeNOx
         endif
         if(changeHNO3.lt.0)then
          rlossN=rlossN+changeHNO3
         else
          rprodN=rprodN+changeHNO3
         endif
         if(changeN2O5.lt.0)then
          rlossN=rlossN+2*changeN2O5
         else
          rprodN=rprodN+2*changeN2O5
         endif
         if(changeAlkylNit.lt.0)then
          rlossN=rlossN+changeAlkylNit
         else
          rprodN=rprodN+changeAlkylNit
         endif
         if(rprodN.gt.rlossN)then
          ratioN=-rlossN/rprodN
          if(changeNOx.gt.0.)     changeNOx     =changeNOx     *ratioN
          if(changeHNO3.gt.0.)    changeHNO3    =changeHNO3    *ratioN
          if(changeN2O5.gt.0.)    changeN2O5    =changeN2O5    *ratioN
          if(changeAlkylNit.gt.0.)changeAlkylNit=changeAlkylNit*ratioN
         else
          ratioN=rprodN/(-rlossN)
          if(changeNOx.lt.0.)     changeNOx     =changeNOx     *ratioN
          if(changeHNO3.lt.0.)    changeHNO3    =changeHNO3    *ratioN
          if(changeN2O5.lt.0.)    changeN2O5    =changeN2O5    *ratioN
          if(changeAlkylNit.lt.0.)changeAlkylNit=changeAlkylNit*ratioN
         endif
C
C Apply Alkenes, AlkyNit, and Aldehyde changes here:
C
         y(n_Alkenes,L)=y(n_Alkenes,L)+changeAlkenes
         y(n_AlkylNit,L)=y(n_AlkylNit,L)+changeAlkylNit
         yAldehyde(I,J,L)=yAldehyde(I,J,L)+changeAldehyde
C
C Note: there is a lower limit of 1 placed on the resulting tracer mass
C from the following changes. This is to prevent negative tracer mass:
C
C -- HCHO --
c        Gas phase NO3 + HCHO -> HNO3 + CO yield of HCHO & CO
         change(I,J,L,n_HCHO)=changeHCHO*pfactor*bymass2vol(n_HCHO)
         if(-change(I,J,L,n_HCHO).gt.trm(I,J,L,n_HCHO))then
           change(I,J,L,n_HCHO)=-.95d0*trm(I,J,L,n_HCHO)
           changeHCHO=change(I,J,L,n_HCHO)*mass2vol(n_HCHO)*bypfactor
         endif
         IF((trm(i,j,l,n_HCHO)+change(i,j,l,n_HCHO)).lt.1.) THEN
           change(i,j,l,n_HCHO) = 1. - trm(i,j,l,n_HCHO)
           changeHCHO=change(I,J,L,n_HCHO)*mass2vol(n_HCHO)*bypfactor
         ENDIF
         wprodHCHO=changeHCHO
C -- CO --
         change(I,J,L,n_CO)=gprodHNO3*bymass2vol(n_CO)
         IF((trm(i,j,l,n_CO)+change(i,j,l,n_CO)).lt.1.)
     &   change(i,j,l,n_CO) = 1. - trm(i,j,l,n_CO)
         wprodCO=gwprodHNO3   ! <<< note
C -- HNO3 --  (HNO3 from gas and het phase rxns )
         change(I,J,L,n_HNO3)=changeHNO3*pfactor*bymass2vol(n_HNO3)
         IF((trm(i,j,l,n_HNO3)+change(i,j,l,n_HNO3)).lt.1.) THEN
           change(i,j,l,n_HNO3) = 1. - trm(i,j,l,n_HNO3)
           changeHNO3=change(I,J,L,n_HNO3)*mass2vol(n_HNO3)*bypfactor
         END IF
C -- N2O5 --  (N2O5 from gas and het phase rxns)
         change(I,J,L,n_N2O5)=changeN2O5*pfactor*bymass2vol(n_N2O5)
         IF((trm(i,j,l,n_N2O5)+change(i,j,l,n_N2O5)).lt.1.) THEN
           change(i,j,l,n_N2O5) = 1. - trm(i,j,l,n_N2O5)
           changeN2O5=change(I,J,L,n_N2O5)*mass2vol(n_N2O5)*bypfactor
         END IF
c -- NOx --   (NOx from gas phase rxns)
         change(I,J,L,n_NOx)=changeNOx*pfactor*bymass2vol(n_NOx)
         IF((trm(i,j,l,n_NOx)+change(i,j,l,n_NOx)).lt.1.) THEN
           change(i,j,l,n_NOx) = 1. - trm(i,j,l,n_NOx)
           changeNOx=change(I,J,L,n_NOx)*mass2vol(n_NOx)*bypfactor
         END IF
C -- Alkenes --  (Alkenes from gas phase rxns)
         change(I,J,L,n_Alkenes)=
     &   changeAlkenes*pfactor*bymass2vol(n_Alkenes)
         IF((trm(i,j,l,n_Alkenes)+change(i,j,l,n_Alkenes)).lt.1.)THEN
           change(i,j,l,n_Alkenes) = 1. - trm(i,j,l,n_Alkenes)
           changeAlkenes=change(I,J,L,n_Alkenes)*mass2vol(n_Alkenes)
     &     *bypfactor
         END IF
c -- Isoprene -- (Isoprene from gas phase rxns)
         change(I,J,L,n_Isoprene)=
     &   changeIsoprene*pfactor*bymass2vol(n_Isoprene)
         IF((trm(i,j,l,n_Isoprene)+change(i,j,l,n_Isoprene)).lt.1.)THEN
           change(i,j,l,n_Isoprene) = 1. - trm(i,j,l,n_Isoprene)
           changeIsoprene=change(I,J,L,n_Isoprene)*mass2vol(n_Isoprene)
     &     *bypfactor
         END IF
c -- AlkylNit -- (AlkylNit from gas phase rxns)
         change(I,J,L,n_AlkylNit)=
     &   changeAlkylNit*pfactor*bymass2vol(n_AlkylNit)
         IF((trm(i,j,l,n_AlkylNit)+change(i,j,l,n_AlkylNit)).lt.1.)THEN
           change(i,j,l,n_AlkylNit) = 1. - trm(i,j,l,n_AlkylNit)
           changeAlkylNit=change(I,J,L,n_AlkylNit)*mass2vol(n_AlkylNit)
     &     *bypfactor
         END IF
C
c Some More Chemistry Diagnostics:
         if(prnchg.and.J.eq.jprn.and.I.eq.iprn.and.L.eq.lprn)then
          write(6,*) 'dark, SALBFJ,sza,I,J,L,Itime= ',SALBFJ(I,J),sza,I
     &,J,L,Itime
          write(6,198) ay(n_NOx),': ',
     *         changeNOx,' molecules produced; ',
     *      100.*(changeNOx)/y(n_NOx,L),' percent of'
     *     ,y(n_NOx,L),'(',1.E9*y(n_NOx,L)/y(nM,L),' ppbv)'
          write(6,198) ay(n_HNO3),': ',
     *         changeHNO3,' molecules produced; ',
     *     100.*(changeHNO3)/y(n_HNO3,L),' percent of'
     *     ,y(n_HNO3,L),'(',1.E9*y(n_HNO3,L)/y(nM,L),' ppbv)'
          write(6,198) ay(n_N2O5),': ',
     *         changeN2O5,' net molec produced; ',
     *      100.*(changeN2O5)/y(n_N2O5,L),' percent of'
     *     ,y(n_N2O5,L),'(',1.E9*y(n_N2O5,L)/y(nM,L),' ppbv)'
          write(6,198) ay(n_N2O5),': ',
     *            gwprodN2O5,' molec prod fm gas;  ',
     *      100.*(gwprodN2O5)/y(n_N2O5,L),' percent of'
     *     ,y(n_N2O5,L),'(',1.E9*y(n_N2O5,L)/y(nM,L),' ppbv)'
          write(6,198) ay(n_HCHO),': ',
     *         wprodHCHO,' molecules produced; ',
     *      100.*(wprodHCHO)/y(n_HCHO,L),' percent of'
     *     ,y(n_HCHO,L),'(',1.E9*y(n_HCHO,L)/y(nM,L),' ppbv)'
          write(6,198) 'Aldehyde',': ',
     *         changeAldehyde,' molecules produced; ',
     *      100.*(changeAldehyde)/yAldehyde(I,J,L),' percent of'
     *     ,yAldehyde(I,J,L),'(',1.E9*yAldehyde(I,J,L)/y(nM,L),' ppbv)'
          write(6,198) 'Alkenes ',': ',
     *         changeAlkenes,' molecules produced; ',
     *      100.*(changeAlkenes)/y(n_Alkenes,L),' percent of'
     *     ,y(n_Alkenes,L),'(',1.E9*y(n_Alkenes,L)/y(nM,L),' ppbv)'
          write(6,198) 'Isoprene',': ',
     *         changeIsoprene,' molecules produced; ',
     *      100.*(changeIsoprene)/y(n_Isoprene,L),' percent of'
     *     ,y(n_Isoprene,L),'(',1.E9*y(n_Isoprene,L)/y(nM,L),' ppbv)'
          write(6,198) 'AlkylNit',': ',
     *         changeAlkylNit,' molecules produced; ',
     *      100.*(changeAlkylNit)/y(n_AlkylNit,L),' percent of'
     *     ,y(n_AlkylNit,L),'(',1.E9*y(n_AlkylNit,L)/y(nM,L),' ppbv)'
          write(6,199) 'NO2, NO3  = ',y(nNO2,L),yNO3(I,J,L)
         endif
 198  format(1x,a8,a2,e13.3,a21,f10.0,a11,2x,e13.3,3x,a1,f12.5,a6)
 199  format(1x,a20,2(2x,e13.3))
C
C Make sure nighttime chemistry changes are not too big:
C
         error=.false.
         if(changeNOx.lt.-1.E15.OR.changeNOx.gt.1.E15) then
          write(*,*) 'Big chg@ Itime,I,J,L,NOx ',Itime,I,J,L,changeNOx
          error=.true.
         end if
         if(changeHNO3.lt.-1.E15.OR.changeHNO3.gt.1.E15) then
          write(*,*) 'Big chg@ Itime,I,J,L,HNO3',Itime,I,J,L,changeHNO3
          error=.true.
         end if
         if(changeN2O5.lt.-1.E15.OR.changeN2O5.gt.1.E15) then
          write(*,*) 'Big chg@ Itime,I,J,L,N2O5',Itime,I,J,L,changeN2O5
          error=.true.
         end if
         if(wprodHCHO.lt.-1.E15.OR.wprodHCHO.gt.1.E15) then
          write(*,*)'Big chg@ Itime,I,J,L,HCHO',Itime,I,J,L,wprodHCHO
          error=.true.
         endif
         if(error)call stop_model('nighttime chem: big changes',255)
C
       enddo  ! troposphere loop
c
      endif                         ! >>> END sunlight/darkness IF <<<
C

      DO L=1,LS1-1  ! loop over troposphere again

C >> Lower limit on HO2NO2 of 1.0 <<
        if(trm(i,j,l,n_HO2NO2)+change(i,j,l,n_HO2NO2).lt.1.)
     &  change(i,j,l,n_HO2NO2) = 1. - trm(i,j,l,n_HO2NO2)
C Save chemistry changes for updating tracers in apply_tracer_3Dsource.
        DO N=1,NTM
          tr3Dsource(i,j,l,1,n) = change(i,j,l,n) * bydtsrc
        END DO
c Reset O3DLJI values for radiation (gcm):
        O3DLJI(L,J,I)=
     &  (trm(I,J,L,n_Ox)+change(I,J,L,n_Ox))*BYDXYP(J)*BYO3MULT

      END DO       ! end current troposphere loop
C
c
      END DO ! >>>> MAIN I LOOP ENDS <<<<
c
      END DO ! >>>> MAIN J LOOP ENDS <<<<
!$OMP END PARALLEL DO  

C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
cc    Stratospheric Overwrite of tracers                cc
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C W A R N I N G :Unfortunately, there is still a hardcoded dependence
C                on the verticle resolution in how we deal with the
C                CH4 overwriting... It is set up for the 23 layer model.
C                Not worth fixing, since we'll soon have stratospheric
C                chemistry instead.
C W A R N I N G :If there is ever stratospheric chemistry (i.e. the
C                'change' variable for L>LS1-1 is non-zero at this point
C                in the code), then the stratospheric changes below
C                should be altered.  Currently, they are functions of
C                tracer mass UNCHANGED by chemistry !!!!!!!
C
C     Make sure that change is zero:
      change = 0.
C
c Calculate non-weighted average of layer 6 and 7 CH4, for M23 models,
c in mixing ratio units, for every 3rd latitude:
c
      jj=0
      do j=(1+3),(JM-3),3
       jj=jj+1
       FACTj=(0.5*mass2vol(n_CH4)*1.d6*bydxyp(j))
       do i=1,im
         avg67(i,jj)=FACTj*(trm(i,j,6,n_CH4)*byam(6,i,j) +
     &   trm(i,j,7,n_CH4)*byam(7,i,j))
       end do
      end do
C     0.55866= 1/1.79 is 1/(obs. tropsph. CH4):
      r179m2v=0.55866d0*bymass2vol(n_CH4)

C Use Ox correction factor for the proper month:
      imonth= 1
      DO m=2,12
       IF((JDAY.LE.MDOFM(m)).AND.(JDAY.GT.MDOFM(m-1))) THEN
        imonth=m
        GOTO 217
       END IF
      END DO
 217  CONTINUE
C
      do L=LS1,LM                  ! >> BEGIN LOOP OVER STRATOSPHERE <<
       do j=1,jm
        j3=NINT(float(j)*19.*BYFJM)! index for CO
        IF(J3.eq.0) J3=1           ! index for CO
        jj=NINT((FLOAT(j)-1.)/3.)  ! index for CH4
        if(j.le.2) jj=1            ! index for CH4
        if(j.ge.45)jj=14           ! index for CH4
        do i=1,im
          if(L.le.15) then
            change(I,J,L,n_Ox)=
     &      OxIC(I,J,L)*corrOx(J,L-11,imonth) - trm(I,J,L,n_Ox)
          else
            change(I,J,L,n_Ox)= OxIC(I,J,L)   - trm(I,J,L,n_Ox)
          end if
C         note: layer 14 in M23 is between layers 8 and 9 in M9 model:
          FACT1=2.0d-9*DXYP(J)*am(L,I,J)*byam(14,I,J)
          change(I,J,L,n_NOx)=trm(I,J,L,n_Ox)*2.3d-4 - trm(I,J,L,n_NOx)
C      dts 12/19/01:NOx strat-trop flux too big, alter lower strat NOx:
          if((L.eq.LS1).or.(L.eq.LS1+1)) change(I,J,L,n_NOx)=
     &         0.9d0*trm(I,J,L,n_Ox)*2.3d-4-trm(I,J,L,n_NOx)
          change(I,J,L,n_N2O5)=  FACT1            - trm(I,J,L,n_N2O5)
          change(I,J,L,n_HNO3)=trm(I,J,L,n_Ox)*4.2d-3-trm(I,J,L,n_HNO3)
          change(I,J,L,n_H2O2)=  FACT1            - trm(I,J,L,n_H2O2)
          change(I,J,L,n_CH3OOH)=FACT1            - trm(I,J,L,n_CH3OOH)
          change(I,J,L,n_HCHO)=  FACT1            - trm(I,J,L,n_HCHO)
          change(I,J,L,n_HO2NO2)=FACT1*70.        - trm(I,J,L,n_HO2NO2)
C         note above: 70. = 1.4E-7/2.0E-9
          change(I,J,L,n_CO)=(COlat(J3)*COalt(L)*1.D-9*bymass2vol(n_CO)
     &    *AM(L,I,J)*DXYP(J))                   - trm(I,J,L,n_CO)
          change(I,J,L,n_PAN)     = FACT1*1.d-4 - trm(I,J,L,n_PAN)
          change(I,J,L,n_Isoprene)= FACT1*1.d-4 - trm(I,J,L,n_Isoprene)
          change(I,J,L,n_AlkylNit)= FACT1*1.d-4 - trm(I,J,L,n_AlkylNit)
          change(I,J,L,n_Alkenes) = FACT1*1.d-4 - trm(I,J,L,n_Alkenes)
          change(I,J,L,n_Paraffin)= FACT1*1.d-4 - trm(I,J,L,n_Paraffin)
c
c         Overwrite stratospheric ch4 based on HALOE obs for tropics
c         and extratropics and scale by the ratio of nearby lvls 6 & 7
c         mixing ratios to 1.79. The 12 "stratospheric" levels are
c         grouped into 4 categories: L={12,13,14} based on 100mb obs.
c         L={15,16,17} on 32mb obs, L={18,19,20} on 3.2mb obs and
c         L={21,22,23} on 0.32mb obs (average of mar,jun,sep,dec).
c
          CH4FACT=avg67(i,jj)*r179m2v
          IF((J.LE.16).OR.(J.GT.30)) THEN                ! extratropics
            IF((L.GE.LS1).AND.(L.LE.14)) CH4FACT=CH4FACT* 1.7d0!was1.44
            IF((L.GE.15).AND.(L.LE.17))  CH4FACT=CH4FACT*   1.130d0
            IF((L.GE.18).AND.(L.LE.20))  CH4FACT=CH4FACT*   0.473d0
            IF((L.GE.21).AND.(L.LE.LM))  CH4FACT=CH4FACT*   0.202d0
          ELSE IF((J.GT.16).AND.(J.LE.30)) THEN           ! tropics
            IF((L.GE.LS1).AND.(L.LE.14)) CH4FACT=CH4FACT*   1.620d0
            IF((L.GE.15).AND.(L.LE.17))  CH4FACT=CH4FACT*   1.460d0
            IF((L.GE.18).AND.(L.LE.20))  CH4FACT=CH4FACT*   0.812d0
            IF((L.GE.21).AND.(L.LE.LM))  CH4FACT=CH4FACT*   0.230d0
          END IF
          change(I,J,L,n_CH4)=
     &    (AM(L,I,J)*DXYP(J)*CH4FACT*1.d-6) - trm(I,J,L,n_CH4)
C
C Save stratosph change for updating tracer, apply_tracer_3Dsource:
          DO N=1,NTM
            tr3Dsource(i,j,l,2,n) = change(i,j,l,n) * bydtsrc
          END DO
C
        end do !i
       end do  !j
      end do                        ! >> END LOOP OVER STRATOSPHERE <<
c
      RETURN
      END SUBROUTINE masterchem
c
      SUBROUTINE Crates(I,J)
!@sum Crates calculate chemical reaction rates for each altitude,
!@+   using JPL 00.  Includes special calculations for pressure
!@+   dependent reactions. Specifically:
!@+   #13 CO+OH->HO2+CO2, #15 HO2+HO2->H2O2+O2, #16 OH+HNO3->H2O+NO3,
!@+   and reactions #29, and #42.
!@auth Drew Shindell (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on ds3ch4_master_0C4.2_M23)
c
C**** GLOBAL parameters and variables:
C
      USE MODEL_COM, only : ls1
      USE TRCHEM_Shindell_COM, only: nr2,nr3,nmm,nhet,ta,ea,rr,pe,
     &                          r1,sb,nst,y,nM,nH2O,ro,sn
c
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
C
!@var I,J passed horizontal position indicies
!@var dd,pp,fw,rkp,rk2,rk3M,nb,rrrr dummy "working" variables
!@var L,jj dummy loop variables
!@var byta reciprocal of the local temperature
      REAL*8 byta, dd, pp, fw, rkp, rk2, rk3M, rrrr
      INTEGER L,jj,nb
      INTEGER, INTENT(IN) :: I,J
C
      do L=1,LS1-1            !  >>> BEGIN TROPOSPHERE LOOP <<<
        byta=1./ta(L)
        do jj=1,nr2             ! bimolecular rates start
          IF(ea(jj).ne.0.) THEN
            rr(jj,L)=pe(jj)*exp(-ea(jj)*byta)
          ELSE
            rr(jj,L)=pe(jj)
          END IF
c         for #13, k=pe*(1+0.6*(Patm/1013)) Patm=[M]*(T*1.38E-19)
          if(jj.eq.13) rr(jj,L) =
     &    pe(jj)*(1.+0.6d0*((y(nM,L)*ta(L)*1.38d-19)/1013.))
c         for reaction #15, k=(kc+kp)fw, kc=rr
          if(jj.eq.15)then
            rkp=1.7d-33*y(nM,L)*exp(1000.*byta)
            fw=(1.+1.4d-21*y(nH2O,L)*exp(2200.*byta))
            rr(jj,L)=(rr(jj,L)+rkp)*fw
          endif
c         for #16, k=[pe*exp(-e(jj)/ta(l))]+k3[M]/(1+k3[M]/k2)
          if(jj.eq.16)then
            rk3M=y(nM,l)*1.90d-33*exp(725.*byta)
            rk2=4.10d-16*exp(1440.*byta)
            rr(jj,L)=rr(jj,L)+rk3M/(1.+(rk3M/rk2))
          endif
          if(jj.eq.29)rr(jj,L)=rr(jj,L)/y(nM,L)!PAN+M really PAN
          if(jj.eq.42)rr(jj,L)=rr(jj,L)/y(nM,L)!ROR+M really ROR
        end do                ! bimolecular rates end
c
        do jj=1,nr3           ! trimolecular rates start
         rr(nr2+jj,L)=y(nM,L)*ro(jj)*(300.*byta)**sn(jj)
c        if(r1(jj).ge.1.E-30)then
         if(sb(jj).ge.0.01)then !dts 3/29/02 alteration of above line
           dd=rr(nr2+jj,L)/(r1(jj)*(300.*byta)**sb(jj))
           pp=0.6d0**(1./(1.+(log10(dd))**2.))
           rr(nr2+jj,L)=(rr(nr2+jj,L)/(1.+dd))*pp
         endif
        end do                ! trimolecular rates end
c
        nb=nr2-nmm
        if(nmm.ge.1) then
          do jj=1,nmm         ! monomolecular rates start
           rrrr=exp(0.5*ea(jj+nb)*byta) !0.5 for precision,cor next lin
           rr(jj+nb,L)=rr(nst(jj),L)/
     &     (rrrr*pe(jj+nb)*rrrr*y(nM,l))
          end do              ! monomolecular rates end
        end if
      end do                  !  >>> END TROPOSPHERE LOOP <<<
c
      RETURN
      END SUBROUTINE Crates
c
c
#ifdef TRACERS_SPECIAL_Shindell
      SUBROUTINE checktracer(I,J)
!@sum checktracer for various debugging of tracer chemistry
!@auth Drew Shindell (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on ds3ch4_master_0C4.2_M23)
c
C**** GLOBAL parameters and variables:
C
      USE MODEL_COM, only  : Itime, LM, ls1
      USE TRACER_COM, only : ntm, trname, n_Ox
      USE TRCHEM_Shindell_COM, only: y, nM
c
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
!@var I,J passed horizontal position indicies
!@var L,igas dummy loop variables
!@var tlimit if tracer goes above this limit, model stops
!@var checkOx logical: should I check for large tropospheric Ox?
!@var checkmax logical: should I check for large tracers throughout?
!@var checkNeg logical: should I check for negative tracers?
!@var checkNaN logical: should I check for unreal tracers?
C
      INTEGER L, igas
      INTEGER, INTENT(IN) :: I,J
      REAL*8, DIMENSION(ntm) :: tlimit
      DATA tlimit/9.E-5,1.E-5,1.E-7,3.E-6,1.E-1,1.E-6,3.E-6,1.E-1,1.E-1,
     &1.E-1,1.E-1,1.E-1,1.E-1,1.E-1,1.E-1/
      LOGICAL checkOx, checkmax, checkNeg, checkNan
      DATA checkNeg /.true./
      DATA checkNan /.true./
      DATA checkOx  /.true./
      DATA checkmax /.true./
c
      IF(i.eq.1.and.j.eq.1)
     & WRITE(6,*) 'WARNING: checktracer call is active.'
      IF(checkmax)
     & call stop_model('checktracer: set tlimit for tracers 11->15',255)
C     please (re)set tlimit values for tracers 11 through 15 in the
C     data statement above. Then, please delete the above stop.
C
C check if ozone gets really big in the troposphere:
       IF(checkOx) THEN
       do L=1,LS1-1
         if(y(n_Ox,L)/y(nM,L).gt.1.E-5) then
           write(6,*)'Ox @ I,J,L,Ox,tau:',I,J,L,y(n_Ox,L),Itime
           call stop_model('checktracer: Ox too big in tropo.',255)
         end if
       end do
       END IF
c general check on maximum of tracers:
      IF(checkmax) THEN
      do L=1,LM
       do igas=1,ntm
        if(y(igas,L)/y(nM,L).gt.tlimit(igas)) then
          write(6,*) trname(igas),'@ I,J,L,Ox :',I,J,L,y(igas,L)
          call stop_model('checktracer: tracer upper limit',255)
        end if
       end do
      end do
      END IF
c check for negative tracers:
      IF(checkNeg) THEN
      do L=1,LM
      do igas=1,NTM
       if(y(igas,L).lt.0.) THEN
         write(6,*)trname(igas),
     &   'negative @ tau,I,J,L,y:',Itime,I,J,L,y(igas,L)
         call stop_model('checktracer: tracer is negative',255)
       end if
      enddo
      end do
      END IF
c check for unreal (not-a-number) tracers:
      IF(checkNaN) THEN
      do L=1,LM
      do igas=1,NTM
        if(.NOT.(y(igas,L).gt.0..OR.y(igas,L).le.0.)) THEN
         write(6,*)trname(igas),
     &   'is not a number @ tau,I,J,L,y:',Itime,I,J,L,y(igas,L)
         call stop_model('checktracer: tracer is NaN',255)
        end if
      enddo
      end do
      END IF
      RETURN

      END SUBROUTINE checktracer
#endif
      
      

