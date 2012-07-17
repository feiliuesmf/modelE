#include "rundeck_opts.h"

      subroutine photoj(nslon,nslat)
!@sum from jv_trop.f: FAST J-Value code, troposphere only (mjprather
!@+ 6/96). Uses special wavelength quadrature spectral data
!@+ (jv_spec.dat) that includes only 289 nm - 800 nm (later a single
!@+ 205 nm add-on). Uses special compact Mie code based on
!@+ Feautrier/Auer/Prather version.
!@auth UCI (see note below), GCM incorporation: Drew Shindell,
!@+ modelEifications: Greg Faluvegi
!@calls SET_PROF,JVALUE,PRTATM,JRATET
c  Fastj2 photolysis scheme obtained from H. Bian (UCI) 8/2002.
c  An expanded version of fastJ that includes stratosphere
c  incorporation into the GISS 4x5 GCM with 25 tracer chemistry
c  D. Shindell, Aug. 2002

C**** GLOBAL parameters and variables:

      USE DOMAIN_DECOMP_ATM,only : GRID,getDomainBounds
      USE RESOLUTION, only   : LM
      USE CONSTANT, only     : radian
      USE TRCHEM_Shindell_COM, only: SZA,TFASTJ,JFASTJ,jpnl,jppj,zj,
     &                           szamax,U0,NCFASTJ2,iprn,jprn,prnrts

      IMPLICIT NONE

C**** Local parameters and variables and arguments:
!@var nslon,nslat I and J spatial indicies passed from master chem
!@var i,j,k dummy loop variables
!@var NCFASTJ2 Number of levels in atmosphere
      INTEGER, INTENT(IN) :: nslon, nslat
      INTEGER             :: i,j,k
      logical             :: jay
      INTEGER             :: J_0, J_1 


      call getDomainBounds(grid, J_STRT    =J_0,  J_STOP    =J_1)
      
      jay = (NSLAT >= J_0 .and. NSLAT <= J_1) 
      
      zj(:,:)    =0.d0
      JFASTJ(:,:)=0.d0
      U0 = DCOS(SZA*radian)

      if(SZA <= szamax)then 
        CALL SET_PROF(NSLON,NSLAT)  ! Set up profiles on model levels
        IF(prnrts .and. NSLON == iprn .and. NSLAT == jprn)
     &  CALL PRTATM(2,NSLON,NSLAT,jay) ! Print out atmosphere
        CALL JVALUE(nslon,nslat)    ! Calculate actinic flux
        CALL JRATET(1.d0,NSLAT,NSLON)! Calculate photolysis rates   
        JFASTJ(:,:)= zj(:,:) ! photolysis rates returned to chemistry
      end if
c
      return
      end subroutine photoj



      subroutine set_prof(NSLON,NSLAT)
!@sum set_prof to set up atmospheric profiles required by Fast-J2 using
!@+   a doubled version of the level scheme used in the CTM. First
!@+   pressure and z* altitude are defined, then O3 and T are taken
!@+  from the supplied climatology and integrated to the CTM levels
!@+  (may be overwritten with values directly from the CTM, if desired)
!@+  and then black carbon and aerosol profiles are constructed.
!@+  Oliver Wild (04/07/99)
!@+  Modifications by Apostolos Voulgarakis (Feb 2010) to take aerosol
!@+  tracers from the model into account.
!@auth UCI (see note above), GCM incorporation: Drew Shindell,
!@+ modelEifications: Greg Faluvegi, Apostolos Voulgarakis
c
C**** GLOBAL parameters and variables:
      USE RESOLUTION, only  : JM,LM
      USE GEOM, only: lat2d_dg
      use model_com, only: modelEclock
      USE MODEL_COM, only: Itime
      USE TRCHEM_Shindell_COM, only: TFASTJ,odcol,O3_FASTJ,PFASTJ2,
     &     dlogp,masfac,oref2,tref2,bref2,TJ2,DO32,DBC2,zfastj2,
     &     dmfastj2,NBFASTJ,AER2,MXFASTJ
      USE RAD_COM,only: ttausv_ntrace,ntrix
      USE RADPAR, only : NTRACE
      USE TRACER_COM, only : trname

      IMPLICIT NONE

C**** Local parameters and variables and arguments:
!@var nslon,nslat I and J spatial indicies passed from master chem
!@var pstd Approximate pressures of levels for supplied climatology
!@var skip_tracer logical to not define aer2 for a rad code tracer
      INTEGER, INTENT(IN) :: nslon, nslat
      integer             :: l, k, i, ii, m, j, iclay, n
      real*8, dimension(52) :: pstd
      real*8, dimension(51) :: oref3, tref3
      real*8              :: ydgrd,f0,t0,b0,pb,pc,xc,scaleh
      logical             :: skip_tracer

c  Set up cloud and surface properties
      call CLDSRF(NSLON,NSLAT)

c  Set up pressure levels for O3/T climatology - assume that value
c  given for each 2 km z* level applies from 1 km below to 1 km above,
c  so select pressures at these boundaries. Surface level values at
c  1000 mb are assumed to extend down to the actual P(nslon,nslat).

      pstd(1) = max(PFASTJ2(1),1000.d0)
      pstd(2) = 865.96432336006535d0 !1000.*10.**(-1/16)
      do L=3,51
        pstd(L) = pstd(L-1)*dlogp
      enddo
      pstd(52) = 0.d0

c  Select appropriate monthly and latitudinal profiles:
      ydgrd=lat2d_dg(nslon,nslat)
      m = max(1,min(12,modelEclock%month()))
      l = max(1,min(18,(int(ydgrd)+99)/10))

c  Temporary arrays for climatology data
      oref3(:)=oref2(:,l,m) ! 51
      tref3(:)=tref2(:,l,m) ! 51

c  Apportion O3 and T on supplied climatology z* levels onto CTM levels 
c  with mass (pressure) weighting, assuming constant mixing ratio and
c  temperature half a layer on either side of the point supplied:

      do i = 1,NBFASTJ
        F0 = 0.d0; T0 = 0.d0; B0 = 0.d0
        do k = 1,51
          PC = min(PFASTJ2(i),pstd(k))
          PB = max(PFASTJ2(i+1),pstd(k+1))
          if(PC > PB) then
            XC = (PC-PB)/(PFASTJ2(i)-PFASTJ2(i+1))
            F0 = F0 + oref3(k)*XC
            T0 = T0 + tref3(k)*XC
            B0 = B0 + bref2(k)*XC
          endif
        end do
        TJ2(i) = T0
        DO32(i)= F0*1.d-6
        DBC2(i)= B0
      end do

c Overwrite O3 with GISS chemistry O3:
      DO32(1:LM)=O3_FASTJ(1:LM)
      TJ2(1:LM) =TFASTJ(1:LM)

c  Calculate effective altitudes using scale height at each level
      zfastj2(1) = 0.d0
      do i=1,LM
        scaleh=1.3806d-19*masfac*TFASTJ(i)
        zfastj2(i+1) = zfastj2(i)-(log(PFASTJ2(i+1)/PFASTJ2(i))*scaleh)
      enddo

c  Add Aerosol Column - include aerosol (+cloud) types here. 

      AER2(:,:)=0.d0

c  LAST two are clouds (liquid or ice)
c  Assume limiting temperature for ice of -40 deg C :
      do i=1,LM
        if(TFASTJ(I) > 233.d0) then
          AER2(16,i) = odcol(i)
          AER2(17,i) = 0.d0
        else
          AER2(16,i) = 0.d0
          AER2(17,i) = odcol(i)
        endif
      enddo

#ifndef TRACERS_TOMAS
#ifndef TRACERS_AMP
c Now do the rest of the aerosols
      iclay=0  
      do n=1,NTRACE
        skip_tracer=.false.
        select case (trname(ntrix(n)))
        case ('SO4')      ; j=1 
        case ('seasalt1') ; j=2 
        case ('seasalt2') ; j=3 
        case ('OCIA', 'vbsAm2') ; j=4 ! multiple tracers 
        case ('OCB')      ; j=4 
        case ('isopp1a')  ; j=4 ! multiple tracers
        case ('BCIA')     ; j=5 
        case ('BCB')      ; j=6 
        case ('NO3p')     ; j=7 
        case ('Clay')     ; j=8+iclay ; iclay=iclay+1
        case ('Silt1')    ; j=12
        case ('Silt2')    ; j=13
        case ('Silt3')    ; j=14
        case ('Silt4')    ; j=15
        case default      ; skip_tracer=.true.
        end select
        if(j>15 .or. n>17)
     &  call stop_model("set_prof: too many opt depths",13)
        if(iclay>4)call stop_model("set_prof: too many clays",13)
        if(.not.skip_tracer)then
          do i=1,LM
            if(j.ne.4 .and. AER2(j,i).ne.0.)call stop_model('gsf',13)
            AER2(j,i)=AER2(j,i)+ttausv_ntrace(NSLON,NSLAT,n,i)
          enddo
        endif
      enddo
#endif
#endif

c Top of the atmosphere
      AER2(:,LM+1) = 0.d0

c  Calculate column quantities for Fast-J2:
      do i=1,NBFASTJ
        DMFASTJ2(i)  = (PFASTJ2(i)-PFASTJ2(i+1))*masfac
        DO32(i) = DO32(i)*DMFASTJ2(i)
      enddo
      DO32(NBFASTJ)=DO32(NBFASTJ)*1.E2

      return
      end subroutine set_prof



      SUBROUTINE CLDSRF(NSLON,NSLAT)
!@sum CLDSRF to set cloud and surface properties
!@auth UCI (see note above), GCM incorporation: Drew Shindell,
!@+ modelEifications: Greg Faluvegi

C**** GLOBAL parameters and variables:

      USE RESOLUTION, only   : LM
      USE RAD_COM, only      : ALB
      USE TRCHEM_Shindell_COM, only: RCLOUDFJ,odsum,odmax,
     &            nlbatm,RFLECT,NBFASTJ,AER2,jadsub,dtausub,odcol

      IMPLICIT NONE

C**** Local parameters and variables and arguments:
!@var nslon,nslat I and J spatial indicies passed from master chem
      INTEGER, INTENT(IN) :: nslon, nslat
      integer             :: l, k, j
      real*8              :: odtot

c Default lower photolysis boundary as bottom of level 1
      nlbatm = 1

c Set and limit surface albedo
      RFLECT = max(0.d0,min(1.d0,(1.-ALB(NSLON,NSLAT,1))))

c Zero aerosol column
      AER2(:,:) = 0.d0

c Scale optical depths as appropriate - limit column to 'odmax'
      odsum = 0.d0
      do L=1,LM
        odcol(L) = RCLOUDFJ(L,nslon,nslat)
        odsum = odsum + odcol(L)
      enddo
      if(odsum > odmax) then
        odsum = odmax/odsum
        odcol(:) = odcol(:)*odsum
        odsum = odmax
      endif

c Set sub-division switch if appropriate
      odtot=0.d0
      jadsub(NBFASTJ)=0
      jadsub(NBFASTJ-1)=0
      do L=NBFASTJ-1,1,-1
        k=2*L
        jadsub(k)=0
        jadsub(k-1)=0
        odtot=odtot+odcol(L)
        if(odcol(L) > 0.d0 .and. dtausub > 0.d0) then
          if(odtot <= dtausub) then
            jadsub(k)=1
            jadsub(k-1)=1
          else 
            jadsub(k)=1
            jadsub(k-1)=0
            jadsub(1:2*(L-1))=0
            EXIT
          endif
        endif
      enddo

      return
      end SUBROUTINE CLDSRF



      SUBROUTINE JRATET(SOLF,NSLAT,NSLON)
!@sum JRATET Calculate and print J-values. Note that the loop in
!@+   this routine only covers the jpnl levels actually needed by
!@+   the CTM.
!@auth UCI (see note above), GCM incorporation: Drew Shindell,
!@+ modelEifications: Greg Faluvegi

C**** GLOBAL parameters and variables:

      USE RESOLUTION, only: LM
      USE TRCHEM_Shindell_COM, only: jpnl,TFASTJ,VALJ,NW1,NW2,NJVAL,
     &                               QQQ,JPPJ,ZJ,jfacta,FFF,TQQ,JIND

      IMPLICIT NONE

C**** Local parameters and variables and arguments:
c     FFF    Actinic flux at each level for each wavelength bin
c     QQQ    Cross sections for species (read in in RD_TJPL)
c     SOLF   Solar distance factor, for scaling; normally given by:
c                      1.0-(0.034*cos(real(iday-172)*2.0*pi/365.))
c     TQQ    Temperatures at which QQQ cross sections supplied
      integer :: i, j, k, l, nslon, nslat,jgas
      real*8  :: qo2tot, qo3tot, qo31d, qo33p, qqqt, xseco2, xseco3,
     &           xsec1d, solf, tfact
      real*8, dimension(LM) :: Tx

      Tx(1:LM)=TFASTJ(1:LM)

      DO I=1,jpnl 
        VALJ(1) = 0.d0
        VALJ(2) = 0.d0
        VALJ(3) = 0.d0
        DO K=NW1,NW2
          QO2TOT= XSECO2(K,Tx(I))
          VALJ(1) = VALJ(1) + QO2TOT*FFF(K,I)
          QO3TOT= XSECO3(K,Tx(I))
          QO31D = XSEC1D(K,Tx(I))*QO3TOT
          QO33P = QO3TOT - QO31D
          VALJ(2) = VALJ(2) + QO33P*FFF(K,I)
          VALJ(3) = VALJ(3) + QO31D*FFF(K,I)
        ENDDO
C------ Calculate remaining J-values with T-dep X-sections
        DO J=4,NJVAL !was NJVAL, add -2 for CFC & O2, ds4
          VALJ(J) = 0.d0
          TFACT = 0.d0
          IF(TQQ(2,J) > TQQ(1,J)) TFACT = DMAX1(0.D0,DMIN1(1.D0,
     &    (Tx(I)-TQQ(1,J))/(TQQ(2,J)-TQQ(1,J)) ))
          DO K=NW1,NW2
            QQQT = QQQ(K,1,J-3) + (QQQ(K,2,J-3) - QQQ(K,1,J-3))*TFACT 
            VALJ(J) = VALJ(J) + QQQT*FFF(K,I)
          ENDDO
        ENDDO

        zj(i,1:jppj)=VALJ(jind(1:jppj))*jfacta(1:jppj)*solf
        
      ENDDO

      RETURN
      END SUBROUTINE JRATET



      SUBROUTINE PRTATM(NFASTJq,NSLON,NSLAT,jay)
!@sum PRTATM Print out the atmosphere and calculate appropriate columns
!@auth UCI (see note above), GCM incorporation: Drew Shindell,
!@+ modelEifications: Greg Faluvegi

C**** GLOBAL parameters and variables:
      USE DOMAIN_DECOMP_ATM, only: write_parallel
      USE RESOLUTION, only  : JM
      USE GEOM, only: lat2d_dg
      use model_com, only: modelEclock
      USE TRCHEM_Shindell_COM, only: SZA,NBFASTJ,MXFASTJ,DMFASTJ2,TJ2,
     &             masfac,dlogp2,oref2,tref2,DO32,AER2,PFASTJ2,ZFASTJ2

      IMPLICIT NONE

C**** Local parameters and variables and arguments:
!@var nslon,nslat I and J spatial indicies passed from master chem
!@var NFASTJq Print out 1=column totals only, 2=
!@+   full columns, 3=full columns and climatology
      INTEGER, INTENT(IN) :: nslon, nslat, nfastjq
      INTEGER             :: I, K, M, L
      character(len=300)  :: out_line
      logical             :: jay
      REAL*8, DIMENSION(NBFASTJ)         :: COLO2,COLO3
      REAL*8, DIMENSION(MXFASTJ,NBFASTJ) :: COLAX
      REAL*8, DIMENSION(9)               :: climat
      REAL*8                             :: ZKM,ZSTAR,PJC,ydgrd
     
      if(NFASTJq == 0) return

C---Calculate columns, for diagnostic output only:
      COLO3(NBFASTJ) = DO32(NBFASTJ)
      COLO2(NBFASTJ) = DMFASTJ2(NBFASTJ)*0.20948d0
      COLAX(:,NBFASTJ) = AER2(:,NBFASTJ)
      do I=NBFASTJ-1,1,-1
        COLO3(i) = COLO3(i+1)+DO32(i)
        COLO2(i) = COLO2(i+1)+DMFASTJ2(i)*0.20948d0
        COLAX(:,i) = COLAX(:,i+1)+AER2(:,i)
      enddo
      write(out_line,1200) '  SZA=',sza
      call write_parallel(trim(out_line),crit=jay)
      write(out_line,1200) ' O3-column(DU)=',COLO3(1)/2.687d16,
     &'  column aerosol @1000nm=',(COLAX(K,1),K=1,MXFASTJ)
      call write_parallel(trim(out_line),crit=jay)

C---Print out atmosphere:
      if(NFASTJq > 1) then
        write(out_line,1000) (' AER-X ','col-AER',k=1,mxfastj)
        call write_parallel(trim(out_line),crit=jay)
        do I=NBFASTJ,1,-1
          PJC = PFASTJ2(I)
          ZKM =1.d-5*ZFASTJ2(I)
          ZSTAR = 16.d0*DLOG10(1000.d0/PJC)
          write(out_line,1100) I,ZKM,ZSTAR,DMFASTJ2(I),DO32(I),
     &    1.d6*DO32(I)/DMFASTJ2(I),TJ2(I),PJC,COLO3(I),COLO2(I),
     &    (AER2(K,I),COLAX(K,I),K=1,MXFASTJ)    
          call write_parallel(trim(out_line),crit=jay)
        enddo
      endif            

C---Print out climatology:
      if(NFASTJq > 2) then
        climat(:)=0.d0
        ydgrd=lat2d_dg(nslon,nslat)
        m = max(1,min(12,modelEclock%month()))
        l = max(1,min(18,(int(ydgrd)+99)/10))
        write(out_line,*) 'Specified Climatology'
        call write_parallel(trim(out_line),crit=jay)
        write(out_line,1000)
        call write_parallel(trim(out_line),crit=jay)
        do i=51,1,-1
          PJC = 1000.d0*dlogp2**(2*i-2)
          climat(1) = 16.d0*DLOG10(1000.D0/PJC)
          climat(2) = climat(1)
          climat(3) = PJC*(1.d0/dlogp2-dlogp2)*masfac
          if(i == 1) climat(3)=PJC*(1.d0-dlogp2)*masfac
          climat(4)=climat(3)*oref2(i,l,m)*1.d-6
          climat(5)=oref2(i,l,m)
          climat(6)=tref2(i,l,m)
          climat(7)=PJC
          climat(8)=climat(8)+climat(4)
          climat(9)=climat(9)+climat(3)*0.20948d0
          write(out_line,1100) I,(climat(k),k=1,9)
          call write_parallel(trim(out_line),crit=jay)
        enddo
        write(out_line,1200) ' O3-column(DU)=',climat(8)/2.687d16
        call write_parallel(trim(out_line),crit=jay)
      endif
      
 1000 format(5X,'Zkm',3X,'Z*',8X,'M',8X,'O3',6X,'f-O3',5X,'T',7X,'P',6x,
     &    'col-O3',3X,'col-O2',2X,10(a7,2x))
 1100 format(1X,I2,0P,2F6.2,1P,2E10.3,0P,F7.3,F8.2,F10.4,1P,10E9.2)
 1200 format(A,F8.1,A,10(1pE10.3))
      return
      end SUBROUTINE PRTATM

 

      SUBROUTINE JVALUE(nslon,nslat)
!@sum JVALUE Calculate the actinic flux at each level for the current
!@+   SZA value. 
!@auth UCI (see note above), GCM incorporation: Drew Shindell,
!@+ modelEifications: Greg Faluvegi

C**** GLOBAL parameters and variables:

      USE RESOLUTION, only: LM
      USE TRCHEM_Shindell_COM, only: NW1,NW2,NBFASTJ,WL,FL,FFF,JPNL,TJ2

      IMPLICIT NONE

C**** Local parameters and variables and arguments:
!@var XQO3_2   fastj2 Absorption cross-section of O3
!@var XQO2_2   fastj2 Absorption cross-section of O2
!@var WAVE Effective wavelength of each wavelength bin
!@var AVGF Attenuation of beam at each level for each wavelength
      INTEGER                    :: K,J
      INTEGER, INTENT(IN)        :: NSLON, NSLAT
      REAL*8, DIMENSION(NBFASTJ) :: XQO3_2, XQO2_2
      REAL*8, DIMENSION(JPNL)    :: AVGF
      REAL*8                     :: WAVE
      REAL*8 XSECO2,XSECO3 ! >>> FUNCTIONS <<<

      AVGF(:) = 0.d0   ! JPNL
      FFF(NW1:NW2,:) = 0.d0 ! JPNL
        
C---Calculate spherical weighting functions:
      CALL SPHERE

C---Loop over all wavelength bins:
      DO K=NW1,NW2
        WAVE = WL(K)
        DO J=1,NBFASTJ
          XQO3_2(J) = XSECO3(K,TJ2(J))
          XQO2_2(J) = XSECO2(K,TJ2(J))
        END DO
        CALL OPMIE(K,WAVE,XQO2_2,XQO3_2,AVGF)
        FFF(K,:) = FFF(K,:) + FL(K)*AVGF(:) ! 1,JPNL
      END DO

      RETURN
      END SUBROUTINE JVALUE



      FUNCTION XSECO3(K,TTT)
!@sum XSECO3  O3 Cross-sections for all processes interpolated across
!@+   3 temps
!@auth UCI (see note above), GCM incorporation: Drew Shindell,
!@+ modelEifications: Greg Faluvegi
!@calls FLINT

C**** GLOBAL parameters and variables:

      USE TRCHEM_Shindell_COM, only: TQQ,QO3

      IMPLICIT NONE

C**** Local parameters and variables and arguments:
!@var k passed index for wavelength bin
!@var TTT returned termperature profile
      INTEGER, INTENT(IN) :: k
      real*8              :: TTT
      REAL*8 xseco3,FLINT ! >>> FUNCTIONS <<<
      
      XSECO3  =
     &FLINT(TTT,TQQ(1,2),TQQ(2,2),TQQ(3,2),QO3(K,1),QO3(K,2),QO3(K,3))
      RETURN
      END FUNCTION XSECO3



      FUNCTION XSEC1D(K,TTT)
!@sum XSEC1D  Quantum yields for O3 --> O2 + O(1D) interpolated across
!@+   3 temps
!@auth UCI (see note above), GCM incorporation: Drew Shindell,
!@+ modelEifications: Greg Faluvegi
!@calls FLINT

C**** GLOBAL parameters and variables:

      USE TRCHEM_Shindell_COM, only: TQQ,Q1D

      IMPLICIT NONE

C**** Local parameters and variables and arguments:
!@var k passed index for wavelength bin
!@var TTT returned termperature profile
      INTEGER, INTENT(IN) :: k
      real*8              :: TTT
      REAL*8 xsec1d,FLINT ! >>> FUNCTIONS <<<
      
      XSEC1D =
     &FLINT(TTT,TQQ(1,3),TQQ(2,3),TQQ(3,3),Q1D(K,1),Q1D(K,2),Q1D(K,3))
      RETURN
      END FUNCTION XSEC1D



      FUNCTION XSECO2(K,TTT)
!@sum XSECO2 Cross-sections for O2 interpolated across 3 temps; No
!@+   S_R Bands yet!
!@auth UCI (see note above), GCM incorporation: Drew Shindell,
!@+ modelEifications: Greg Faluvegi
!@calls FLINT

C**** GLOBAL parameters and variables:

      USE TRCHEM_Shindell_COM, only: TQQ,QO2

      IMPLICIT NONE

C**** Local parameters and variables and arguments:
!@var k passed index for wavelength bin
!@var TTT returned termperature profile
      INTEGER, INTENT(IN) :: k
      real*8              :: TTT
      REAL*8 xseco2,FLINT ! >>> FUNCTIONS <<<
      
      XSECO2 =
     &FLINT(TTT,TQQ(1,1),TQQ(2,1),TQQ(3,1),QO2(K,1),QO2(K,2),QO2(K,3))
      RETURN
      END FUNCTION XSECO2



      REAL*8 FUNCTION FLINT(TINT,T1,T2,T3,F1,F2,F3)
!@sum FLINT Three-point linear interpolation function
!@auth UCI (see note above), GCM incorporation: Drew Shindell,
!@+ modelEifications: Greg Faluvegi

      IMPLICIT NONE

C**** Local parameters and variables and arguments:
!@var T1,T2,T3 passed temperature variables
!@var F1,F2,F3 passed X-section variables
!@var TINT returned temperature profile?
      REAL*8, INTENT(IN) :: T1,T2,T3,F1,F2,F3
      real*8             :: TINT

      IF (TINT  <=  T2)  THEN
        IF (TINT  <=  T1)  THEN
          FLINT  = F1
        ELSE
          FLINT = F1 + (F2 - F1)*(TINT -T1)/(T2 -T1)
        ENDIF
      ELSE
        IF (TINT  >=  T3)  THEN
          FLINT  = F3
        ELSE
          FLINT = F2 + (F3 - F2)*(TINT -T2)/(T3 -T2)
        ENDIF
      ENDIF
      RETURN
      END FUNCTION FLINT


      
      SUBROUTINE SPHERE
!@sum SPHERE Calculation of spherical geometry; derive tangent
!@+   heights, slant path lengths and air mass factor for each
!@+   layer. Beyond 90 degrees, include treatment of emergent
!@+   beam (where tangent height is below altitude J-value desired at). 
!@auth UCI (see note above), GCM incorporation: Drew Shindell,
!@+ modelEifications: Greg Faluvegi

C**** GLOBAL parameters and variables:

      USE CONSTANT, only: radius
      USE RESOLUTION, only: LM
      USE TRCHEM_Shindell_COM, only: U0,NBFASTJ,ZFASTJ2,ZZHT,TANHT,
     & nlbatm,AMF

      IMPLICIT NONE

C**** Local parameters and variables and arguments:
!@var AIRMAS Inlined air mass factor function for top of atmosphere
!@var Ux, Htemp dummy arguments to airmas function
!@var GMU MU, cos(solar zenith angle)
!@var RZ Distance from centre of Earth to each point (cm)
!@var RQ Square of radius ratios
!@var XL Slant path between points
      INTEGER :: II, I, J, K
      REAL*8  :: Ux, Htemp, AIRMAS, GMU, ZBYR, xmu1, xmu2, xl, DIFF
      REAL*8, DIMENSION(NBFASTJ) :: RZ, RQ
  
      AIRMAS(Ux,Htemp) = (1.0d0+Htemp)/SQRT(Ux*Ux+2.0d0*Htemp*(1.0d0-
     & 0.6817d0*EXP(-57.3d0*ABS(Ux)/SQRT(1.0d0+5500.d0*Htemp))/
     & (1.0d0+0.625d0*Htemp)))

      GMU = U0
      RZ(1)=radius+ZFASTJ2(1)
      ZBYR = ZZHT/radius
      DO II=2,NBFASTJ
        RZ(II) = radius + ZFASTJ2(II)
        RQ(II-1) = (RZ(II-1)/RZ(II))**2
      END DO
      IF (GMU < 0.d0) THEN
        TANHT = RZ(nlbatm)/DSQRT(1.0d0-GMU**2)
      ELSE
        TANHT = RZ(nlbatm)
      ENDIF

c Go up from the surface calculating the slant paths between each level
c and the level above, and deriving the appropriate Air Mass Factor:
      DO J=1,NBFASTJ
        AMF(:,J)=0.D0 ! K=1,NBFASTJ

c Air Mass Factors all zero if below the tangent height:
        IF (RZ(J) < TANHT) CYCLE

c Ascend from layer J calculating Air Mass Factors (AMFs): 
        XMU1=ABS(GMU)
        DO I=J,LM
          XMU2=DSQRT(1.0d0-RQ(I)*(1.0d0-XMU1**2))
          XL=RZ(I+1)*XMU2-RZ(I)*XMU1
          AMF(I,J)=XL/(RZ(I+1)-RZ(I))
          XMU1=XMU2
        END DO

c Use function and scale height to provide AMF above top of model:
        AMF(NBFASTJ,J)=AIRMAS(XMU1,ZBYR)

c Twilight case - Emergent Beam:
        IF (GMU >= 0.0d0) CYCLE
        XMU1=ABS(GMU)

c Descend from layer J :
        DO II=J-1,1,-1
          DIFF=RZ(II+1)*DSQRT(1.0d0-XMU1**2)-RZ(II)
          if(II == 1) DIFF=max(DIFF,0.d0)   ! filter
c Tangent height below current level - beam passes through twice:
          IF (DIFF < 0.d0) THEN
            XMU2=DSQRT(1.0d0-(1.0d0-XMU1**2)/RQ(II))
            XL=ABS(RZ(II+1)*XMU1-RZ(II)*XMU2)
            AMF(II,J)=2.d0*XL/(RZ(II+1)-RZ(II))
            XMU1=XMU2
c Lowest level intersected by emergent beam;
          ELSE
            XL=RZ(II+1)*XMU1*2.0d0
            AMF(II,J)=XL/(RZ(II+1)-RZ(II))
            CYCLE
          ENDIF
        END DO

      END DO
      RETURN
      END SUBROUTINE SPHERE



      SUBROUTINE OPMIE(KW,WAVEL,XQO2_2,XQO3_2,FMEAN)
!@sum OPMIE NEW Mie code for Js, only uses 8-term expansion, 
!@+   4-Gauss pts.
!@auth UCI (see note above), GCM incorporation: Drew Shindell,
!@+ modelEifications: Greg Faluvegi
c
C Currently allow up to NP aerosol phase functions (at all altitudes)
C to be associated with optical depth AER2(1:NC) = aerosol opt.depth
C @ 1000 nm
C  
C  Pick Mie-wavelength with phase function and Qext:
C
c 01 RAYLE  = Rayleigh phase
c 02 ISOTR  = isotropic
c 03 W_H01 = water haze (H1/Deirm.) (n=1.335, gamma:  r-mode=0.1um / alpha=2)
c 04 W_H04 = water haze (H1/Deirm.) (n=1.335, gamma:  r-mode=0.4um / alpha=2)
c 05 W_C02 = water cloud (C1/Deirm.) (n=1.335, gamma:  r-mode=2.0um / alpha=6)
c 06 W_C04 = water cloud (C1/Deirm.) (n=1.335, gamma:  r-mode=4.0um / alpha=6)
c 07 W_C08 = water cloud (C1/Deirm.) (n=1.335, gamma:  r-mode=8.0um / alpha=6)
c 08 W_C13 = water cloud (C1/Deirm.) (n=1.335, gamma:  r-mode=13.3um / alpha=6)
c 09 W_L06 = water cloud (Lacis) (n=1.335, r-mode=5.5um / alpha=11/3)
c 10 Ice-H = hexagonal ice cloud (Mishchenko)
c 11 Ice-I = irregular ice cloud (Mishchenko)
c 12 SO4   RH=00%  (A. Voulgarakis - Jan 2010)
c 13 SO4   RH=30%  (A. Voulgarakis - Jan 2010)
c 14 SO4   RH=50%  (A. Voulgarakis - Jan 2010)
c 15 SO4   RH=70%  (A. Voulgarakis - Jan 2010)
c 16 SO4   RH=80%  (A. Voulgarakis - Jan 2010)
c 17 SO4   RH=90%  (A. Voulgarakis - Jan 2010)
c 18 SO4   RH=95%  (A. Voulgarakis - Jan 2010)
c 19 SO4   RH=99%  (A. Voulgarakis - Jan 2010)
c 20 SS1   RH=00%  (A. Voulgarakis - Jan 2010)
c 21 SS1   RH=30%  (A. Voulgarakis - Jan 2010)
c 22 SS1   RH=50%  (A. Voulgarakis - Jan 2010)
c 23 SS1   RH=70%  (A. Voulgarakis - Jan 2010)
c 24 SS1   RH=80%  (A. Voulgarakis - Jan 2010)
c 25 SS1   RH=90%  (A. Voulgarakis - Jan 2010)
c 26 SS1   RH=95%  (A. Voulgarakis - Jan 2010)
c 27 SS1   RH=99%  (A. Voulgarakis - Jan 2010)
c 28 SS2   RH=00%  (A. Voulgarakis - Jan 2010)
c 29 SS2   RH=30%  (A. Voulgarakis - Jan 2010)
c 30 SS2   RH=50%  (A. Voulgarakis - Jan 2010)
c 31 SS2   RH=70%  (A. Voulgarakis - Jan 2010)
c 32 SS2   RH=80%  (A. Voulgarakis - Jan 2010)
c 33 SS2   RH=90%  (A. Voulgarakis - Jan 2010)
c 34 SS2   RH=95%  (A. Voulgarakis - Jan 2010)
c 35 SS2   RH=99%  (A. Voulgarakis - Jan 2010)
c 36 O_C   RH=00%  (A. Voulgarakis - Jan 2010)
c 37 O_C   RH=30%  (A. Voulgarakis - Jan 2010)
c 38 O_C   RH=50%  (A. Voulgarakis - Jan 2010)
c 39 O_C   RH=70%  (A. Voulgarakis - Jan 2010)
c 40 O_C   RH=80%  (A. Voulgarakis - Jan 2010)
c 41 O_C   RH=90%  (A. Voulgarakis - Jan 2010)
c 42 O_C   RH=95%  (A. Voulgarakis - Jan 2010)
c 43 O_C   RH=99%  (A. Voulgarakis - Jan 2010)
c 44 B_C  (A. Voulgarakis - Jan 2010)
c 45 NO3   RH=00%  (A. Voulgarakis - Jan 2010)
c 46 NO3   RH=30%  (A. Voulgarakis - Jan 2010)
c 47 NO3   RH=50%  (A. Voulgarakis - Jan 2010)
c 48 NO3   RH=70%  (A. Voulgarakis - Jan 2010)
c 49 NO3   RH=80%  (A. Voulgarakis - Jan 2010)
c 50 NO3   RH=90%  (A. Voulgarakis - Jan 2010)
c 51 NO3   RH=95%  (A. Voulgarakis - Jan 2010)
c 52 NO3   RH=99%  (A. Voulgarakis - Jan 2010)
c 53 DU1  (A. Voulgarakis - Jan 2010)
c 54 DU2  (A. Voulgarakis - Jan 2010)
c 55 DU3  (A. Voulgarakis - Jan 2010)
c 56 DU4  (A. Voulgarakis - Jan 2010)
c 57 DU5  (A. Voulgarakis - Jan 2010)
c 58 DU6  (A. Voulgarakis - Jan 2010)
c 59 DU7  (A. Voulgarakis - Jan 2010)
c 60 DU8  (A. Voulgarakis - Jan 2010)
C
C---------------------------------------------------------------------
C  FUNCTION RAYLAY(WAVE)---RAYLEIGH CROSS-SECTION for wave > 170 nm
C       WSQI = 1.E6/(WAVE*WAVE)
C       REFRM1 = 1.0E-6*(64.328+29498.1/(146.-WSQI)+255.4/(41.-WSQI))
C       RAYLAY = 5.40E-21*(REFRM1*WSQI)**2
C--------------------------------------------------------------------
C
C**** GLOBAL parameters and variables:
C
      USE DOMAIN_DECOMP_ATM, only: write_parallel
      USE RESOLUTION, only : LM
      USE MODEL_COM, only: itime
      USE TRCHEM_Shindell_COM, only: NBFASTJ,POMEGA,NCFASTJ2,
     & POMEGAJ,MIEDX2,QAAFASTJ,SSA,NLBATM,DO32,DMFASTJ2,QRAYL,
     & AMF,PAA,jaddlv,dtaumax,dtausub,dsubdiv,U0,RFLECT,MXFASTJ,
     & NLFASTJ,ZTAU,jadsub,N__,ZU0,ZREFL,ZFLUX,FZ,jaddto,jndlev,
     & FJFASTJ,M__,AER2,MFIT
                          
      IMPLICIT NONE

C**** Local parameters and variables and arguments:
!@var DTAUX Local optical depth of each CTM level
!@var PIRAY2 Contribution of Rayleigh scattering to extinction
!@var PIAER2   Contribution of Aerosol scattering to extinction
!@var TTAU Optical depth of air vertically above each point
!@+   (to top of atm)
!@var FTAU Attenuation of solar beam
!@var FMEAN Mean actinic flux at desired levels
!@var zk fractional increment in level
!@var dttau change in ttau per increment    (linear, positive)
!@var dpomega change in pomega per increment  (linear)
!@var ftaulog change in ftau per increment (exponential, normally < 1)
! POMEGA=Scattering phase function
! jaddlv(i)=Number of new levels to add between (i) and (i+1)
! jaddto(i)=Total number of new levels to add to and above level (i)
! jndlev(j)=Level needed for J-value for CTM layer (j)

      integer :: KW,km,i,j,k,l,ix,j1,ND
      character(len=300) :: out_line
      REAL*8, DIMENSION(NBFASTJ) :: DTAUX,PIRAY2
      REAL*8, INTENT(IN), DIMENSION(NBFASTJ) :: XQO2_2,XQO3_2
      REAL*8, DIMENSION(MXFASTJ,NBFASTJ) :: PIAER2
      REAL*8, DIMENSION(NCFASTJ2+1) :: TTAU,FTAU
      REAL*8, INTENT(OUT), DIMENSION(LM) :: FMEAN
      REAL*8, DIMENSION(MXFASTJ,NBFASTJ) :: QXMIE,SSALB
      REAL*8, DIMENSION(MXFASTJ) :: XLAER
      REAL*8, INTENT(IN) :: WAVEL
      REAL*8, DIMENSION(2*M__) :: dpomega,dpomega2
      REAL*8 xlo2,xlo3,xlray,xltau2,zk,zk2,taudn,tauup,
     & ftaulog,dttau,ftaulog2,dttau2

C---Pick nearest Mie wavelength, no interpolation--------------
                             KM=1
      if( WAVEL  >  355.d0 ) KM=2
      if( WAVEL  >  500.d0 ) KM=3
      if( WAVEL  >  800.d0 ) KM=4 

C---For Mie code scale extinction at 1000 nm to wavelength WAVEL(QXMIE)
      do j=1,NBFASTJ
        QXMIE(1:MXFASTJ,j) =
     &  QAAFASTJ(KM,MIEDX2(j,1:MXFASTJ))/
     &  QAAFASTJ(4,MIEDX2(j,1:MXFASTJ))
        SSALB(1:MXFASTJ,j) = SSA(KM,MIEDX2(j,1:MXFASTJ))
      enddo

C---Reinitialize arrays: ! loop 1,NCFASTJ2+1
      ttau(:)=0.d0
      ftau(:)=0.d0

C---Set up total optical depth over each CTM level, DTAUX:
      J1 = NLBATM
      do J=J1,NBFASTJ
        XLO3=DO32(J)*XQO3_2(J)
        XLO2=DMFASTJ2(J)*XQO2_2(J)*0.20948d0
        XLRAY=DMFASTJ2(J)*QRAYL(KW)
        if(WAVEL <= 291.d0) XLRAY=XLRAY * 0.57d0
        XLAER(:)=AER2(:,J)*QXMIE(:,J) ! MXFASTJ
c Total optical depth from all elements:
        DTAUX(J)=XLO3+XLO2+XLRAY
        do I=1,MXFASTJ
          DTAUX(J)=DTAUX(J)+XLAER(I)
        enddo
c Fractional extinction for Rayleigh scattering and each aerosol type:
        PIRAY2(J)=XLRAY/DTAUX(J)
        PIAER2(:,J)=SSALB(:,J)*XLAER(:)/DTAUX(J) ! MXFASTJ
      enddo ! J

C---Calculate attenuated incident beam EXP(-TTAU/U0) & flux on surface:
      do J=J1,NBFASTJ
        if(AMF(J,J) > 0.0d0) then
          XLTAU2=0.0d0
          do I=1,NBFASTJ
            XLTAU2=XLTAU2 + DTAUX(I)*AMF(I,J)
          enddo
          if(XLTAU2 > 450.d0) then 
            FTAU(j)=0.d0 ! compilers with no underflow trapping
          else
            FTAU(J)=DEXP(-XLTAU2)
          endif
        else
          FTAU(J)=0.0d0
        endif
      enddo

C---in UV region, use pseudo-Rayleigh absorption instead of scattering:
      if (WAVEL <= 291.d0) then
C---Accumulate attenuation for level centers:
        do j=1,LM
          if (j < J1) then
            FMEAN(J) = 0.d0
          else
            FMEAN(J) = sqrt(FTAU(J)*FTAU(J+1))
          endif
        enddo
        return
C---In visible region, consider scattering. Define the scattering
C---phase function with mix of Rayleigh(1) & Mie(MIEDX2).
C No. of quadrature pts fixed at 4 (M__), expansion of phase fn @ 8
      else 
       do j=j1,NBFASTJ
        do i=1,MFIT
         pomegaj(i,j) = PIRAY2(J)*PAA(i,KM,1)
         do k=1,MXFASTJ
          pomegaj(i,j)=pomegaj(i,j)+PIAER2(K,j)*PAA(i,KM,MIEDX2(j,K))
         enddo
        enddo
       enddo
        
C--------------------------------------------------------------------
c  Take optical properties on GCM layers and convert to a photolysis
c  level grid corresponding to layer centres and boundaries. This is
c  required so that J-values can be calculated for the centre of GCM
c  layers; the index of these layers is kept in the jndlev array.
C--------------------------------------------------------------------

c Set lower boundary and levels to calculate J-values at: 
        J1=2*J1-1
        do j=1,LM
          jndlev(j)=2*j
        enddo

c Calculate column optical depths above each level, TTAU:
        TTAU(NCFASTJ2+1)=0.0D0
        do J=NCFASTJ2,J1,-1
          I=(J+1)/2
          TTAU(J)=TTAU(J+1) + 0.5d0*DTAUX(I)
          jaddlv(j)=int(0.5d0*DTAUX(I)/dtaumax)
          ! Subdivide cloud-top levels if required:
          if(jadsub(j) > 0) then
            jadsub(j)=min(jaddlv(j)+1,nint(dtausub))*(nint(dsubdiv)-1)
            jaddlv(j)=jaddlv(j)+jadsub(j)
          endif
        enddo

C---reflect flux from surface

        if(U0 > 0.d0) then
          ZFLUX = U0*FTAU(J1)*RFLECT/(1.d0+RFLECT)
        else
          ZFLUX = 0.d0
        endif

c Calculate attenuated beam, FTAU, level boundaries then level centres:
        FTAU(NCFASTJ2+1)=1.0d0
        do J=NCFASTJ2-1,J1,-2
          I=(J+1)/2
          FTAU(J)=FTAU(I)
        enddo
        do J=NCFASTJ2,J1,-2
          FTAU(J)=sqrt(FTAU(J+1)*FTAU(J-1))
        enddo

c Calculate scattering properties, level centres then level boundaries,
c using an inverse interpolation to give correctly-weighted values:
        do j=NCFASTJ2,J1,-2
          pomegaj(1:MFIT,j) = pomegaj(1:MFIT,j/2)
        enddo
        do j=J1+2,NCFASTJ2,2
          taudn = ttau(j-1)-ttau(j)
          tauup = ttau(j)-ttau(j+1)
          pomegaj(1:MFIT,j) = (pomegaj(1:MFIT,j-1)*taudn + 
     &    pomegaj(1:MFIT,j+1)*tauup) / (taudn+tauup)
        enddo
c Define lower and upper boundaries:
        pomegaj(1:MFIT,J1) = pomegaj(1:MFIT,J1+1)
        pomegaj(1:MFIT,NCFASTJ2+1) = pomegaj(1:MFIT,NCFASTJ2)

C--------------------------------------------------------------------
c  Calculate cumulative total and define levels at which we want
c  the J-values.  Sum upwards for levels, and then downwards for Mie
c  code readjustments.
C--------------------------------------------------------------------

c Reinitialize level arrays:
        jaddto(:)=0 ! NLFASTJ+1

        jaddto(J1)=jaddlv(J1)
        do j=J1+1,NCFASTJ2
          jaddto(j)=jaddto(j-1)+jaddlv(j)
        enddo
        if((jaddto(NCFASTJ2)+NCFASTJ2) > NLFASTJ) then
          write(out_line,1500)
     &    jaddto(NCFASTJ2)+NCFASTJ2,'NLFASTJ',NLFASTJ
          call write_parallel(trim(out_line),crit=.true.)
          call stop_model('problem in fastj2 with jaddto; This may
     &    reflect issues with the optical depth of a tracer or 
     &    clouds being unreasonable. Please check!!
     &    If in 1997 (very anomalous fire year), it is OK to 
     &    increase NLFASTJ to 1200 in TRCHEM_Shindell_COM.
     &    But for other years, please check
     &    for possible issues.',255)
        endif
        jndlev(:)=jndlev(:)+jaddto(jndlev(:)-1) ! LM
        jaddto(NCFASTJ2)=jaddlv(NCFASTJ2)
        do j=NCFASTJ2-1,J1,-1
          jaddto(j)=jaddto(j+1)+jaddlv(j)
        enddo

C---------------------SET UP FOR MIE CODE--------------------------
c
c  Transpose the ascending TTAU grid to a descending ZTAU grid.
c  Double the resolution - TTAU points become the odd points on the
c  ZTAU grid, even points needed for asymm phase fn soln, contain 'h'.
c  Odd points added at top of grid for unattenuated beam (Z='inf')
c  
c        Surface:   TTAU(1)   now use ZTAU(2*NCFASTJ2+1)
c        Top:       TTAU(NCFASTJ2)  now use ZTAU(3)
c        Infinity:            now use ZTAU(1)
c
c  Mie scattering code only used from surface to level NCFASTJ2
C---------------------------------------------------------------------

c Initialise all Fast-J2 optical property arrays:
        pomega(:,:) = 0.d0 ! (2*M__,N__)
        ztau(:)     = 0.d0 ! N__
        fz(:)       = 0.d0 ! N__

c Ascend through atmosphere transposing grid and adding extra points:
        do j=J1,NCFASTJ2+1
          k = 2*(NCFASTJ2+1-j)+2*jaddto(j)+1
          ztau(k)= ttau(j)
          fz(k)  = ftau(j)
          pomega(:,k) = pomegaj(:,j) ! MFIT
        enddo

C---------------------------------------------------------------------
c  Insert new levels, working downwards from the top of the atmosphere
c  to the surface (down in 'j', up in 'k'). This allows ztau and pomega
c  to be incremented linearly (in a +ve sense), and the flux fz to be
c  attenuated top-down (avoiding problems where lower level fluxes are
c  zero).
C---------------------------------------------------------------------

        do j=NCFASTJ2,J1,-1
          zk = 0.5d0/(1.d0+dble(jaddlv(j)-jadsub(j)))
          dttau = (ttau(j)-ttau(j+1))*zk
          dpomega(:) = (pomegaj(:,j)-pomegaj(:,j+1))*zk ! MFIT
c  Filter attenuation factor - set minimum at 1.0d-05
          if(ftau(j+1) == 0.d0) then
            ftaulog=0.d0
          else
            ftaulog = ftau(j)/ftau(j+1)
            if(ftaulog < 1.d-150) then
              ftaulog=1.0d-05
            else
              ftaulog=exp(log(ftaulog)*zk)
            endif
          endif
          k = 2*(NCFASTJ2-j+jaddto(j)-jaddlv(j))+1   !  k at level j+1
          l = 0          
          
c Additional subdivision of first level if required:
          if(jadsub(j) /= 0) then
            l=jadsub(j)/nint(dsubdiv-1)
            zk2=1.d0/dsubdiv
            dttau2=dttau*zk2
            ftaulog2=ftaulog**zk2
            dpomega2(:)=dpomega(:)*zk2 ! MFIT
            do ix=1,2*(jadsub(j)+l)
              ztau(k+1) = ztau(k) + dttau2
              fz(k+1) = fz(k)*ftaulog2
              pomega(:,k+1) = pomega(:,k) + dpomega2(:) ! MFIT
              k = k+1
              if(k > 1800)then
                write(out_line,*) 'k fault:',k,NCFASTJ2,j,jaddto(j),
     &          jadsub(j),dsubdiv,jaddlv(j)
                call write_parallel(trim(out_line),crit=.true.)
              endif
            enddo
          endif
          l = 2*(jaddlv(j)-jadsub(j)-l)+1

c Add values at all intermediate levels:
          do ix=1,l
            ztau(k+1) = ztau(k) + dttau
            fz(k+1) = fz(k)*ftaulog
            pomega(:,k+1) = pomega(:,k) + dpomega(:) ! MFIT
            k = k+1
          enddo
        enddo

C---Update total number of levels and check does not exceed N__
        ND = 2*(NCFASTJ2+jaddto(J1)-J1)  + 3
        if(nd > N__) then
          write(out_line,1500) ND, 'N__',N__
          call write_parallel(trim(out_line),crit=.true.)
          call stop_model('problem in fastj2 with ND; This may
     &    reflect issues with the optical depth of a tracer or 
     &    clouds being unreasonable. Please check!!
     &    If in 1997 (very anomalous fire year), it is OK to 
     &    increase N__ to 3000 in TRCHEM_Shindell_COM.
     &    But for other years, please check
     &    for possible issues.',255)
        endif

C---Add boundary/ground layer to ensure no negative Js caused by
C---too large a TTAU-step in the 2nd-order lower b.c.
        ZTAU(ND+1) = ZTAU(ND)*1.000005d0
        ZTAU(ND+2) = ZTAU(ND)*1.000010d0
        zk=max(abs(U0),0.01d0)
        zk=dexp(-ZTAU(ND)*5.d-6/zk)
        FZ(ND+1) = FZ(ND)*zk
        FZ(ND+2) = FZ(ND+1)*zk
        POMEGA(:,ND+1)   = POMEGA(:,ND) ! MFIT
        POMEGA(:,ND+2)   = POMEGA(:,ND) ! MFIT
        ND = ND+2
        ZU0 = U0
        ZREFL = RFLECT

C-----------------------------------------
        CALL MIESCT(ND)
C-----------------------------------------

c Accumulate attenuation for selected levels:
        l=2*(NCFASTJ2+jaddto(J1))+3
        do j=1,LM
          k=l-(2*jndlev(j))
          if(k > ND-2) then
            FMEAN(j) = 0.d0
          else
            FMEAN(j) = FJFASTJ(k)
          endif
        enddo

      endif ! WAVEL

      return
 1000 format(1x,i3,3(2x,1pe10.4),1x,i3)
 1300 format(1x,50(i3))
 1500 format(' Too many levels in photolysis code: need ',i5,' but ',a,
     $       ' dimensioned as ',i4)
      END SUBROUTINE OPMIE      



      SUBROUTINE MIESCT(ND)
!@sum MIESCT This is an adaption of the Prather rad transfer code. 
!@+  see comments. 
!@auth UCI (see note above), GCM incorporation: Drew Shindell,
!@+ modelEifications: Greg Faluvegi
!@calls BLKSLV, GAUSSP, LEGND0
C
C-------------------------------------------------------------------
C   This is an adaption of the Prather rad transfer code, (mjp, 10/95)
C     Prather, 1974, Astrophys. J. 192, 787-792.
C         Soln of inhomogeneous Rayleigh scattering atmosphere.
C         (original Rayleigh w/ polarization)
C     Cochran and Trafton, 1978, Ap.J., 219, 756-762.
C         Raman scattering in the atmospheres of the major planets.
C         (first use of anisotropic code)
C     Jacob, Gottlieb and Prather,89, J.G..Res., 94, 12975-13002.
C         Chemistry of a polluted cloudy boundary layer,
C         (documentation of extension to anisotropic scattering)
C
C    takes atmospheric structure and source terms from std J-code
C    ALSO limited to 4 Gauss points, only calculates mean field!
C
C   mean rad. field ONLY (M=1)
C   initialize variables FIXED/UNUSED in this special version:
C   FTOP=1.0=astrophys flux (unit of pi) at SZA, -ZU0, use for scaling
C   FBOT=0.0=ext isotropic flux on lower boundary

C**** GLOBAL parameters and variables:

      USE TRCHEM_Shindell_COM, only: NFASTJ,EMU,WTFASTJ,ZFLUX,ZU0,
     & MFASTJ,MFIT,ZREFL,FZ,PM0,PM,FJFASTJ
                          
      IMPLICIT NONE

C**** Local parameters and variables and arguments:

      integer           :: I, id, imm, ND
      real*8, parameter :: cmeq1 = 0.25D0

C Fix scattering to 4 Gausss pts = 8-stream.
C Solve eqn of R.T. only for first-order M=1
      ZFLUX = (ZU0*FZ(ND)*ZREFL)/(1.0d0+ZREFL)
      DO I=1,NFASTJ
        CALL LEGND0 (EMU(I),PM0,MFIT)
        PM(I,MFASTJ:MFIT) = PM0(MFASTJ:MFIT)
      ENDDO

      CALL LEGND0 (-ZU0,PM0,MFIT)
      PM0(MFASTJ:MFIT) = CMEQ1*PM0(MFASTJ:MFIT)

      CALL BLKSLV(ND)

      DO ID=1,ND,2
        FJFASTJ(ID) = 4.0d0*FJFASTJ(ID) + FZ(ID)
      ENDDO
      RETURN
      END SUBROUTINE MIESCT



      SUBROUTINE BLKSLV(ND)
!@sum BLKSLV Solves the block tri-diagonal system:
!@+   A(I)*X(I-1) + B(I)*X(I) + C(I)*X(I+1) = H(I)
!@auth UCI (see note above), GCM incorporation: Drew Shindell,
!@+ modelEifications: Greg Faluvegi

C**** GLOBAL parameters and variables:

      USE TRCHEM_Shindell_COM, only: BFASTJ,NFASTJ,CC,HFASTJ,
     & AFASTJ,DD,RR2,C1,FJFASTJ,WTFASTJ,EMU,AAFASTJ
                           
      IMPLICIT NONE

C**** Local parameters and variables and arguments:
      integer :: i, j, k, id, ND
      real*8  :: sum

C-----------UPPER BOUNDARY ID=1
      CALL GEN(1,ND)
      CALL MATIN4(BFASTJ)
      DO I=1,NFASTJ
        RR2(I,1) = 0.0d0
        DO J=1,NFASTJ
          SUM = 0.0d0
          DO K=1,NFASTJ
            SUM = SUM - BFASTJ(I,K)*CC(K,J)
          ENDDO
          DD(I,J,1) = SUM
          RR2(I,1) = RR2(I,1) + BFASTJ(I,J)*HFASTJ(J)
        ENDDO
      ENDDO
C----------CONTINUE THROUGH ALL DEPTH POINTS ID=2 TO ID=ND-1
      DO ID=2,ND-1
        CALL GEN(ID,ND)
        DO I=1,NFASTJ
          DO J=1,NFASTJ
            BFASTJ(I,J) = BFASTJ(I,J) + AFASTJ(I)*DD(I,J,ID-1)
          ENDDO
          HFASTJ(I) = HFASTJ(I) - AFASTJ(I)*RR2(I,ID-1)
        ENDDO
        CALL MATIN4 (BFASTJ)
        DO I=1,NFASTJ
          RR2(I,ID) = 0.0d0
          DO J=1,NFASTJ
            RR2(I,ID) = RR2(I,ID) + BFASTJ(I,J)*HFASTJ(J)
            DD(I,J,ID) = - BFASTJ(I,J)*C1(J)
          ENDDO
        ENDDO
      ENDDO
C---------FINAL DEPTH POINT: ND
      CALL GEN(ND,ND)
      DO I=1,NFASTJ
        DO J=1,NFASTJ
          SUM = 0.0d0
          DO K=1,NFASTJ
            SUM = SUM + AAFASTJ(I,K)*DD(K,J,ND-1)
          ENDDO
          BFASTJ(I,J) = BFASTJ(I,J) + SUM
          HFASTJ(I) = HFASTJ(I) - AAFASTJ(I,J)*RR2(J,ND-1)
        ENDDO
      ENDDO
      CALL MATIN4 (BFASTJ)
      DO I=1,NFASTJ
        RR2(I,ND) = 0.0d0
        DO J=1,NFASTJ        
          RR2(I,ND) = RR2(I,ND) + BFASTJ(I,J)*HFASTJ(J)
        ENDDO
      ENDDO
C-----------BACK SOLUTION
      DO ID=ND-1,1,-1
        DO I=1,NFASTJ
          DO J=1,NFASTJ
            RR2(I,ID) = RR2(I,ID) + DD(I,J,ID)*RR2(J,ID+1)
          ENDDO
        ENDDO
      ENDDO
C----------MEAN J & H
      DO ID=1,ND,2
        FJFASTJ(ID) = 0.0d0
        DO I=1,NFASTJ
          FJFASTJ(ID) = FJFASTJ(ID) + RR2(I,ID)*WTFASTJ(I)
        ENDDO
      ENDDO
      DO ID=2,ND,2
        FJFASTJ(ID) = 0.0d0
        DO I=1,NFASTJ
          FJFASTJ(ID) = FJFASTJ(ID) + RR2(I,ID)*WTFASTJ(I)*EMU(I)
        ENDDO
      ENDDO
      RETURN
      END SUBROUTINE BLKSLV



      SUBROUTINE GEN(ID,ND)
!@sum GEN Generates coefficient matrices for the block tri-diagonal
!@+    system:  A(I)*X(I-1) + B(I)*X(I) + C(I)*X(I+1) = H(I)
!@auth UCI (see note above), GCM incorporation: Drew Shindell,
!@+ modelEifications: Greg Faluvegi

C**** GLOBAL parameters and variables:

      USE TRCHEM_Shindell_COM, only: NFASTJ,MFASTJ,MFIT,POMEGA,PM,PM0,
     & HFASTJ,FZ,SFASTJ,WTFASTJ,WFASTJ,U1,V1,BFASTJ,CC,AFASTJ,EMU,C1,
     & ZTAU,AAFASTJ,ZFLUX,ZREFL
                           
      IMPLICIT NONE

C**** Local parameters and variables and arguments:
      integer :: id, id0, id1, im, i, j, k, mstart, ND
      real*8  :: sum0, sum1, sum2, sum3, deltau, d1, d2, surfac
      REAL*8, DIMENSION(8) ::  TTMP

C---------------------------------------------
      IF(ID == 1 .OR. ID == ND) THEN
C---------calculate generic 2nd-order terms for boundaries
       ID0 = ID
       ID1 = ID+1
       IF(ID >= ND) ID1 = ID-1
       DO I=1,NFASTJ
         SUM0 = 0.0d0
         SUM1 = 0.0d0
         SUM2 = 0.0d0
         SUM3 = 0.0d0
         DO IM=MFASTJ,MFIT,2
           SUM0 = SUM0 + POMEGA(IM,ID0)*PM(I,IM)*PM0(IM)
           SUM2 = SUM2 + POMEGA(IM,ID1)*PM(I,IM)*PM0(IM)
         ENDDO
         DO IM=MFASTJ+1,MFIT,2
           SUM1 = SUM1 + POMEGA(IM,ID0)*PM(I,IM)*PM0(IM)
           SUM3 = SUM3 + POMEGA(IM,ID1)*PM(I,IM)*PM0(IM)
         ENDDO
         HFASTJ(I) = 0.5d0*(SUM0*FZ(ID0) + SUM2*FZ(ID1))
         AFASTJ(I) = 0.5d0*(SUM1*FZ(ID0) + SUM3*FZ(ID1))
         DO J=1,I
           SUM0 = 0.0d0
           SUM1 = 0.0d0
           SUM2 = 0.0d0
           SUM3 = 0.0d0
           DO IM=MFASTJ,MFIT,2
             SUM0 = SUM0 + POMEGA(IM,ID0)*PM(I,IM)*PM(J,IM)
             SUM2 = SUM2 + POMEGA(IM,ID1)*PM(I,IM)*PM(J,IM)
           ENDDO
           DO IM=MFASTJ+1,MFIT,2
             SUM1 = SUM1 + POMEGA(IM,ID0)*PM(I,IM)*PM(J,IM)
             SUM3 = SUM3 + POMEGA(IM,ID1)*PM(I,IM)*PM(J,IM)
           ENDDO
           SFASTJ(I,J) = - SUM2*WTFASTJ(J)
           SFASTJ(J,I) = - SUM2*WTFASTJ(I)
           WFASTJ(I,J) = - SUM1*WTFASTJ(J)
           WFASTJ(J,I) = - SUM1*WTFASTJ(I)
           U1(I,J) = - SUM3*WTFASTJ(J)
           U1(J,I) = - SUM3*WTFASTJ(I)
           SUM0 = 0.5d0*(SUM0 + SUM2)
           BFASTJ(I,J) = - SUM0*WTFASTJ(J)
           BFASTJ(J,I) = - SUM0*WTFASTJ(I)
         ENDDO
         SFASTJ(I,I) = SFASTJ(I,I) + 1.0d0
         WFASTJ(I,I) = WFASTJ(I,I) + 1.0d0
         U1(I,I) = U1(I,I) + 1.0d0
         BFASTJ(I,I) = BFASTJ(I,I) + 1.0d0
       END DO ! I
      
       DO I=1,NFASTJ
         SUM0 = 0.0d0
         DO J=1,NFASTJ
           SUM0 = SUM0 + SFASTJ(I,J)*AFASTJ(J)/EMU(J)
         ENDDO
         C1(I) = SUM0
       ENDDO
       DO I=1,NFASTJ
         DO J=1,NFASTJ
           SUM0 = 0.0d0
           SUM2 = 0.0d0
           DO K=1,NFASTJ
             SUM0 = SUM0 + SFASTJ(J,K)*WFASTJ(K,I)/EMU(K)
             SUM2 = SUM2 + SFASTJ(J,K)*U1(K,I)/EMU(K)
           ENDDO
           AFASTJ(J) = SUM0
           V1(J) = SUM2
         ENDDO
         DO J=1,NFASTJ
           WFASTJ(J,I) = AFASTJ(J)
           U1(J,I) = V1(J)
         ENDDO
       ENDDO

       IF (ID == 1) THEN
C-------------upper boundary, 2nd-order, C-matrix is full (CC)
        DELTAU = ZTAU(2) - ZTAU(1)
        D2 = 0.25d0*DELTAU
        DO I=1,NFASTJ
          D1 = EMU(I)/DELTAU
          DO J=1,NFASTJ
            BFASTJ(I,J) = BFASTJ(I,J) + D2*WFASTJ(I,J)
            CC(I,J) = D2*U1(I,J)
          ENDDO
          BFASTJ(I,I) = BFASTJ(I,I) + D1
          CC(I,I) = CC(I,I) - D1
          HFASTJ(I) = HFASTJ(I) + 2.0d0*D2*C1(I)
          AFASTJ(I) = 0.0d0
        ENDDO
       ELSE
C-------------lower boundary, 2nd-order, A-matrix is full (AAFASTJ)
        DELTAU = ZTAU(ND) - ZTAU(ND-1)
        D2 = 0.25d0*DELTAU
        SURFAC = 4.0d0*ZREFL/(1.0d0 + ZREFL)
        DO I=1,NFASTJ
          D1 = EMU(I)/DELTAU
          HFASTJ(I) = HFASTJ(I) - 2.0d0*D2*C1(I)
          SUM0 = 0.0d0
          DO J=1,NFASTJ
            SUM0 = SUM0 + WFASTJ(I,J)
          ENDDO
          SUM0 = D1 + D2*SUM0
          SUM1 = SURFAC*SUM0
          DO J=1,NFASTJ
            BFASTJ(I,J)=BFASTJ(I,J) + 
     &      D2*WFASTJ(I,J)-SUM1*EMU(J)*WTFASTJ(J)
          ENDDO
          BFASTJ(I,I) = BFASTJ(I,I) + D1
          HFASTJ(I) = HFASTJ(I) + SUM0*ZFLUX
          DO J=1,NFASTJ
            AAFASTJ(I,J) = - D2*U1(I,J)
          ENDDO
          AAFASTJ(I,I) = AAFASTJ(I,I) + D1
          C1(I) = 0.0d0
        ENDDO
       ENDIF

C------------intermediate points:  can be even or odd, A & C diagonal
      ELSE

        DELTAU = ZTAU(ID+1) - ZTAU(ID-1)
        MSTART = MFASTJ + MOD(ID+1,2)
        DO I=1,NFASTJ
          AFASTJ(I) = EMU(I)/DELTAU
          C1(I) = -AFASTJ(I)
          SUM0 = 0.0d0
          DO IM=MSTART,MFIT,2
            TTMP(IM) = POMEGA(IM,ID)*PM(I,IM) 
            SUM0 = SUM0 + TTMP(IM)*PM0(IM)  
          ENDDO
          HFASTJ(I) = SUM0*FZ(ID)
          DO J=1,I
            SUM0 = 0.0d0
            DO IM=MSTART,MFIT,2
              SUM0 = SUM0 + TTMP(IM)*PM(J,IM)  
            ENDDO
            BFASTJ(I,J) =  - SUM0*WTFASTJ(J)
            BFASTJ(J,I) =  - SUM0*WTFASTJ(I)
          ENDDO
          BFASTJ(I,I) = BFASTJ(I,I) + 1.0d0
        ENDDO
      ENDIF

      RETURN
      END SUBROUTINE GEN



      SUBROUTINE LEGND0(X,PL,NFASTJ)
!@sum LEGND0 Calculates ORDINARY LEGENDRE fns of X (real)
!@+   from P[0] = PL(1) = 1,  P[1] = X, .... P[N-1] = PL(N)
!@auth UCI (see note above), GCM incorporation: Drew Shindell,
!@+ modelEifications: Greg Faluvegi
!@calls 
                           
      IMPLICIT NONE

C**** Local parameters and variables and arguments:

      INTEGER :: NFASTJ,I
      REAL*8  :: X,PL(NFASTJ),DEN
      
C---Always does PL(2) = P[1]
      PL(1) = 1.D0
      PL(2) = X
      DO I=3,NFASTJ
        DEN = (I-1)
        PL(I) = PL(I-1)*X*(2.d0-1.D0/DEN) - PL(I-2)*(1.d0-1.D0/DEN)
      ENDDO
      RETURN
      END SUBROUTINE LEGND0



      SUBROUTINE MATIN4(AFASTJ)
!@sum MATIN4 invert 4x4 matrix A(4,4) in place with L-U decomp
!@+   (mjp, old...)
!@auth UCI (see note above), GCM incorporation: Drew Shindell,
!@+ modelEifications: Greg Faluvegi

      IMPLICIT NONE

C**** Local parameters and variables and arguments:
!@var AFASTJ passed (actually BFASTJ)
      REAL*8 :: AFASTJ(4,4)

C---SETUP L AND U
      AFASTJ(2,1) = AFASTJ(2,1)/AFASTJ(1,1)
      AFASTJ(2,2) = AFASTJ(2,2)-AFASTJ(2,1)*AFASTJ(1,2)
      AFASTJ(2,3) = AFASTJ(2,3)-AFASTJ(2,1)*AFASTJ(1,3)
      AFASTJ(2,4) = AFASTJ(2,4)-AFASTJ(2,1)*AFASTJ(1,4)
      AFASTJ(3,1) = AFASTJ(3,1)/AFASTJ(1,1)
      AFASTJ(3,2) = (AFASTJ(3,2)-AFASTJ(3,1)*AFASTJ(1,2))/AFASTJ(2,2)
      AFASTJ(3,3) = AFASTJ(3,3)-AFASTJ(3,1)*AFASTJ(1,3)-AFASTJ(3,2)
     $*AFASTJ(2,3)
      AFASTJ(3,4) = AFASTJ(3,4)-AFASTJ(3,1)*AFASTJ(1,4)-AFASTJ(3,2)
     $*AFASTJ(2,4)
      AFASTJ(4,1) = AFASTJ(4,1)/AFASTJ(1,1)
      AFASTJ(4,2) = (AFASTJ(4,2)-AFASTJ(4,1)*AFASTJ(1,2))/AFASTJ(2,2)
      AFASTJ(4,3) = (AFASTJ(4,3)-AFASTJ(4,1)*AFASTJ(1,3)-AFASTJ(4,2)
     $*AFASTJ(2,3))/AFASTJ(3,3)
      AFASTJ(4,4) = AFASTJ(4,4)-AFASTJ(4,1)*AFASTJ(1,4)-AFASTJ(4,2)*
     $AFASTJ(2,4)-AFASTJ(4,3)*AFASTJ(3,4)
C---INVERT L
      AFASTJ(4,3) = -AFASTJ(4,3)
      AFASTJ(4,2) = -AFASTJ(4,2)-AFASTJ(4,3)*AFASTJ(3,2)
      AFASTJ(4,1) = -AFASTJ(4,1)-AFASTJ(4,2)*AFASTJ(2,1)-AFASTJ(4,3)
     $*AFASTJ(3,1)
      AFASTJ(3,2) = -AFASTJ(3,2)
      AFASTJ(3,1) = -AFASTJ(3,1)-AFASTJ(3,2)*AFASTJ(2,1)
      AFASTJ(2,1) = -AFASTJ(2,1)
C---INVERT U
      AFASTJ(4,4) = 1.D0/AFASTJ(4,4)
      AFASTJ(3,4) = -AFASTJ(3,4)*AFASTJ(4,4)/AFASTJ(3,3)
      AFASTJ(3,3) = 1.D0/AFASTJ(3,3)
      AFASTJ(2,4) = -(AFASTJ(2,3)*AFASTJ(3,4)+AFASTJ(2,4)*
     $AFASTJ(4,4))/AFASTJ(2,2)
      AFASTJ(2,3) = -AFASTJ(2,3)*AFASTJ(3,3)/AFASTJ(2,2)
      AFASTJ(2,2) = 1.D0/AFASTJ(2,2)
      AFASTJ(1,4) = -(AFASTJ(1,2)*AFASTJ(2,4)+AFASTJ(1,3)*AFASTJ(3,4)
     $+AFASTJ(1,4)*AFASTJ(4,4))/AFASTJ(1,1)
      AFASTJ(1,3) = -(AFASTJ(1,2)*AFASTJ(2,3)+AFASTJ(1,3)*
     $AFASTJ(3,3))/AFASTJ(1,1)
      AFASTJ(1,2) = -AFASTJ(1,2)*AFASTJ(2,2)/AFASTJ(1,1)
      AFASTJ(1,1) = 1.D0/AFASTJ(1,1)
C---MULTIPLY (U-INVERSE)*(L-INVERSE)
      AFASTJ(1,1) = AFASTJ(1,1)+AFASTJ(1,2)*AFASTJ(2,1)+AFASTJ(1,3)
     $*AFASTJ(3,1)+AFASTJ(1,4)*AFASTJ(4,1)
      AFASTJ(1,2) = AFASTJ(1,2)+AFASTJ(1,3)*AFASTJ(3,2)+AFASTJ(1,4)
     $*AFASTJ(4,2)
      AFASTJ(1,3) = AFASTJ(1,3)+AFASTJ(1,4)*AFASTJ(4,3)
      AFASTJ(2,1) = AFASTJ(2,2)*AFASTJ(2,1)+AFASTJ(2,3)*AFASTJ(3,1)
     $+AFASTJ(2,4)*AFASTJ(4,1)
      AFASTJ(2,2) = AFASTJ(2,2)+AFASTJ(2,3)*AFASTJ(3,2)+AFASTJ(2,4)
     $*AFASTJ(4,2)
      AFASTJ(2,3) = AFASTJ(2,3)+AFASTJ(2,4)*AFASTJ(4,3)
      AFASTJ(3,1) = AFASTJ(3,3)*AFASTJ(3,1)+AFASTJ(3,4)*AFASTJ(4,1)
      AFASTJ(3,2) = AFASTJ(3,3)*AFASTJ(3,2)+AFASTJ(3,4)*AFASTJ(4,2)
      AFASTJ(3,3) = AFASTJ(3,3)+AFASTJ(3,4)*AFASTJ(4,3)
      AFASTJ(4,1) = AFASTJ(4,4)*AFASTJ(4,1)
      AFASTJ(4,2) = AFASTJ(4,4)*AFASTJ(4,2)
      AFASTJ(4,3) = AFASTJ(4,4)*AFASTJ(4,3)
      RETURN
      END SUBROUTINE MATIN4
