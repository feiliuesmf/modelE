#include "rundeck_opts.h"
c
      subroutine photoj(nslon,nslat)
!@sum from jv_trop.f: FAST J-Value code, troposphere only (mjprather
!@+ 6/96). Uses special wavelength quadrature spectral data
!@+ (jv_spec.dat) that includes only 289 nm - 800 nm (later a single
!@+ 205 nm add-on). Uses special compact Mie code based on
!@+ Feautrier/Auer/Prather version.
!@auth UCI (see note below), GCM incorporation: Drew Shindell,
!@+ modelEifications: Greg Faluvegi
!@ver  1.0 (based on ds4p_fastj2_M23.S)
!@calls SET_PROF,JVALUE,PRTATM,JRATET
c  Fastj2 photolysis scheme obtained from H. Bian (UCI) 8/2002.
c  An expanded version of fastJ that includes stratosphere
c  incorporation into the GISS 4x5 GCM with 25 tracer chemistry
c  D. Shindell, Aug. 2002
c
C**** GLOBAL parameters and variables:
C
      USE MODEL_COM, only    : LM
      USE CONSTANT, only     : radian
      USE TRCHEM_Shindell_COM, only: SZA,TFASTJ,JFASTJ,jpnl,jppj,zj,
     &                            szamax,U0,NCFASTJ2,iprn,jprn,prnrts
c
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
!@var nslon,nslat I and J spatial indicies passed from master chem
!@var i,j,k dummy loop variables
!@var NCFASTJ2 Number of levels in atmosphere
      INTEGER, INTENT(IN) :: nslon, nslat
      INTEGER i,j,k
C
#ifndef SHINDELL_STRAT_CHEM
      call stop_model(
     &'only use TRCHEM_fastj2 ifdef SHINDELL_STRAT_CHEM',255)
#endif
C
      zj(:,:)    =0.d0
      JFASTJ(:,:)=0.d0
      U0 = DCOS(SZA*radian)
c
      if(SZA.le.szamax)then 
        CALL SET_PROF(NSLON,NSLAT)  ! Set up profiles on model levels
        IF(prnrts.and.NSLON.eq.iprn.and.NSLAT.eq.jprn)
     &  CALL PRTATM(2,NSLON,NSLAT)  ! Print out atmosphere
        call sys_flush(6)
        CALL JVALUE(nslon,nslat)    ! Calculate actinic flux
        CALL JRATET(1.d0,NSLAT,NSLON)! Calculate photolysis rates   
        JFASTJ(:,:)= zj(:,:) ! photolysis rates returned to chemistry
      end if
c
      return
      end subroutine photoj
C
C
      subroutine set_prof(NSLON,NSLAT)
!@sum set_prof to set up atmospheric profiles required by Fast-J2 using
!@+   a doubled version of the level scheme used in the CTM. First
!@+   pressure and z* altitude are defined, then O3 and T are taken
!@+  from the supplied climatology and integrated to the CTM levels
!@+  (may be overwritten with values directly from the CTM, if desired)
!@+  and then black carbon and aerosol profiles are constructed.
!@+  Oliver Wild (04/07/99)
!@auth UCI (see note above), GCM incorporation: Drew Shindell,
!@+ modelEifications: Greg Faluvegi
!@ver  1.0 (based on ds4p_fastj2_M23)
c
C**** GLOBAL parameters and variables:
      USE MODEL_COM, only: IM,JM,LM,month=>JMON
      USE TRCHEM_Shindell_COM, only: TFASTJ,odcol,O3_FASTJ,PFASTJ2,
     &     dlogp,masfac,oref2,tref2,bref2,TJ2,DO32,DBC2,zfastj2,
     &     dmfastj2,NBFASTJ,AER2
c
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
!@var nslon,nslat I and J spatial indicies passed from master chem
!@var pstd Approximate pressures of levels for supplied climatology
      INTEGER, INTENT(IN) :: nslon, nslat
      integer l, k, i, m
      real*8 pstd(52),ydgrd,oref3(51),tref3(51),f0,t0,b0,pb,pc,
     & xc,scaleh
c
c  Set up cloud and surface properties
      call CLDSRF(NSLON,NSLAT)
c
c  Set up pressure levels for O3/T climatology - assume that value
c  given for each 2 km z* level applies from 1 km below to 1 km above,
c  so select pressures at these boundaries. Surface level values at
c  1000 mb are assumed to extend down to the actual P(nslon,nslat).
c
      pstd(1) = max(PFASTJ2(1),1000.d0)
      pstd(2) = 865.96432336006535d0 !1000.*10.**(-1/16)
      do L=3,51
        pstd(L) = pstd(L-1)*dlogp
      enddo
      pstd(52) = 0.d0
c
c  Select appropriate monthly and latitudinal profiles
      ydgrd=-(90.d0-(real(nslat)*90.d0/(real(JM)/2.d0)))!lat in degrees
      m = max(1,min(12,month))
      l = max(1,min(18,(int(ydgrd)+99)/10))
c
c  Temporary arrays for climatology data
      oref3(:)=oref2(:,l,m) ! 51
      tref3(:)=tref2(:,l,m) ! 51
c      
c  Apportion O3 and T on supplied climatology z* levels onto CTM levels 
c  with mass (pressure) weighting, assuming constant mixing ratio and
c  temperature half a layer on either side of the point supplied.
c
      do i = 1,NBFASTJ
        F0 = 0.d0; T0 = 0.d0; B0 = 0.d0
        do k = 1,51
          PC = min(PFASTJ2(i),pstd(k))
          PB = max(PFASTJ2(i+1),pstd(k+1))
          if(PC.gt.PB) then
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
c
c Overwrite O3 with GISS GCM O3
      do i=1,LM
       DO32(i)=O3_FASTJ(i)
       TJ2(i)=TFASTJ(i)
      enddo
c
c  Calculate effective altitudes using scale height at each level
      zfastj2(1) = 0.d0
      do i=1,LM
        scaleh=1.3806d-19*masfac*TFASTJ(i)
        zfastj2(i+1) = zfastj2(i)-(log(PFASTJ2(i+1)/PFASTJ2(i))*scaleh)
      enddo
c
c  Add Aerosol Column - include aerosol types here. Currently use soot
c  water and ice; assume black carbon x-section of 10 m2/g, independent
c  of wavelength; assume limiting temperature for ice of -40 deg C.
c
      do i=1,LM
        AER2(1,i) = DBC2(i)*10.d0*(zfastj2(i+1)-zfastj2(i))
        if(TFASTJ(I).gt.233.d0) then
          AER2(2,i) = odcol(i)
          AER2(3,i) = 0.d0
        else
          AER2(2,i) = 0.d0
          AER2(3,i) = odcol(i)
        endif          
      enddo
      AER2(:,LM+1) = 0.d0
c
c  Calculate column quantities for Fast-J2
      do i=1,NBFASTJ
        DMFASTJ2(i)  = (PFASTJ2(i)-PFASTJ2(i+1))*masfac
        DO32(i) = DO32(i)*DMFASTJ2(i)
      enddo
      DO32(NBFASTJ)=DO32(NBFASTJ)*1.E2
c
      return
      end subroutine set_prof
c
c
      SUBROUTINE CLDSRF(NSLON,NSLAT)
!@sum CLDSRF to set cloud and surface properties
!@auth UCI (see note above), GCM incorporation: Drew Shindell,
!@+ modelEifications: Greg Faluvegi
!@ver  1.0 (based on ds4p_fastj2_M23)
c
C**** GLOBAL parameters and variables:
C
      USE MODEL_COM, only    : IM,JM,LM
      USE TRCHEM_Shindell_COM, only: SALBFJ,RCLOUDFJ,odsum,odmax,
     &nlbatm,RFLECT,NBFASTJ,AER2,jadsub,dtausub,odcol
c
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
!@var nslon,nslat I and J spatial indicies passed from master chem
      INTEGER, INTENT(IN) :: nslon, nslat
      integer l, k, j
      real*8 odtot
c
c Default lower photolysis boundary as bottom of level 1
      nlbatm = 1
c
c Set and limit surface albedo
      RFLECT = max(0.d0,min(1.d0,(1.-SALBFJ(NSLON,NSLAT))))
c
c Zero aerosol column
      AER2(:,:) = 0.d0
c
c Scale optical depths as appropriate - limit column to 'odmax'
      odsum = 0.d0
      do L=1,LM
        odcol(L) = RCLOUDFJ(nslon,nslat,L)
        odsum = odsum + odcol(L)
      enddo
      if(odsum.gt.odmax) then
        odsum = odmax/odsum
        odcol(:) = odcol(:)*odsum
        odsum = odmax
      endif
c
c Set sub-division switch if appropriate
      odtot=0.d0
      jadsub(NBFASTJ)=0
      jadsub(NBFASTJ-1)=0
      do L=NBFASTJ-1,1,-1
        k=2*L
        jadsub(k)=0
        jadsub(k-1)=0
        odtot=odtot+odcol(L)
        if(odcol(L).gt.0.d0.and.dtausub.gt.0.d0) then
          if(odtot.le.dtausub) then
            jadsub(k)=1
            jadsub(k-1)=1
          else if(odtot.gt.dtausub) then
            jadsub(k)=1
            jadsub(k-1)=0
            do j=1,2*(L-1)
              jadsub(j)=0
            enddo
            go to 20
          endif
        endif
      enddo
 20   continue
c
      return
      end SUBROUTINE CLDSRF
C
C
      SUBROUTINE JRATET(SOLF,NSLAT,NSLON)
!@sum JRATET Calculate and print J-values. Note that the loop in
!@+   this routine only covers the jpnl levels actually needed by
!@+   the CTM.
!@auth UCI (see note above), GCM incorporation: Drew Shindell,
!@+ modelEifications: Greg Faluvegi
!@ver  1.0 (based on ds4p_fastj2_M23)
c
C**** GLOBAL parameters and variables:
C
      USE MODEL_COM, only: IM,JM,LM
      USE TRCHEM_Shindell_COM, only: jpnl,TFASTJ,VALJ,NW1,NW2,NJVAL,
     &   QQQ,JPPJ,ZJ,jfacta,FFF,TQQ,JIND
c
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:

c     FFF    Actinic flux at each level for each wavelength bin
c     QQQ    Cross sections for species (read in in RD_TJPL)
c     SOLF   Solar distance factor, for scaling; normally given by:
c                      1.0-(0.034*cos(real(iday-172)*2.0*pi/365.))
c     TQQ    Temperatures at which QQQ cross sections supplied

      integer i, j, k, l, nslon, nslat,jgas
      real*8 qo2tot, qo3tot, qo31d, qo33p, qqqt
      real*8 xseco2, xseco3, xsec1d, solf, tfact
      real*8, dimension(IM,JM,LM) :: T

      T(nslon,nslat,:)=TFASTJ(:) ! LM

      DO I=1,jpnl ! how can this line crash the model.
       VALJ(1) = 0.d0
       VALJ(2) = 0.d0
       VALJ(3) = 0.d0
       DO K=NW1,NW2
         QO2TOT= XSECO2(K,dble(T(nslon,nslat,I)))
         VALJ(1) = VALJ(1) + QO2TOT*FFF(K,I)
         QO3TOT= XSECO3(K,dble(T(nslon,nslat,I)))
         QO31D = XSEC1D(K,dble(T(nslon,nslat,I)))*QO3TOT
         QO33P = QO3TOT - QO31D
         VALJ(2) = VALJ(2) + QO33P*FFF(K,I)
         VALJ(3) = VALJ(3) + QO31D*FFF(K,I)
       ENDDO
C
C------Calculate remaining J-values with T-dep X-sections
       DO J=4,NJVAL !was NJVAL, add -2 for CFC & O2, ds4
         VALJ(J) = 0.d0
         TFACT = 0.d0
         IF(TQQ(2,J).GT.TQQ(1,J)) TFACT = DMAX1(0.D0,DMIN1(1.D0,
     $      (T(nslon,nslat,I)-TQQ(1,J))/(TQQ(2,J)-TQQ(1,J)) ))
         DO K=NW1,NW2
           QQQT = QQQ(K,1,J-3) + (QQQ(K,2,J-3) - QQQ(K,1,J-3))*TFACT 
           VALJ(J) = VALJ(J) + QQQT*FFF(K,I)
         ENDDO
       ENDDO
C
       DO jgas=1,jppj  !was NJVAL, add -2 for CFC & O2, ds4        
         zj(i,jgas)=VALJ(jind(jgas))*jfacta(jgas)*solf
       ENDDO
      ENDDO

      RETURN
      END SUBROUTINE JRATET
c
c
      SUBROUTINE PRTATM(NFASTJq,NSLON,NSLAT)
!@sum PRTATM Print out the atmosphere and calculate appropriate columns
!@auth UCI (see note above), GCM incorporation: Drew Shindell,
!@+ modelEifications: Greg Faluvegi
!@ver  1.0 (based on ds4p_fastj2_M23)
c
C**** GLOBAL parameters and variables:
C
      USE MODEL_COM, only: JM, month=>JMON 
      USE TRCHEM_Shindell_COM, only: SZA,NBFASTJ,MXFASTJ,DMFASTJ2,TJ2,
     & masfac,dlogp2,oref2,tref2,DO32,AER2,PFASTJ2,ZFASTJ2
c
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
!@var nslon,nslat I and J spatial indicies passed from master chem
!@var NFASTJq Print out 1=column totals only, 2=
!@+   full columns, 3=full columns and climatology
      INTEGER, INTENT(IN) :: nslon, nslat, nfastjq
      INTEGER I, K, M, L
      REAL*8, DIMENSION(NBFASTJ)         :: COLO2,COLO3
      REAL*8, DIMENSION(MXFASTJ,NBFASTJ) :: COLAX
      REAL*8, DIMENSION(9)               :: climat
      REAL*8 ZKM,ZSTAR,PJC,ydgrd
c     
      if(NFASTJq.eq.0) return
c
C---Calculate columns, for diagnostic output only:
      COLO3(NBFASTJ) = DO32(NBFASTJ)
      COLO2(NBFASTJ) = DMFASTJ2(NBFASTJ)*0.20948d0
      COLAX(:,NBFASTJ) = AER2(:,NBFASTJ)
      do I=NBFASTJ-1,1,-1
        COLO3(i) = COLO3(i+1)+DO32(i)
        COLO2(i) = COLO2(i+1)+DMFASTJ2(i)*0.20948d0
        COLAX(:,i) = COLAX(:,i+1)+AER2(:,i)
      enddo
      write(6,1200) '  SZA=',sza
      write(6,1200) ' O3-column(DU)=',COLO3(1)/2.687d16,
     $            '  column aerosol @1000nm=',(COLAX(K,1),K=1,MXFASTJ)
C
C---Print out atmosphere
      if(NFASTJq.gt.1) then
        write(6,1000) (' AER-X ','col-AER',k=1,mxfastj)
        do I=NBFASTJ,1,-1
          PJC = PFASTJ2(I)
          ZKM =1.d-5*ZFASTJ2(I)
          ZSTAR = 16.d0*DLOG10(1000.d0/PJC)
          write(6,1100) I,ZKM,ZSTAR,DMFASTJ2(I),DO32(I),
     $    1.d6*DO32(I)/DMFASTJ2(I),TJ2(I),PJC,COLO3(I),COLO2(I),
     $    (AER2(K,I),COLAX(K,I),K=1,MXFASTJ)    
        enddo
      endif            
c
C---Print out climatology
      if(NFASTJq.gt.2) then
        climat(:)=0.d0
        ydgrd=-(90.d0-(real(nslat)*90.d0/(real(JM)/2.d0)))!lat degrees
        m = max(1,min(12,month))
        l = max(1,min(18,(int(ydgrd)+99)/10))
        write(6,*) 'Specified Climatology'
        write(6,1000)
        do i=51,1,-1
          PJC = 1000.d0*dlogp2**(2*i-2)
          climat(1) = 16.d0*DLOG10(1000.D0/PJC)
          climat(2) = climat(1)
          climat(3) = PJC*(1.d0/dlogp2-dlogp2)*masfac
          if(i.eq.1) climat(3)=PJC*(1.d0-dlogp2)*masfac
          climat(4)=climat(3)*oref2(i,l,m)*1.d-6
          climat(5)=oref2(i,l,m)
          climat(6)=tref2(i,l,m)
          climat(7)=PJC
          climat(8)=climat(8)+climat(4)
          climat(9)=climat(9)+climat(3)*0.20948d0
          write(6,1100) I,(climat(k),k=1,9)
        enddo
        write(6,1200) ' O3-column(DU)=',climat(8)/2.687d16
      endif
      return
c
 1000 format(5X,'Zkm',3X,'Z*',8X,'M',8X,'O3',6X,'f-O3',5X,'T',7X,'P',6x,
     $    'col-O3',3X,'col-O2',2X,10(a7,2x))
 1100 format(1X,I2,0P,2F6.2,1P,2E10.3,0P,F7.3,F8.2,F10.4,1P,10E9.2)
 1200 format(A,F8.1,A,10(1pE10.3))
      end SUBROUTINE PRTATM
C      
C
      SUBROUTINE JVALUE(nslon,nslat)
!@sum JVALUE Calculate the actinic flux at each level for the current
!@+   SZA value. 
!@auth UCI (see note above), GCM incorporation: Drew Shindell,
!@+ modelEifications: Greg Faluvegi
!@ver  1.0 (based on ds4p_fastj2_M23)
c
C**** GLOBAL parameters and variables:
C
      USE MODEL_COM, only: LM
      USE TRCHEM_Shindell_COM, only: NW1,NW2,NBFASTJ,WL,FL,FFF,JPNL,TJ2
c
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
!@var XQO3_2   fastj2 Absorption cross-section of O3
!@var XQO2_2   fastj2 Absorption cross-section of O2
!@var WAVE Effective wavelength of each wavelength bin
!@var AVGF Attenuation of beam at each level for each wavelength
      INTEGER K
      INTEGER, INTENT(IN) :: NSLON, NSLAT
      REAL*8, DIMENSION(NBFASTJ) :: XQO3_2, XQO2_2
      REAL*8, DIMENSION(JPNL)    :: AVGF
      REAL*8 WAVE
      REAL*8 XSECO2,XSECO3 ! >>> FUNCTIONS <<<
C
      DO K=NW1,NW2
       FFF(K,:) = 0.d0 ! J=1,JPNL
      ENDDO
C        
C---Calculate spherical weighting functions
      CALL SPHERE
C
C---Loop over all wavelength bins
      DO K=NW1,NW2
        WAVE = WL(K)
        XQO3_2(:) = XSECO3(K,TJ2(:)) ! J=1,NBFASTJ
        XQO2_2(:) = XSECO2(K,TJ2(:)) ! J=1,NBFASTJ
        CALL OPMIE(K,WAVE,XQO2_2,XQO3_2,AVGF)
        FFF(K,:) = FFF(K,:) + FL(K)*AVGF(:) ! J=1,JPNL
      ENDDO
c
      RETURN
      END SUBROUTINE JVALUE
C
C
      FUNCTION XSECO3(K,TTT)
!@sum XSECO3  O3 Cross-sections for all processes interpolated across
!@+   3 temps
!@auth UCI (see note above), GCM incorporation: Drew Shindell,
!@+ modelEifications: Greg Faluvegi
!@ver  1.0 (based on ds4p_fastj2_M23)
!@calls FLINT
c
C**** GLOBAL parameters and variables:
C
      USE TRCHEM_Shindell_COM, only: TQQ,QO3
c
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
!@var k passed index for wavelength bin
!@var TTT returned termperature profile
      INTEGER, INTENT(IN) :: k
      REAL*8 xseco3,FLINT ! >>> FUNCTIONS <<<
      real*8 TTT
c
      XSECO3  =
     F  FLINT(TTT,TQQ(1,2),TQQ(2,2),TQQ(3,2),QO3(K,1),QO3(K,2),QO3(K,3))
      RETURN
      END FUNCTION XSECO3
C
C
      FUNCTION XSEC1D(K,TTT)
!@sum XSEC1D  Quantum yields for O3 --> O2 + O(1D) interpolated across
!@+   3 temps
!@auth UCI (see note above), GCM incorporation: Drew Shindell,
!@+ modelEifications: Greg Faluvegi
!@ver  1.0 (based on ds4p_fastj2_M23)
!@calls FLINT
c
C**** GLOBAL parameters and variables:
C
      USE TRCHEM_Shindell_COM, only: TQQ,Q1D
c
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
!@var k passed index for wavelength bin
!@var TTT returned termperature profile
      INTEGER, INTENT(IN) :: k
      REAL*8 xsec1d,FLINT ! >>> FUNCTIONS <<<
      real*8 TTT
c
      XSEC1D =
     F  FLINT(TTT,TQQ(1,3),TQQ(2,3),TQQ(3,3),Q1D(K,1),Q1D(K,2),Q1D(K,3))
      RETURN
      END FUNCTION XSEC1D
C
C
      FUNCTION XSECO2(K,TTT)
!@sum XSECO2 Cross-sections for O2 interpolated across 3 temps; No
!@+   S_R Bands yet!
!@auth UCI (see note above), GCM incorporation: Drew Shindell,
!@+ modelEifications: Greg Faluvegi
!@ver  1.0 (based on ds4p_fastj2_M23)
!@calls FLINT
c
C**** GLOBAL parameters and variables:
C
      USE TRCHEM_Shindell_COM, only: TQQ,QO2
c
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
!@var k passed index for wavelength bin
!@var TTT returned termperature profile
      INTEGER, INTENT(IN) :: k
      REAL*8 xseco2,FLINT ! >>> FUNCTIONS <<<
      real*8 TTT
c
      XSECO2 =
     & FLINT(TTT,TQQ(1,1),TQQ(2,1),TQQ(3,1),QO2(K,1),QO2(K,2),QO2(K,3))
      RETURN
      END FUNCTION XSECO2
c
c
      REAL*8 FUNCTION FLINT(TINT,T1,T2,T3,F1,F2,F3)
!@sum FLINT Three-point linear interpolation function
!@auth UCI (see note above), GCM incorporation: Drew Shindell,
!@+ modelEifications: Greg Faluvegi
!@ver  1.0 (based on ds4p_fastj2_M23)
c
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
!@var T1,T2,T3 passed temperature variables
!@var F1,F2,F3 passed X-section variables
!@var TINT returned temperature profile?
      REAL*8, INTENT(IN) :: T1,T2,T3,F1,F2,F3
      real*8 TINT
!!    REAL*8 FLINT ! >>> FUNCTIONS <<<

      IF (TINT .LE. T2)  THEN
        IF (TINT .LE. T1)  THEN
          FLINT  = F1
        ELSE
          FLINT = F1 + (F2 - F1)*(TINT -T1)/(T2 -T1)
        ENDIF
      ELSE
        IF (TINT .GE. T3)  THEN
          FLINT  = F3
        ELSE
          FLINT = F2 + (F3 - F2)*(TINT -T2)/(T3 -T2)
        ENDIF
      ENDIF
      RETURN
      END FUNCTION FLINT
c
c      
      SUBROUTINE SPHERE
!@sum SPHERE Calculation of spherical geometry; derive tangent
!@+   heights, slant path lengths and air mass factor for each
!@+   layer. Beyond 90 degrees, include treatment of emergent
!@+   beam (where tangent height is below altitude J-value desired at). 
!@auth UCI (see note above), GCM incorporation: Drew Shindell,
!@+ modelEifications: Greg Faluvegi
!@ver  1.0 (based on ds4p_fastj2_M23)
c
C**** GLOBAL parameters and variables:
C
      USE CONSTANT, only: radius
      USE MODEL_COM, only: LM
      USE TRCHEM_Shindell_COM, only: U0,NBFASTJ,ZFASTJ2,ZZHT,TANHT,
     & nlbatm,AMF
c
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
!@var AIRMAS Inlined air mass factor function for top of atmosphere
!@var Ux, Htemp dummy arguments to airmas function
!@var GMU MU, cos(solar zenith angle)
!@var RZ Distance from centre of Earth to each point (cm)
!@var RQ Square of radius ratios
!@var XL Slant path between points
      REAL*8 Ux, Htemp, AIRMAS, GMU, ZBYR, xmu1, xmu2, xl, DIFF
      REAL*8, DIMENSION(NBFASTJ) :: RZ, RQ
      INTEGER II, I, J, K
c  
      AIRMAS(Ux,Htemp) = (1.0d0+Htemp)/SQRT(Ux*Ux+2.0d0*Htemp*(1.0d0-
     $         0.6817d0*EXP(-57.3d0*ABS(Ux)/SQRT(1.0d0+5500.d0*Htemp))/
     $                                         (1.0d0+0.625d0*Htemp)))
c
      GMU = U0
      RZ(1)=radius+ZFASTJ2(1)
      ZBYR = ZZHT/radius
      DO II=2,NBFASTJ
        RZ(II) = radius + ZFASTJ2(II)
        RQ(II-1) = (RZ(II-1)/RZ(II))**2
      END DO
      IF (GMU.LT.0.d0) THEN
        TANHT = RZ(nlbatm)/DSQRT(1.0d0-GMU**2)
      ELSE
        TANHT = RZ(nlbatm)
      ENDIF
c
c  Go up from the surface calculating the slant paths between each level
c  and the level above, and deriving the appropriate Air Mass Factor
      DO J=1,NBFASTJ
        AMF(:,J)=0.D0 ! K=1,NBFASTJ
c
c  Air Mass Factors all zero if below the tangent height
        IF (RZ(J).LT.TANHT) CYCLE
C
c  Ascend from layer J calculating AMFs
        XMU1=ABS(GMU)
        DO I=J,LM
          XMU2=DSQRT(1.0d0-RQ(I)*(1.0d0-XMU1**2))
          XL=RZ(I+1)*XMU2-RZ(I)*XMU1
          AMF(I,J)=XL/(RZ(I+1)-RZ(I))
          XMU1=XMU2
        END DO
c
c  Use function and scale height to provide AMF above top of model
        AMF(NBFASTJ,J)=AIRMAS(XMU1,ZBYR)
c
c  Twilight case - Emergent Beam
        IF (GMU.GE.0.0d0) CYCLE
        XMU1=ABS(GMU)
c
c  Descend from layer J
        DO II=J-1,1,-1
          DIFF=RZ(II+1)*DSQRT(1.0d0-XMU1**2)-RZ(II)
          if(II.eq.1) DIFF=max(DIFF,0.d0)   ! filter
c  Tangent height below current level - beam passes through twice
          IF (DIFF.LT.0.d0) THEN
            XMU2=DSQRT(1.0d0-(1.0d0-XMU1**2)/RQ(II))
            XL=ABS(RZ(II+1)*XMU1-RZ(II)*XMU2)
            AMF(II,J)=2.d0*XL/(RZ(II+1)-RZ(II))
            XMU1=XMU2
c  Lowest level intersected by emergent beam
          ELSE
            XL=RZ(II+1)*XMU1*2.0d0
            AMF(II,J)=XL/(RZ(II+1)-RZ(II))
            CYCLE
          ENDIF
        END DO
c
      END DO
      RETURN
      END SUBROUTINE SPHERE
C
C
      SUBROUTINE OPMIE(KW,WAVEL,XQO2_2,XQO3_2,FMEAN)
!@sum OPMIE NEW Mie code for Js, only uses 8-term expansion, 
!@+   4-Gauss pts.
!@auth UCI (see note above), GCM incorporation: Drew Shindell,
!@+ modelEifications: Greg Faluvegi
!@ver  1.0 (based on ds4p_fastj2_M23)
c
C Currently allow up to NP aerosol phase functions (at all altitudes)
C to be associated with optical depth AER2(1:NC) = aerosol opt.depth
C @ 1000 nm
C  
C  Pick Mie-wavelength with phase function and Qext:
C
C  01 RAYLE = Rayleigh phase
C  02 ISOTR = isotropic
C  03 ABSRB = fully absorbing 'soot', wavelength indep.
C  04 S_Bkg = backgrnd stratospheric sulfate
C             (n=1.46,log-norm:r=.09um/sigma=.6)
C  05 S_Vol = volcanic stratospheric sulfate
C             (n=1.46,log-norm:r=.08um/sigma=.8)
C  06 W_H01 = water haze (H1/Deirm.)
C             (n=1.335, gamma:  r-mode=0.1um /alpha=2)
C  07 W_H04 = water haze (H1/Deirm.)
C             (n=1.335, gamma:  r-mode=0.4um /alpha=2)
C  08 W_C02 = water cloud (C1/Deirm.)
C             (n=1.335, gamma:  r-mode=2.0um /alpha=6)
C  09 W_C04 = water cloud (C1/Deirm.)
C             (n=1.335, gamma:  r-mode=4.0um /alpha=6)
C  10 W_C08 = water cloud (C1/Deirm.)
C             (n=1.335, gamma:  r-mode=8.0um /alpha=6)
C  11 W_C13 = water cloud (C1/Deirm.)
C             (n=1.335, gamma:  r-mode=13.3um /alpha=6)
C  12 W_L06 = water cloud (Lacis) (n=1.335, r-mode=5.5um / alpha=11/3)
C  13 Ice-H = hexagonal ice cloud (Mishchenko)
C  14 Ice-I = irregular ice cloud (Mishchenko)
C
C  Choice of aerosol index MIEDX2 is made in TRCHEM_Shindell_COM.f
C  Optical depths are apportioned to the AER2 array in SET_PROF
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
      USE MODEL_COM, only: LM
      USE TRCHEM_Shindell_COM, only: NBFASTJ,POMEGA,NCFASTJ2,
     & POMEGAJ,MIEDX2,QAAFASTJ,SSA,NLBATM,DO32,DMFASTJ2,QRAYL,
     & AMF,PAA,jaddlv,dtaumax,dtausub,dsubdiv,U0,RFLECT,MXFASTJ,
     & NLFASTJ,ZTAU,jadsub,N__,ZU0,ZREFL,ZFLUX,FZ,jaddto,jndlev,
     & FJFASTJ,M__,AER2,MFIT
c                           
      IMPLICIT NONE
c
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
C
      integer KW,km,i,j,k,l,ix,j1,ND
      REAL*8, DIMENSION(NBFASTJ) :: DTAUX,PIRAY2
      REAL*8, INTENT(IN), DIMENSION(NBFASTJ) :: XQO2_2,XQO3_2
      REAL*8, DIMENSION(MXFASTJ,NBFASTJ) :: PIAER2
      REAL*8, DIMENSION(NCFASTJ2+1) :: TTAU,FTAU
      REAL*8, INTENT(OUT), DIMENSION(LM) :: FMEAN
      REAL*8, DIMENSION(MXFASTJ) :: QXMIE,XLAER,SSALB
      REAL*8, INTENT(IN) :: WAVEL
      REAL*8, DIMENSION(2*M__) :: dpomega,dpomega2
      REAL*8 xlo2,xlo3,xlray,xltau2,zk,zk2,taudn,tauup,
     & ftaulog,dttau,ftaulog2,dttau2
c
C---Pick nearest Mie wavelength, no interpolation--------------
                              KM=1
      if( WAVEL .gt. 355.d0 ) KM=2
      if( WAVEL .gt. 500.d0 ) KM=3
C     if( WAVEL .gt. 800.d0 ) KM=4  !drop the 1000 nm wavelength
c
C---For Mie code scale extinction at 1000 nm to wavelength
C---WAVEL (QXMIE) 
      do k=1,MXFASTJ
        QXMIE(k) = QAAFASTJ(KM,MIEDX2(k))/QAAFASTJ(4,MIEDX2(k))
        SSALB(k) = SSA(KM,MIEDX2(k))
      end do
c
C---Reinitialize arrays ! loop 1,NCFASTJ2+1
        ttau(:)=0.d0
        ftau(:)=0.d0
c
C---Set up total optical depth over each CTM level, DTAUX
      J1 = NLBATM
      do J=J1,NBFASTJ
        XLO3=DO32(J)*XQO3_2(J)
        XLO2=DMFASTJ2(J)*XQO2_2(J)*0.20948d0
        XLRAY=DMFASTJ2(J)*QRAYL(KW)
        if(WAVEL.le.291.d0) XLRAY=XLRAY * 0.57d0
c  Zero absorption for testing purposes
         XLAER(:)=AER2(:,J)*QXMIE(:) ! MXFASTJ
c  Total optical depth from all elements
        DTAUX(J)=XLO3+XLO2+XLRAY
        do I=1,MXFASTJ
          DTAUX(J)=DTAUX(J)+XLAER(I)
        enddo
c  Fractional extinction for Rayleigh scattering and each aerosol type
        PIRAY2(J)=XLRAY/DTAUX(J)
        PIAER2(:,J)=SSALB(:)*XLAER(:)/DTAUX(J) ! MXFASTJ
      enddo ! J
c
C---Calculate attenuated incident beam EXP(-TTAU/U0) and flux on
C---surface
      do J=J1,NBFASTJ
        if(AMF(J,J).gt.0.0d0) then
          XLTAU2=0.0d0
          do I=1,NBFASTJ
            XLTAU2=XLTAU2 + DTAUX(I)*AMF(I,J)
          enddo
          if(XLTAU2.gt.450.d0) then 
            FTAU(j)=0.d0 ! compilers with no underflow trapping
          else
            FTAU(J)=DEXP(-XLTAU2)
          endif
        else
          FTAU(J)=0.0d0
        endif
      enddo
c
C---in UV region, use pseudo-Rayleigh absorption instead of scattering
      if (WAVEL.le.291.d0) then
c
C---Accumulate attenuation for level centers
c
      do j=1,LM
        if (j.lt.J1) then
          FMEAN(J) = 0.d0
        else
          FMEAN(J) = sqrt(FTAU(J)*FTAU(J+1))
        endif
      enddo
      return
c
C---in visible region, consider scattering
C---Define the scattering phase fn. with mix of Rayleigh(1) &
C---Mie(MIEDX2)
C No. of quadrature pts fixed at 4 (M__), expansion of phase fn @ 8
c
      else
c      
      do j=j1,NBFASTJ
        do i=1,MFIT
          pomegaj(i,j) = PIRAY2(J)*PAA(i,KM,1)
          do k=1,MXFASTJ
            pomegaj(i,j)=pomegaj(i,j)+PIAER2(K,j)*PAA(i,KM,MIEDX2(K))
          enddo
        enddo
      enddo
c
C--------------------------------------------------------------------
c  Take optical properties on CTM layers and convert to a photolysis
c  level grid corresponding to layer centres and boundaries. This is
c  required so that J-values can be calculated for the centre of CTM
c  layers; the index of these layers is kept in the jndlev array.
C--------------------------------------------------------------------
c
c  Set lower boundary and levels to calculate J-values at 
      J1=2*J1-1
      do j=1,LM
        jndlev(j)=2*j
      enddo
c
c  Calculate column optical depths above each level, TTAU
      TTAU(NCFASTJ2+1)=0.0D0
      do J=NCFASTJ2,J1,-1
        I=(J+1)/2
        TTAU(J)=TTAU(J+1) + 0.5d0*DTAUX(I)
        jaddlv(j)=int(0.5d0*DTAUX(I)/dtaumax)
c       
c  Subdivide cloud-top levels if required
        if(jadsub(j).gt.0) then
          jadsub(j)=min(jaddlv(j)+1,nint(dtausub))*(nint(dsubdiv)-1)
          jaddlv(j)=jaddlv(j)+jadsub(j)
        endif
      enddo
c
C---reflect flux from surface
c
      if(U0.gt.0.d0) then
        ZFLUX = U0*FTAU(J1)*RFLECT/(1.d0+RFLECT)
      else
        ZFLUX = 0.d0
      endif
c
c  Calculate attenuated beam, FTAU, level boundaries then level centres
      FTAU(NCFASTJ2+1)=1.0d0
      do J=NCFASTJ2-1,J1,-2
        I=(J+1)/2
        FTAU(J)=FTAU(I)
      enddo
      do J=NCFASTJ2,J1,-2
        FTAU(J)=sqrt(FTAU(J+1)*FTAU(J-1))
      enddo
c
c  Calculate scattering properties, level centres then level boundaries
c  using an inverse interpolation to give correctly-weighted values
      do j=NCFASTJ2,J1,-2
        do i=1,MFIT
          pomegaj(i,j) = pomegaj(i,j/2)
        enddo
      enddo
      do j=J1+2,NCFASTJ2,2
        taudn = ttau(j-1)-ttau(j)
        tauup = ttau(j)-ttau(j+1)
        do i=1,MFIT
          pomegaj(i,j) = (pomegaj(i,j-1)*taudn + 
     $                    pomegaj(i,j+1)*tauup) / (taudn+tauup)
        enddo
      enddo
c  Define lower and upper boundaries
      do i=1,MFIT
        pomegaj(i,J1)   = pomegaj(i,J1+1)
        pomegaj(i,NCFASTJ2+1) = pomegaj(i,NCFASTJ2)
      enddo
c
C--------------------------------------------------------------------
c  Calculate cumulative total and define levels at which we want
c  the J-values.  Sum upwards for levels, and then downwards for Mie
c  code readjustments.
C--------------------------------------------------------------------
c
c  Reinitialize level arrays
      jaddto(:)=0 ! NLFASTJ+1
c
      jaddto(J1)=jaddlv(J1)
      do j=J1+1,NCFASTJ2
        jaddto(j)=jaddto(j-1)+jaddlv(j)
      enddo
      if((jaddto(NCFASTJ2)+NCFASTJ2).gt.NLFASTJ) then
         write(6,1500)  jaddto(NCFASTJ2)+NCFASTJ2, 'NLFASTJ',NLFASTJ
         call stop_model('problem in fastj2 with jaddto',255)
      endif
      jndlev(:)=jndlev(:)+jaddto(jndlev(:)-1) ! LM
      jaddto(NCFASTJ2)=jaddlv(NCFASTJ2)
      do j=NCFASTJ2-1,J1,-1
        jaddto(j)=jaddto(j+1)+jaddlv(j)
      enddo
c
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
C
c  Initialise all Fast-J2 optical property arrays
      pomega(:,:) = 0.d0 ! (2*M__,N__)
      ztau(:)     = 0.d0 ! N__
      fz(:)       = 0.d0 ! N__
c
c  Ascend through atmosphere transposing grid and adding extra points
      do j=J1,NCFASTJ2+1
        k = 2*(NCFASTJ2+1-j)+2*jaddto(j)+1
        ztau(k)= ttau(j)
        fz(k)  = ftau(j)
        pomega(:,k) = pomegaj(:,j) ! MFIT
      enddo
c
c  Check profiles if desired
c      ND = 2*(NCFASTJ2+jaddto(J1)-J1)  + 3
c      if(kw.eq.1) call CH_PROF(ND)
c
C---------------------------------------------------------------------
c  Insert new levels, working downwards from the top of the atmosphere
c  to the surface (down in 'j', up in 'k'). This allows ztau and pomega
c  to be incremented linearly (in a +ve sense), and the flux fz to be
c  attenuated top-down (avoiding problems where lower level fluxes are
c  zero).
C---------------------------------------------------------------------
c
      do j=NCFASTJ2,J1,-1
          zk = 0.5d0/(1.d0+dble(jaddlv(j)-jadsub(j)))
          dttau = (ttau(j)-ttau(j+1))*zk
          dpomega(:) = (pomegaj(:,j)-pomegaj(:,j+1))*zk ! MFIT
c  Filter attenuation factor - set minimum at 1.0d-05
          if(ftau(j+1).eq.0.d0) then
            ftaulog=0.d0
          else
            ftaulog = ftau(j)/ftau(j+1)
            if(ftaulog.lt.1.d-150) then
              ftaulog=1.0d-05
            else
              ftaulog=exp(log(ftaulog)*zk)
            endif
          endif
          k = 2*(NCFASTJ2-j+jaddto(j)-jaddlv(j))+1   !  k at level j+1
          l = 0
c          
c  Additional subdivision of first level if required
          if(jadsub(j).ne.0) then
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
              if(k.gt.1800)then
               write(6,*) 'k fault:',k,NCFASTJ2,j,jaddto(j),jadsub(j),
     *          dsubdiv,jaddlv(j)
               call sys_flush(6)
              endif
            enddo
          endif
          l = 2*(jaddlv(j)-jadsub(j)-l)+1
c
c  Add values at all intermediate levels
          do ix=1,l
            ztau(k+1) = ztau(k) + dttau
            fz(k+1) = fz(k)*ftaulog
            pomega(:,k+1) = pomega(:,k) + dpomega(:) ! MFIT
            k = k+1
          enddo
      enddo
c
C---Update total number of levels and check does not exceed N__
      ND = 2*(NCFASTJ2+jaddto(J1)-J1)  + 3
      if(nd.gt.N__) then
        write(6,1500) ND, 'N__',N__
        call stop_model('problem in fastj2 with ND',255)
      endif
c
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
c
      ZU0 = U0
      ZREFL = RFLECT
c
C-----------------------------------------
      CALL MIESCT(ND)
C-----------------------------------------
C
c  Accumulate attenuation for selected levels
      l=2*(NCFASTJ2+jaddto(J1))+3
      do j=1,LM
        k=l-(2*jndlev(j))
        if(k.gt.ND-2) then
          FMEAN(j) = 0.d0
        else
          FMEAN(j) = FJFASTJ(k)
        endif
      enddo
c
      endif                                     ! WAVEL
c
      return
 1000 format(1x,i3,3(2x,1pe10.4),1x,i3)
 1300 format(1x,50(i3))
 1500 format(' Too many levels in photolysis code: need ',i5,' but ',a,
     $       ' dimensioned as ',i3)
      END SUBROUTINE OPMIE
C
C
       subroutine CH_PROF(ND)
!@sum CH_PROF Check profiles to be passed to MIESCT
!@auth UCI (see note above), GCM incorporation: Drew Shindell,
!@+ modelEifications: Greg Faluvegi
!@ver  1.0 (based on ds4p_fastj2_M23)
!@calls 
C
C**** GLOBAL parameters and variables:
C
      USE TRCHEM_Shindell_COM, only: ZTAU,FZ,POMEGA
c                           
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:

      integer i,j,ND
      write(6,1100) 'lev','ztau','fz  ','pomega( )'
      do i=1,ND
        if(ztau(i).ne.0.d0) then
          write(6,1200) i,ztau(i),fz(i),(pomega(j,i),j=1,8)
        endif
      enddo
      return
 1100 format(1x,a3,4(a9,2x))
 1200 format(1x,i3,11(1x,1pe9.3))
      end subroutine CH_PROF       
C
C
      SUBROUTINE MIESCT(ND)
!@sum MIESCT This is an adaption of the Prather rad transfer code. 
!@+  see comments. 
!@auth UCI (see note above), GCM incorporation: Drew Shindell,
!@+ modelEifications: Greg Faluvegi
!@ver  1.0 (based on ds4p_fastj2_M23)
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
C
C**** GLOBAL parameters and variables:
C
      USE TRCHEM_Shindell_COM, only: NFASTJ,EMU,WTFASTJ,ZFLUX,ZU0,
     & MFASTJ,MFIT,ZREFL,FZ,PM0,PM,FJFASTJ
c                           
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:

      integer I, id, imm, ND
      real*8, parameter ::  cmeq1 = 0.25D0
C
C Fix scattering to 4 Gausss pts = 8-stream.
C Solve eqn of R.T. only for first-order M=1
      ZFLUX = (ZU0*FZ(ND)*ZREFL)/(1.0d0+ZREFL)
      DO I=1,NFASTJ
        CALL LEGND0 (EMU(I),PM0,MFIT)
        DO IMM=MFASTJ,MFIT
          PM(I,IMM) = PM0(IMM)
        ENDDO
      ENDDO
C
      CALL LEGND0 (-ZU0,PM0,MFIT)
      DO IMM=MFASTJ,MFIT
        PM0(IMM) = CMEQ1*PM0(IMM)
      ENDDO
C
      CALL BLKSLV(ND)
C
      DO ID=1,ND,2
        FJFASTJ(ID) = 4.0d0*FJFASTJ(ID) + FZ(ID)
      ENDDO
      RETURN
      END SUBROUTINE MIESCT
C
C
      SUBROUTINE BLKSLV(ND)
!@sum BLKSLV Solves the block tri-diagonal system:
!@+   A(I)*X(I-1) + B(I)*X(I) + C(I)*X(I+1) = H(I)
!@auth UCI (see note above), GCM incorporation: Drew Shindell,
!@+ modelEifications: Greg Faluvegi
!@ver  1.0 (based on ds4p_fastj2_M23)
C
C**** GLOBAL parameters and variables:
C
      USE TRCHEM_Shindell_COM, only: BFASTJ,NFASTJ,RR2,CC,DD,HFASTJ,
     & AFASTJ,DD,RR2,C1,FJFASTJ,WTFASTJ,EMU,AAFASTJ
c                           
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
      integer i, j, k, id, ND
      real*8  sum
C
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
C
C
      SUBROUTINE GEN(ID,ND)
!@sum GEN Generates coefficient matrices for the block tri-diagonal
!@+    system:  A(I)*X(I-1) + B(I)*X(I) + C(I)*X(I+1) = H(I)
!@auth UCI (see note above), GCM incorporation: Drew Shindell,
!@+ modelEifications: Greg Faluvegi
!@ver  1.0 (based on ds4p_fastj2_M23)
C
C**** GLOBAL parameters and variables:
C
      USE TRCHEM_Shindell_COM, only: NFASTJ,MFASTJ,MFIT,POMEGA,PM,PM0,
     & HFASTJ,FZ,SFASTJ,WTFASTJ,WFASTJ,U1,V1,BFASTJ,CC,AFASTJ,EMU,C1,
     & ZTAU,AAFASTJ,ZFLUX,ZREFL
c                           
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
C
      integer id, id0, id1, im, i, j, k, mstart, ND
      real*8  sum0, sum1, sum2, sum3, deltau, d1, d2, surfac
      REAL*8, DIMENSION(8)::  TTMP
C
C---------------------------------------------
      IF(ID.EQ.1 .OR. ID.EQ.ND) THEN
C---------calculate generic 2nd-order terms for boundaries
       ID0 = ID
       ID1 = ID+1
       IF(ID.GE.ND) ID1 = ID-1
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
c       
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
       IF (ID.EQ.1) THEN
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
C
C
      SUBROUTINE LEGND0(X,PL,NFASTJ)
!@sum LEGND0 Calculates ORDINARY LEGENDRE fns of X (real)
!@+   from P[0] = PL(1) = 1,  P[1] = X, .... P[N-1] = PL(N)
!@auth UCI (see note above), GCM incorporation: Drew Shindell,
!@+ modelEifications: Greg Faluvegi
!@ver  1.0 (based on ds4p_fastj2_M23)
!@calls 
c                           
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
c
      INTEGER NFASTJ,I
      REAL*8 X,PL(NFASTJ),DEN
C      
C---Always does PL(2) = P[1]
        PL(1) = 1.D0
        PL(2) = X
        DO I=3,NFASTJ
         DEN = (I-1)
         PL(I) = PL(I-1)*X*(2.d0-1.D0/DEN) - PL(I-2)*(1.d0-1.D0/DEN)
        ENDDO
      RETURN
      END SUBROUTINE LEGND0
c
c
      SUBROUTINE MATIN4(AFASTJ)
!@sum MATIN4 invert 4x4 matrix A(4,4) in place with L-U decomp
!@+   (mjp, old...)
!@auth UCI (see note above), GCM incorporation: Drew Shindell,
!@+ modelEifications: Greg Faluvegi
!@ver  1.0 (based on ds4p_fastj2_M23)
c
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
!@var AFASTJ passed (actually BFASTJ)
      REAL*8 AFASTJ(4,4)
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




