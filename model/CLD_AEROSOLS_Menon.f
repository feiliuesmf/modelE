!     MODULE CLOUD_DROP_PRED
!     USE TRACER_COM
!     USE CLOUDS_COM
!     END MODULE CLOUD_DROP_PRED

      SUBROUTINE GET_CDNC(L,LHX,WCONST,WMUI,AIRM,WMX,DXYPJ,
     *FCLD,CAREA,CLDSAVL,DSS,SMFPML,OLDCDO,OLDCDL,
     *DSU,CDNL1,CDNO1)
      USE CLOUDS_COM
      USE TRACER_COM
      USE CONSTANT,only:mb2kg,LHE,LHS
      IMPLICIT NONE
      real*8 CAREA,CLDSAVL,AIRM,WMX,
     * SMFPML,OLDCDO,OLDCDL

      real*8,dimension(9)::DSS,DSU
      real*8 EXPL,EXPO,WCDNO,WCDNL,CDNO0,CDNL0,
     *CCLD0,CCLD1,DCLD,dfn,CDNL1,CDNO1,amass,tams,smalphaf
     *,FCLD,WCONST,LHX,WMUI,DXYPJ
      real*8 SSM1,SSM2,SSM3,SSM4,SSM5,SSMAL,SSMAO,SSML,SSMO
      
      integer L,n
      do n = 1,9
        DSU(n)=0                         
      end do

!stop "CLD_AEROSOLS_MENON"


!add in terms for AMASS from other program to get aerosol mass conc.
      amass=AIRM*mb2kg*DXYPJ  
      tams=1.d0/amass*1.292d0
      do n = 1,9
         DSU(n) =1.d9*DSS(n)*tams
c        if(DSU(n).gt.90.d0)write(6,*)"MASS",DSU(n),n,l,tams,DSS(n)
      enddo

C** use CTEI effect for CDNC as in Menon et al. 2002 JAS
      smalphaf=(1.d0+ 2.d0*SMFPML)*0.5d0 
C** Converting aerosol mass to number
      SSM1=DSU(1)/1769.d0                       !sulfate

      SSM2=(DSU(2)/2169.)/(0.004189*(.44**3.))  !seasalt 0.1-1 um range
      SSM3=(DSU(3)/2169.)/(0.004189*(1.7**3.))  !seasalt 1-4 um range    

      SSM4=((1.0d0*DSU(4))+(1.0d0*DSU(6)))/1000.d0   !OCI and BCI aged use 100% & 100% 
      SSM5=((0.8d0*DSU(5))+(0.6d0*DSU(7)))/1000.d0   !OCB and BCB      use 80%  & 60%
   
      SSMAL=SSM1+SSM4+SSM5
      SSMAO=SSMAL+SSM2!+SSM3   !add in larger size when including TKE effects

      SSML=SSMAL/(0.004189d0*(.050d0**3.d0))
      SSMO=SSMAO/(0.004189d0*(.085d0**3.d0))

      IF(SSML.le.100.d0) SSML=100.d0
      IF(SSMO.le.50.d0) SSMO=50.d0
C** Here we use Gultepe's paramet for CDNC = f(Na)
      EXPL=(298.d0*log10(SSML))-595.d0
      EXPO=(162.d0*log10(SSMO))-273.d0

      IF (EXPO.LT.10.d0) EXPO=10.d0
      IF (EXPL.LT.10.d0) EXPL=10.d0
      WCDNO= EXPO*smalphaf
      WCDNL= EXPL*smalphaf
      CDNL0=OLDCDL   !term initialised to 10. CLOUDS_DRV                  
      CDNO0=OLDCDO   

C** Using the new CDNC scheme where we calculate it as a function of
C** gas phase sulfate and cloud area changes
      CCLD0 = CLDSAVL   !cld frac from previous time step   was 3D
      CCLD1 = 1.d0-CAREA  !CLDSS(L)cld frac from present time step 
      DCLD = CCLD1-CCLD0     ! cloud fraction change
      dfn = 0.0d0
      IF(LHX.EQ.LHE.AND.WMX/FCLD.GE.WCONST*1.d-3) dfn=0.1d0  
      IF(LHX.EQ.LHS.AND.WMX/FCLD.GE.WMUI*1.d-3) dfn=0.1d0

C** If previous time step is clear sky
      if (CCLD0.eq.0.d0) then
        CDNL1 = WCDNL        !EXPL1 = N from present time step
        CDNO1 = WCDNO
C** If previous time step is cloudy then depending on cld frac change
      elseif (DCLD.le.0.d0) then
        CDNL1= CDNL0 - dfn*CDNL0  !EXPL =  N from previous time step
        CDNO1= CDNO0 - dfn*CDNO0  !dfn=fraction of cloud area raining
      elseif (DCLD.gt.0.d0) then
        CDNL1 = (((CDNL0*CCLD0)+(WCDNL*DCLD))/CCLD1) - dfn*CDNL0
        CDNO1 = (((CDNO0*CCLD0)+(WCDNO*DCLD))/CCLD1) - dfn*CDNO0
      endif

      IF (CDNL1.le.10.d0) CDNL1=10.d0
      IF (CDNO1.le.10.d0) CDNO1=10.d0
!     if(CDNL1.gt.1400.d0.or.CDNO1.gt.1400.d0) 
!    *write(6,*) "Nc1st",CDNL1,CDNO1,DSU(1),DSU(2),DSU(3),DSU(4),DSU(5),
!    *DSU(6),DSU(7),L,SSML,SSMO

      RETURN
      
      END SUBROUTINE GET_CDNC 


      SUBROUTINE GET_QAUT(L,TL,FCLD,WMX,SCDNCW,RHO,QCRIT,QAUT)
      USE TRACER_COM
      USE CONSTANT,only:TWOPI,GRAV
      IMPLICIT NONE
      real*8 TL,WMX,dynvis                         
      real*8 QCRIT,rcr,SCDNCW,QAUT,RHO,FCLD
      integer L
   
C**sm Tripoli and Cotton (1980), adapted from Jones Tech rpt.
C     using Qcrit for QAUT
      dynvis=2.4784d0*1.d-6*TL/(TL+393.16d0)
      rcr=7.0d0
      QCRIT=(2.d0/3.d0)*TWOPI*1.d9*((rcr*1.d-06)**3.d0)*SCDNCW/
     *(.001d0*RHO)
   
      QAUT=(0.104d0*GRAV*0.55d0*((.001d0*RHO)**(4.d0/3.d0))*
     * ((WMX/(FCLD+1.d-20))**(7.d0/3.d0)))/
     * (dynvis*((SCDNCW*1.d09)**(1.d0/3.d0)))

      if(QAUT.le.0.) write(6,*)"QCR",QAUT,SCDNCW,WMX,l

      RETURN
    
      END SUBROUTINE GET_QAUT

      SUBROUTINE GET_CDNC_UPD(L,LHX,WCONST,WMUI,WMX,FCLD,
     *CLDSSL,CLDSAVL,DSU,SMFPML,OLDCDO,OLDCDL,CDNL1,CDNO1)
      USE CLOUDS_COM
      USE CONSTANT,only:LHE,LHS
      USE TRACER_COM
      IMPLICIT NONE
      real*8 ::CLDSSL,CLDSAVL,WMX
     *,SMFPML,OLDCDO,OLDCDL
      real*8 EXPL,EXPO,WCDNO,WCDNL,CDNO0,CDNL0,
     *CCLD0,CCLD1,DCLD,dfn,CDNL1,CDNO1,smalfaf,FCLD
     *,LHX,WMUI,WCONST
      real*8,dimension(9)::DSU

      real*8 SSM1,SSM2,SSM3,SSM4,SSM5,SSMAL,SSMAO,SSML,SSMO
       integer L

      smalfaf=(1.d0+ 2.d0*SMFPML)*0.5d0   !SMFPM was 3D
C** Converting aerosol mass to number
      SSM1=DSU(1)/1769.d0                       !sulfate

      SSM2=(DSU(2)/2169.)/(0.004189*(.44**3.))  !seasalt 0.1-1 um range
      SSM3=(DSU(3)/2169.)/(0.004189*(1.7**3.))  !seasalt 1-4 um range    

      SSM4=((1.0d0*DSU(4))+(1.0d0*DSU(6)))/1000.d0   !OCI and BCI aged use 100% & 100% 
      SSM5=((0.8d0*DSU(5))+(0.6d0*DSU(7)))/1000.d0   !OCB and BCB      use 80%  & 60%
   
      SSMAL=SSM1+SSM4+SSM5
      SSMAO=SSMAL+SSM2!+SSM3   !add in larger size when including TKE effects

      SSML=SSMAL/(0.004189d0*(.050d0**3.d0))
      SSMO=SSMAO/(0.004189d0*(.085d0**3.d0))

      IF(SSML.le.100.d0) SSML=100.d0
      IF(SSMO.le.50.d0) SSMO=50.d0
C** Here we use Gultepe's paramet for CDNC = f(Na)
      EXPL=(298.d0*log10(SSML))-595.d0
      EXPO=(162.d0*log10(SSMO))-273.d0

      IF (EXPO.LT.10.d0) EXPO=10.d0
      IF (EXPL.LT.10.d0) EXPL=10.d0
      WCDNO= EXPO*smalfaf
      WCDNL= EXPL*smalfaf
      CDNL0=OLDCDL
      CDNO0=OLDCDO 

C** Using the new CDNC scheme where we calculate it as a function of
C** gas phase sulfate and cloud area changes
      CCLD0 = CLDSAVL   ! cld frac from previous time step  was 3D
      CCLD1 = CLDSSL    !cld frac from present time step  was 3D
      DCLD = CCLD1-CCLD0   ! cloud fraction change
      dfn = 0.d0
      IF(LHX.EQ.LHE.AND.WMX/FCLD.GE.WCONST*1.d-3) dfn=0.1d0
      IF(LHX.EQ.LHS.AND.WMX/FCLD.GE.WMUI*1.d-3) dfn=0.1d0

C** If previous time step is clear sky
      if (CCLD0.eq.0.d0) then
        CDNL1 = WCDNL        !EXPL1 = N from present time step
        CDNO1 = WCDNO
C** If previous time step is cloudy then depending on cld frac change
      elseif (DCLD.le.0.d0) then
        CDNL1= CDNL0 - dfn*CDNL0  !EXPL =  N from previous time step
        CDNO1= CDNO0 - dfn*CDNO0  !dfn=fraction of cloud area raining
      elseif (DCLD.gt.0.d0) then
        CDNL1 = (((CDNL0*CCLD0)+(WCDNL*DCLD))/CCLD1) - dfn*CDNL0
        CDNO1 = (((CDNO0*CCLD0)+(WCDNO*DCLD))/CCLD1) - dfn*CDNO0
      endif
c     if(CDNL1.gt.1400.d0.or.CDNO1.gt.1400.d0) 
c    *write(6,*) "Nc2nd",CDNL1,CDNO1,DSU(1),DSU(2),DSU(3),DSU(4),DSU(5),
c    *DSU(6),DSU(7),L,SSML,SSMO
      IF (CDNL1.le.10.d0) CDNL1=10.d0
      IF (CDNO1.le.10.d0) CDNO1=10.d0

      RETURN

      END SUBROUTINE GET_CDNC_UPD
