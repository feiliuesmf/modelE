!     MODULE CLOUD_DROP_PRED
!     USE TRACER_COM
!     USE CLOUDS_COM
!     END MODULE CLOUD_DROP_PRED

      SUBROUTINE GET_CDNC(L,LHX,WCONST,WMUI,AIRM,WMX,DXYPJ,
     *FCLD,CAREA,CLDSAVL,DSGL,SMFPML,OLDCDO,OLDCDL,DSU,CDNL1,CDNO1)
      USE CLOUDS_COM
      USE TRACER_COM
      USE CONSTANT,only:mb2kg,LHE,LHS
      IMPLICIT NONE
      real*8 CAREA,CLDSAVL,AIRM,WMX,DSGL,
     * SMFPML,OLDCDO,OLDCDL

      real*8 EXPL,EXPO,SSM7,SSM8,SSM1,WCDNO,WCDNL,CDNO0,CDNL0,
     *CCLD0,CCLD1,DCLD,dfn,CDNL1,CDNO1,amass,tams,smalphaf
     *,FCLD,DSU,WCONST,LHX,WMUI,DXYPJ
      integer L

!add in terms for AMASS from other program to get aerosol mass conc.
      amass=AIRM*mb2kg*DXYPJ  
      tams=1.d0/amass*1.292d0
      DSU =(1.d9*DSGL)*tams          
C** use CTEI effect for CDNC as in Menon et al. 2002 JAS
      smalphaf=(1.d0+ 2.d0*SMFPML)*0.5d0 
C** Converting aerosol mass to number
      SSM1=DSU/1769.d0
      SSM7=SSM1/(0.004189d0*(.050d0**3.d0))
      SSM8=SSM1/(0.004189d0*(.085d0**3.d0))
      IF(SSM7.le.100.d0) SSM7=100.d0
      IF(SSM8.le.50.d0) SSM8=50.d0
C** Here we use Gultepe's paramet for CDNC = f(Na)
      EXPL=(298.d0*log10(SSM7))-595.d0
      EXPO=(162.d0*log10(SSM8))-273.d0

      IF (EXPO.LT.10.d0) EXPO=10.0d0
      IF (EXPL.LT.10.d0) EXPL=10.0d0
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
!     if(CDNL1.gt.1000.d0) 
!    *write(6,*) "CDNC1st",CDNL1,DSU,smalphaf,L

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

c     if(QAUT.eq.0.) write(6,*)"QCR",QAUT,SCDNCW,WMX(L),l

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
      real*8 EXPL,EXPO,SSM7,SSM8,SSM1,WCDNO,WCDNL,CDNO0,CDNL0,
     *CCLD0,CCLD1,DCLD,dfn,CDNL1,CDNO1,smalfaf,FCLD,DSU
     *,LHX,WMUI,WCONST
       integer L

      smalfaf=(1.d0+ 2.d0*SMFPML)*0.5d0   !SMFPM was 3D
C** Converting aerosol mass to number
      SSM1=DSU/1769.d0
      SSM7=SSM1/(0.004189d0*(.050d0**3.d0))
      SSM8=SSM1/(0.004189d0*(.085d0**3.d0))
      IF(SSM7.le.100.d0) SSM7=100.d0
      IF(SSM8.le.50.d0) SSM8=50.d0
C** Here we use Gultepe's paramet for CDNC = f(Na)
      EXPL=(298.d0*log10(SSM7))-595.d0
      EXPO=(162.d0*log10(SSM8))-273.d0

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
!     if (CDNL1.gt.1000.) 
!    *write(6,*) "CDNCUPD",CDNL1,DSU,smalfaf
      IF (CDNL1.le.10.d0) CDNL1=10.d0
      IF (CDNO1.le.10.d0) CDNO1=10.d0

      RETURN

      END SUBROUTINE GET_CDNC_UPD
