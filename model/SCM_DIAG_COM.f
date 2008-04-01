!     SCM_DIAG_COM.f    
!     save diagnostics for run of MODELE SCM

      MODULE SCMDIAG  

      USE RESOLUTION , only : LM 

      IMPLICIT NONE


      real*8 LHPSAV(LM),LHPMC(LM),PRESAV(LM),PREMC(LM)

      real*8   CLCVSS(LM),CLCVMC(LM),CLTHCK(LM),CSIZE(LM,2),   
     *         EFFRAD(LM),TAUSSC(LM),TAUMCC(LM),CUMFLX(LM),
     *         CUMHET(LM),CUMOST(LM),PRCSS,PRCMC,EVPFLX,SHFLX,
     *         SOILMS,CLDFLG(LM),DWNFLX(LM),RHC(LM)
      real*8   clsav(LM)
      real*8 SRDFLBTOP,SRNFLBTOP,TRUFLBTOP,SRDFLBBOT,SRNFLBBOT,
     *       TRUFLBBOT,TRDFLBBOT
      real*8, DIMENSION(LM) :: SRFHRLCOL,TRFCRLCOL
c
c     cumulus updraft speed
      real*8   WCUSCM(LM,2),WCUALL(LM,2,LM)
      real*8   WCUDEEP(LM,2)
C--- Added by J.W. starting ---C
      real*8   MPLUMESCM(LM,2),MPLUMEALL(LM,2,LM)
      real*8   MPLUMEDEEP(LM,2)
      real*8   ENTSCM(LM,2),ENTALL(LM,2,LM)
      real*8   ENTDEEP(LM,2)
      real*8   DETRAINDEEP(LM,2,LM)
C--- Added by J.W. ending ---C
c
c     precipitating and non-precipitating convective condensate for
c     deep convection 
      real*8 PRCCDEEP(LM,2,LM),NPRCCDEEP(LM,2,LM),TPALL(LM,2,LM)   
      real*8 PRCCGRP(LM,2,LM),PRCCICE(LM,2,LM)
     
c     condensate for all convection
      real*8 mccond(LM,2,LM)


      END MODULE SCMDIAG  
