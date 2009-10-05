!     SCMDATA_TWPICE_COM.f
!@sum set up variables for SCM forcing data 
!@auth  Audrey Wolf
!@ver  Version to read TWP ICE  IOP                          
!
C--------------------------------------------------------------------------------
      Module SCMCOM
      USE RESOLUTION , ONLY : LM
      IMPLICIT NONE
      SAVE

C SCM DATA as provided from ARM variational analysis 
!@var SG_P Pressure at GCM sigma levels (mb)
      REAL*8 SG_P(LM)
      REAL*8 SGE_P(LM)
!@var SG_HGT Height of layer at GCM sigma levels (m)
      REAL*8 SG_HGT(LM)
!@var SG_T ARM Temperature at GCM sigma levels (K)
      REAL*8 SG_T(LM)
!@var SG_Q ARM Specific Humidity at GCM sigma levels (g/kg)
      REAL*8 SG_Q(LM)
!@var SG_U ARM U wind at GCM Sigma levels (m/s)
      REAL*8 SG_U(LM)
!@var SG_V ARM V wind at GCM sigma levels (m/s)
      REAL*8 SG_V(LM)
!@var SG_OMEGA ARM Omega at GCM sigma levels (mb/hr)
      REAL*8 SG_OMEGA(LM)
!@var SG_WINDIV ARM Wind Divergence at GCM sigma levels (1/s)
      REAL*8 SG_WINDIV(LM)
!@var SG_CONV   convergence as calculated in FCONV subr from windiv for clouds   
      REAL*8 SG_CONV(LM)   
!@var SG_HOR_TMP_ADV ARM Horizontal temperature advection at GCM sigma levels (K/s)
      REAL*8 SG_HOR_TMP_ADV(LM)
!@var SG_VER_S_ADV ARM Vertical S advection at GCM sigma levels (K/s)
      REAL*8 SG_VER_S_ADV(LM)
!@var SG_HOR_Q_ADV ARM Horizontal Q advection at GCM sigma levels( kg/kg/s)
      REAL*8 SG_HOR_Q_ADV(LM)
!@var SG_VER_Q_ADV ARM Vertical Q advection at GCM sigma levels (kg/kg/s) 
      REAL*8 SG_VER_Q_ADV(LM)
!@var SG_ARSCL ARSCL Cloud amounts at GCM sigma levels (%)
      REAL*8 SG_ARSCL(LM)
!@var ASTIME ARM TIME STAMP - surface data 
      REAL*8 ASTIME
!@var ALTIME ARM TIME STAMP - layer data
      REAL*8 ALTIME
!@var APREC ARM precipitation (mm/hour)
      REAL*8 APREC
!@var ALH Latent Heat Flux (W/m**2)
      REAL*8 ALH
!@var ASH Sensible Heat Flux (W/m**2)
      REAL*8 ASH
!@var AMEANPS ARM Mean Surface Pressure (mb) 
      REAL*8 AMEANPS
!@var ATSAIR ARM surface air temperature (C)
      REAL*8 ATSAIR
!@var ATSKIN surface skin temperature (C)
      REAL*8 ATSKIN
!@var ARHSAIR ARM Surface Air Relative Humidity (%)
      REAL*8 ARHSAIR
!@var ASWINDSPD 10m  wind speed (m/s)
      REAL*8 ASWINDSPD
!@var AQS 2m water vapor mixing ration(kg/kg)
      REAL*8 AQS
!@var AUS 10m u component (m/s)
      REAL*8 AUS
!@var AVS 10m v component (m/s)
      REAL*8 AVS 
!@var AUSRF surface u wind (m/s)
      REAL*8 AUSRF
!@var AVSRF surface v wind (m/s) 
      REAL*8 AVSRF
!@var ALWP  ARM MWR Cloud Liquid Water Path (cm)
      REAL*8 ALWP
!@var ADWDT ARM d(Column H20)/dt  (mm/hr)
      REAL*8 ADWDT
!@var ADWADV ARM Column_H20_Advection_ (mm/hr)
      REAL*8 ADWADV
!@var ATLWUP ARM TOA LW UP (W/m**2)
      REAL*8 ATLWUP
!@var ATSWDN ARM TOA SW DOWN (W/m**2)
      REAL*8 ATSWDN
!@var ATSWIN ARM TOA SW INS (W/m**2)
      REAL*8 ATSWIN
  
!@var SCM_SURFACE_FLAG 0-use GCM calculated surface fluxes,
!                      1-use SCM prescribed surface fluxes  
!                      2-use SCM surf temps and GCM calculated surface fluxes
      INTEGER SCM_SURFACE_FLAG
!@var SCM_ATURB_FLAG   0-run with dry convection routine
!                      1-run with aturb routine
      INTEGER SCM_ATURB_FLAG
!@var NARM #of GCM time steps per ARM time step
      INTEGER NARM
!@var NRINIT #of GCM time steps between reinitializing T,Q 
      INTEGER NRINIT
!@var TAUARM starting TAU of ARM DATA  (should this be starting date and time)
      INTEGER TAUARM
!@var IKT index to arm data interpolated to time steps
      INTEGER IKT
      INTEGER iu_scm_prt,iu_scm_diag    
   
      
      INTEGER MCT

      parameter (MCT=1500)

      REAL*8 HTA_HR(LM,MCT)    
      REAL*8 VSA_HR(LM,MCT)
      REAL*8 HQA_HR(LM,MCT)
      REAL*8 VQA_HR(LM,MCT)
    
      REAL*4 STMSTEP(MCT)
      REAL*4 STMSTEPL(MCT)
      REAL*4 AMPS(MCT)
      REAL*4 ATSK(MCT)
      REAL*4 ATSA(MCT)
      REAL*4 ARHHR(MCT)
      REAL*4 THR(LM,MCT)
      REAL*4 QHR(LM,MCT)
      REAL*4 APRCHR(MCT)
      REAL*4 ALHHR(MCT)
      REAL*4 ASHHR(MCT)
      REAL*4 AWSHR(MCT)
      REAL*4 AQSHR(MCT)
      REAL*4 AUSHR(MCT)
      REAL*4 AVSHR(MCT)
      REAL*4 UHR(LM,MCT)
      REAL*4 VHR(LM,MCT)
      REAL*4 OMGHR(LM,MCT)
      REAL*4 WDHR(LM,MCT)
      REAL*4 ACLDHR(LM,MCT)
      REAL*4 ALWPHR(MCT)
      REAL*4 ADWDTHR(MCT)
      REAL*4 ADWADVHR(MCT)
      REAL*4 ATLWUPHR(MCT) 
      REAL*4 ATSWDNHR(MCT)
      REAL*4 ATSWINHR(MCT)

      REAL*8 SCM_SAVE_T(LM),SCM_SAVE_Q(LM),SCM_DEL_T(LM),
     &       SCM_DEL_Q(LM)


      end module SCMCOM
c
c    
      subroutine ALLOC_SCM_COM()
   
      USE SCMCOM, only : SCM_SURFACE_FLAG,SCM_ATURB_FLAG


!@var SCM_SURFACE FLAG   0-run with GCM calculated surface fluxes
!                       1-run with ARM prescribed surface fluxes
!                       2-RUN WITH ARM srf tmps and GCM calc srf fluxes
      SCM_SURFACE_FLAG = 1     

!@var SCM_ATURB_FLAG     0-run with DRYCNV dry convection routine 
!                        1-run with ATURB turbulence routine    
      SCM_ATURB_FLAG = 1
      if (SCM_ATURB_FLAG.eq.0) then
          write(0,*) 'RUN with DRYCNV routine '
      elseif (SCM_ATURB_FLAG.eq.1) then
          write(0,*) 'RUN with ATURB routine '
      endif    
  
      return

      end subroutine ALLOC_SCM_COM
