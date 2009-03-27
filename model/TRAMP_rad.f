#include "rundeck_opts.h"

c -----------------------------------------------------------------


      SUBROUTINE SETAMP(EXT,SCT,GCB,TAB)
!@sum Calculation of extinction, asymmetry and scattering for AMP Aerosols
!@sum Calculation of absorption in the longwave
!@sum Called in SETAER / RCOMPX
!@auth Susanne Bauer

      USE AMP_AEROSOL, only: AMP_EXT, AMP_ASY, AMP_SCA,
     &                       Reff_LEV, NUMB_LEV, RindexAMP, AMP_Q55, dry_Vf_LEV
      USE AERO_CONFIG, only: NMODES

      USE MODEL_COM,   only: lm,itime,itimeI
      USE TRACER_COM,  only: TRM
      USE RADPAR,      only: TTAUSV,aesqex,aesqsc,aesqcb,FSTOPX,FTTOPX !Diagnostics

      IMPLICIT NONE

      ! Arguments: Optical Parameters dimension(lm,wavelength)
      REAL(8), INTENT(OUT) :: EXT(LM,6)       ! Extinction, SW
      REAL(8), INTENT(OUT) :: SCT(LM,6)       ! Single Scattering Albedo, SW
      REAL(8), INTENT(OUT) :: GCB(LM,6)       ! Asymmetry Factor, SW
      REAL(8), INTENT(OUT) :: TAB(LM,33)      ! Thermal absorption Cross section, LW


      ! Local
      
      INTEGER l,n,w,MA,MB,MD,NA,NS
      REAL*8 Size(23), Mie_IM(17), Mie_RE(15), HELP, AMP_TAB(33), CORE_CLASS(nmodes), SHELL_CLASS(nmodes),Vf(6)
      DATA Size/0.002, 0.005,0.01,0.05,0.08,0.1,0.13,0.17,0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.8,1.0,1.2,1.5,2.,3.,5.,10./
      DATA Mie_RE/1.25,1.3,1.35,1.4,1.45,1.5,1.55,1.6,1.65,1.7,1.75,1.8,1.85,1.9,1.95/
      DATA Mie_IM/0.0,0.00001,0.00002,0.00005,0.0001,0.0002,0.0005,0.001,0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.5,1.0/
c                        AKK  ACC  DD1  DS1  DD2  DS2  SSA  SSC  OCC  BC1  BC2  BC3  DBC  BOC  BCS  MXX
c                        1    2    3    4    5    6    7    8    9    10   11   12   13   14   15   16
      DATA CORE_CLASS   /1,   1,   6,   6,   6,   6,   2,   2,   4,   5,   5,   5,   6,   4,   5,   6/
      DATA SHELL_CLASS  /0,   0,   0,   0,   0,   0,   0,   0,   1,   0,   0,   0,   0,   1,   1,   2/
c                   NA1= SO4  NA2=SS  NA3=NO3 NA4=OC NA5=BC NA6=DU
      EXT(:,:)    = 0.d0
      SCT(:,:)    = 0.d0
      GCB(:,:)    = 0.d0
      TAB(:,:)    = 0.d0
      TTAUSV(:,:) = 0.d0

      if (itime.ne.itimeI) then 
c Shortwave: ---------------------------------------------------------------------------------------------    

      DO l = 1,lm
      DO n = 1,nmodes
         w = 6    ! aot at 550
         do MA = 1,17
            if ( aimag   (RindexAMP(l,n,w)) .le. Mie_IM(MA)) goto 100
         enddo
 100      continue
         do MB = 1,15
            if ( real    (RindexAMP(l,n,w)) .le. Mie_RE(MB)) goto 200
         enddo
 200      continue
         do MD = 1,23
            if (Reff_LEV(l,n) .le. Size(md)) goto 300             
         enddo
 300      continue  
          MA = min(17,MA)
          MB = min(15,MB)
          MD = min(23,MD)
          TTAUSV(l,n) = NUMB_LEV(l,n) * AMP_Q55(MA,MB,MD) 

      DO w = 1,6  !wavelength
         do MA = 1,17
            if ( aimag   (RindexAMP(l,n,w)) .le. Mie_IM(MA)) goto 101
         enddo
 101     continue
         do MB = 1,15
            if ( real    (RindexAMP(l,n,w)) .le. Mie_RE(MB)) goto 201
         enddo
 201     continue
         do MD = 1,23
            if (Reff_LEV(l,n) .le. Size(md)) goto 301
         enddo
 301     continue               !  [#/m2]         [um2] ??
          MA = min(17,MA)
          MB = min(15,MB)
          MD = min(23,MD)

          EXT(l,w) = EXT(l,w) + ( AMP_EXT(MA,MB,w,MD) * TTAUSV(l,n)) * FSTOPX(n)
          HELP     = ((GCB(l,w) *  SCT(l,w) ) + (AMP_ASY(MA,MB,w,MD) * AMP_SCA(MA,MB,w,MD) * TTAUSV(l,n)) ) 
          SCT(l,w) = SCT(l,w) +  (AMP_SCA(MA,MB,w,MD) * TTAUSV(l,n)) * FSTOPX(n)
          GCB(l,w) = HELP / (SCT(l,w)+ 1.D-10)
          GCB(l,w) = GCB(l,w) * FSTOPX(n)

          aesqex(l,w,n)= AMP_EXT(MA,MB,w,MD) * TTAUSV(l,n) 
          aesqsc(l,w,n)= AMP_SCA(MA,MB,w,MD) * TTAUSV(l,n) 
          aesqcb(l,w,n)= AMP_ASY(MA,MB,w,MD) * aesqsc(l,w,n)

      ENDDO   ! wave
C Longwave: ---------------------------------------------------------------------------------------------
      NA = CORE_CLASS(n)
      NS = SHELL_CLASS(n)
      Vf(:)=dry_Vf_LEV(l,n,1:6)
      CALL GET_LW(NA,NS,Reff_LEV(l,n),AMP_TAB,Vf)
         TAB(l,:) = TAB(l,:) + AMP_TAB(:) *  TTAUSV(l,n) * FTTOPX(n)
      ENDDO   ! modes
      ENDDO   ! level

c    write ss diagnostic on ds1 and ds2
        TTAUSV(:,4) =  TTAUSV(:,7) 
        TTAUSV(:,6) =  TTAUSV(:,8) 
        aesqex(:,:,4)= aesqex(:,:,7)
        aesqex(:,:,6)= aesqex(:,:,8)
        aesqsc(:,:,4)= aesqsc(:,:,7)
        aesqsc(:,:,6)= aesqsc(:,:,8)
        aesqcb(:,:,4)= aesqcb(:,:,7)
        aesqcb(:,:,6)= aesqcb(:,:,8)


      endif
  
      RETURN
      END SUBROUTINE SETAMP
c -----------------------------------------------------------------

c -----------------------------------------------------------------
      SUBROUTINE SETAMP_LEV(i,j,l)
!@sum Calulates effective Radius and Refractive Index for Mixed Aerosols
!@sum Puts AMP Aerosols in 1 dimension CALLED in RADIA
!@auth Susanne Bauer

      USE AMP_AEROSOL, only: DIAM, Reff_LEV, NUMB_LEV, RindexAMP,NUMB_SS,dry_Vf_LEV
      USE TRACER_COM,  only: TRM, ntmAMP, AMP_AERO_MAP,AMP_NUMB_MAP,AMP_MODES_MAP,trname
      USE AERO_CONFIG, only: NMODES
      USE AERO_SETUP,  only: SIG0, CONV_DPAM_TO_DGN   !(nmodes * npoints) lognormal parameters for each mode
      USE GEOM,        only: BYAXYP ! inverse area of gridbox [m-2]
      IMPLICIT NONE

      ! Arguments: 
      INTEGER, INTENT(IN) :: i,j,l

      ! Local
      INTEGER n,w,s
      REAL*8,     DIMENSION(nmodes,7) :: VolFrac, Mass
      REAL*8                          :: H2O, NO3 
      REAL(8), PARAMETER :: TINYNUMER = 1.0D-30 
      COMPLEX*8, DIMENSION(6,7)      :: Ri
c Andies data incl Solar weighting - integral over 6 radiation band
      DATA Ri/(1.46099,    0.0764233)  ,(1.48313,  0.000516502),    !Su
     &        (1.49719,  1.98240e-05)  ,(1.50793,  1.64469e-06),
     &        (1.52000,  1.00000e-07)  ,(1.52815,  1.00000e-07),

     &        (1.80056,     0.605467)  ,(1.68622,     0.583112),    !Bc
     &        (1.63586,     0.551897)  ,(1.59646,     0.515333), 
     &        (1.57466,     0.484662)  ,(1.56485,     0.487992),

     &        (1.46099,    0.0761930)  ,(1.48313,   0.00470000),    !Oc
     &        (1.49719,   0.00470000)  ,(1.50805,   0.00480693),
     &        (1.52000,   0.00540000)  ,(1.52775,    0.0144927), 

     &        (1.47978,    0.0211233)  ,(1.50719,   0.00584169),    !Du
     &        (1.51608,   0.00378434)  ,(1.52998,   0.00178703),
     &        (1.54000,  0.000800000)  ,(1.56448,   0.00221463),

     &        (1.46390,   0.00571719)  ,(1.45000,      0.00000),    !Ss
     &        (1.45000,      0.00000)  ,(1.45000,      0.00000),
     &        (1.45000,      0.00000)  ,(1.45000,      0.00000),
 
     &        (1.46099,  0.0764233)    ,(1.48313,  0.000516502),    !No3
     &        (1.49719,  1.98240e-05)  ,(1.50793,  1.64469e-06),
     &        (1.52000,  1.00000e-07)  ,(1.52815,  1.00000e-07),

     &        (1.26304,    0.0774872)  ,(1.31148,  0.000347758),    !H2O
     &        (1.32283,  0.000115835)  ,(1.32774,  3.67435e-06),
     &        (1.33059,  1.58222e-07)  ,(1.33447,  3.91074e-08)/

    ! + Effective Radius [um] per Mode
      DO n=1,nmodes
      Reff_LEV(l,n) = DIAM(i,j,l,n) * CONV_DPAM_TO_DGN(n)* 0.5e6
      ENDDO

      ! + Mass and Number Concentration
       DO n=1,ntmAMP 
           if(trname(n) .eq.'M_NO3') NO3 =trm(i,j,l,n)
           if(trname(n) .eq.'M_H2O') H2O =trm(i,j,l,n)
          if(AMP_NUMB_MAP(n).eq. 0) then  ! Volume fraction
           if(trname(n)(6:8).eq.'_SU') Mass(AMP_MODES_MAP(n),1) =trm(i,j,l,n)
           if(trname(n)(6:8).eq.'_BC') Mass(AMP_MODES_MAP(n),2) =trm(i,j,l,n)
           if(trname(n)(6:8).eq.'_OC') Mass(AMP_MODES_MAP(n),3) =trm(i,j,l,n)
           if(trname(n)(6:8).eq.'_DU') Mass(AMP_MODES_MAP(n),4) =trm(i,j,l,n)
           if(trname(n)(6:8).eq.'_SS') Mass(AMP_MODES_MAP(n),5) =trm(i,j,l,n)
          else                           ! Number
!          [ - ]                        [#/gb]         [m-2]      
           NUMB_LEV(l,AMP_NUMB_MAP(n)) =trm(i,j,l,n) * byaxyp(i,j)
        endif
         ENDDO

      NUMB_LEV(l,7) = NUMB_SS(i,j,l,1) *  byaxyp(i,j)
      NUMB_LEV(l,8) = NUMB_SS(i,j,l,2) *  byaxyp(i,j)

      ! + Volume Fraction
      DO n=1,nmodes  ![#/m2]         pi/4     [m2]
      NUMB_LEV(l,n) = NUMB_LEV(l,n)* 0.7853 * DIAM(i,j,l,n)**2
       ! NO3   
      Mass(n,6) = Mass(n,1) / Sum(Mass(:,1)) * NO3
      ! H2O
      Mass(n,7) = Mass(n,1) /(Sum(Mass(:,1)) + TINYNUMER)  * H2O
      ENDDO
      DO s=1,7  ! loop over species 
      DO n=1,nmodes  ! loop over modes
      Volfrac(n,s) = Mass(n,s) / (Sum(Mass(n,:)) + TINYNUMER)
      dry_Vf_LEV(l,n,s) = Mass(n,s) / (Sum(Mass(n,1:6)) + TINYNUMER)
      ENDDO
      ENDDO
 
      ! + Refractive Index of Aerosol mix per mode and wavelength
      
      RindexAMP(l,:,:) = 0.d0
      DO s=1,7  ! loop over species 
      DO w=1,6  ! loop over wavelength
      DO n=1,nmodes  ! loop over modes
      RindexAMP(l,n,w) = RindexAMP(l,n,w) + ( Volfrac(n,s) * Ri(w,s))
      ENDDO
      ENDDO
      ENDDO
  
      RETURN
      END SUBROUTINE SETAMP_LEV
c -----------------------------------------------------------------

c -----------------------------------------------------------------
      SUBROUTINE SETUP_RAD
!@sum Initialization for Radiation incl. Aerosol Microphysics
!@auth Susanne Bauer

      USE AMP_AEROSOL, only: AMP_EXT, AMP_ASY, AMP_SCA, AMP_Q55
	  
	  IMPLICIT NONE
      include 'netcdf.inc'
      integer start(4),count(4),count3(3),status
      integer ncid, id1, id2, id3, id4

      real*4, DIMENSION(17,15,6,23) ::  ASY,SCA,EXT
      real*4, DIMENSION(17,15,23) ::    QEX
c -----------------------------------------------------------------
c   Opening of the files to be read
c -----------------------------------------------------------------
          status=NF_OPEN('AMP_MIE_TABLES',NCNOWRIT,ncid)
          status=NF_INQ_VARID(ncid,'ASYM',id1)
          status=NF_INQ_VARID(ncid,'QEXT',id2)
          status=NF_INQ_VARID(ncid,'QSCT',id3)
          status=NF_INQ_VARID(ncid,'Q55E',id4)
c -----------------------------------------------------------------
c   read
c -----------------------------------------------------------------
          start(1)=1
          start(2)=1
          start(3)=1
          start(4)=1

          count(1)=17
          count(2)=15
          count(3)=6 
          count(4)=23
          count3(1)=17
          count3(2)=15
          count3(3)=23
 

          status=NF_GET_VARA_REAL(ncid,id1,start,count,ASY)
          status=NF_GET_VARA_REAL(ncid,id2,start,count,EXT)
          status=NF_GET_VARA_REAL(ncid,id3,start,count,SCA)
          status=NF_GET_VARA_REAL(ncid,id4,start,count3,QEX)

          status=NF_CLOSE('AMP_MIE_TABLES',NCNOWRIT,ncid)

          AMP_ASY = ASY 
          AMP_EXT = EXT 
          AMP_SCA = SCA 
          AMP_Q55 = QEX


      RETURN
      END SUBROUTINE SETUP_RAD
c -----------------------------------------------------------------

c -----------------------------------------------------------------
      SUBROUTINE GET_LW(NA,NS,AREFF,TQAB,Vf)
!@sum Calculation of LW absorption for AMP aerosols
!@sum Called in SETAER / RCOMPX
!@auth Susanne Bauer

      USE RADPAR, only: TRUQEX, TRSQEX, TRDQEX, TRUQSC, TRSQSC, TRDQSC
     *                  , REFU22, REFS25, REFD25
      INTEGER, intent(IN) :: NA,NS
      REAL*8,  intent(in) :: areff,Vf(6)
      REAL*8   TQEX(33),TQSC(33),TQAB(33),TQEX_S(33),TQSC_S(33)
      REAL*8   QXAERN(25),QSAERN(25)
      REAL*8   wts,wta
      INTEGER  n0,k,n,nn

c CORE  
       IF(NA==0) THEN
         TQAB(:)= 0d0
       ENDIF
                                      !                               1   2   3   4
      IF(NA > 0 .and. NA < 5) THEN    !    NA : Aerosol compositions SO4,SEA,NO3,OC
        N0=0
        IF(NA==2) N0=22
        IF(NA==3) N0=44
        IF(NA==4) N0=88

        DO 114 K=1,33
        DO 113 N=1,22
        NN=N0+N
        QXAERN(N)=TRUQEX(K,NN)
        QSAERN(N)=TRUQSC(K,NN)
  113   CONTINUE
        CALL SPLINE(REFU22,QXAERN,22,AREFF,TQEX(K),1.D0,1.D0,1)
        CALL SPLINE(REFU22,QSAERN,22,AREFF,TQSC(K),1.D0,1.D0,1)
        TQAB(K)=TQEX(K)-TQSC(K)
  114   CONTINUE

      ENDIF

                                      !                              5   
      IF(NA==5) THEN                  !   NA : Aerosol compositions BC
        DO 124 K=1,33
        QXAERN(:)=TRSQEX(K,:)    ! 1:25
        QSAERN(:)=TRSQSC(K,:)    ! 1:25
        CALL SPLINE(REFS25,QXAERN,25,AREFF,TQEX(K),1.D0,1.D0,1)
        CALL SPLINE(REFS25,QSAERN,25,AREFF,TQSC(K),1.D0,1.D0,1)
        TQAB(K)=TQEX(K)-TQSC(K)
  124   CONTINUE

      ENDIF

                                      !                              6
      IF(NA==6) THEN                  !   NA : Aerosol composition DST
        DO 134 K=1,33
        QXAERN(:)=TRDQEX(K,:)    ! 1:25
        QSAERN(:)=TRDQSC(K,:)    ! 1:25
        CALL SPLINE(REFD25,QXAERN,25,AREFF,TQEX(K),1.D0,1.D0,1)
        CALL SPLINE(REFD25,QSAERN,25,AREFF,TQSC(K),1.D0,1.D0,1)
        TQAB(K)=TQEX(K)-TQSC(K)
  134   CONTINUE

      ENDIF

c SHELL
         IF(NS > 0 .and. NS < 5) THEN    !    NS : Aerosol compositions SO4,SEA,NO3,OC
         N0=0
         IF(NS==2) N0=22
         IF(NS==3) N0=44
         IF(NS==4) N0=88


         DO K=1,33
         DO N=1,22
         NN=N0+N
         IF (NS==1) WTS=Vf(1)                      ! <- shell fraction of aerosol composition
         IF (NS==2) WTS=Vf(5)                      ! <- shell fraction of aerosol composition
         WTA=1.D0-WTS
         QXAERN(N)=TRUQEX(K,NN)
         QSAERN(N)=TRUQSC(K,NN)
         ENDDO
         CALL SPLINE(REFU22,QXAERN,22,AREFF,TQEX_S(K),1.D0,1.D0,1)
         CALL SPLINE(REFU22,QSAERN,22,AREFF,TQSC_S(K),1.D0,1.D0,1)
         TQAB(K)=(TQEX(K)*WTA + TQEX_S(K)*WTS)-(TQSC(K)*WTA + TQSC_S(K)*WTS)
         ENDDO

      ENDIF

      RETURN
      END SUBROUTINE GET_LW
c -----------------------------------------------------------------

