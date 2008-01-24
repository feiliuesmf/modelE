c
c     save diagnostics for run of MODELE SCM

      SUBROUTINE  SCM_DIAG  


      USE RESOLUTION, only : LM
      USE MODEL_COM , only :  p,u,v,t,q,wm,NSTEPSCM,sige,sig,
     &                        I_TARG,J_TARG           
      USE CLOUDS_COM, only : SVLHX,SVLAT,RHSAV,CLDSAV,tauss,taumc,
     &                cldss,cldmc,csizmc,csizss
      USE SCMCOM
      USE SCMDIAG
      USE RAD_COM, only : srhr,trhr
      USE PBLCOM, only : TSAVG,WSAVG,QSAVG,USAVG,VSAVG
      USE FLUXES, only : GTEMP 
      USE DYNAMICS, only : PK 
      USE CONSTANT, only : SHA, GRAV   
      USE GEOM, only : dxyp 
      USE FILEMANAGER, only : openunit,closeunit
      

      IMPLICIT NONE


C     ALERT!!! LX=57 is for a 53-layer model (i.e., LX=LM+4)
C              LX=40          23-layer model 
C     If LX is changed, need to change data statements in 
C     BLOCK DATA RADPAR that initialize PLB and HLB
ccc   IMPLICIT REAL*8(A-H,O-Z)
CCC   PARAMETER (LX=57,LMOX=12*(1998-1881),MO3X=12*(2050-1950) )
c     PARAMETER (LX=40,LMOX=12*(1998-1881),MO3X=12*(2050-1850+1) )
c     REAL*4 TROAER,VDBCSU,TDUST,EPLMHC,UVLEAN   ! ,FVEG11,FOLGIZ     
c     REAL*4 ATAU98,SIZE98,HTF498,O3CLIM         !  only for offline use
c     CHARACTER*80 ATITLE,VTITLE,DTITLE,TITLE
C
c     COMMON/RADCM1/V6ALB(11,4,7),TAUWC0,TAUIC0,EPSCON,RO3COL,FULGAS(13)
c    A             ,FRAYLE,FALGAE,FCLDTR,FCLDSR,PTLISO,TLGRAD,FGOLDH(13)
c    B             ,WETTRA,WETSRA,FSXAER(5),FTXAER(5),FZSRA(6),FEMTRA(6)
c    C             ,KWVCON,KEEPAL,KEEPRH,KEEP10,KCNORM,KCLDEP,ICE012,NV
c    D             ,MADVES,MRELAY,MOZONE,KO3LON,NO3COL,NORMS0,KSOLAR
c    E             ,KTREND,NTRACE,ITR(8),NL,NLP,MLAT46,MLON72,LASTVC
c    X             ,HLB(LX),RHL(LX),TRACER(LX,8),AERTAU(LX),ZOICE
c    X             ,S00WM2,RATLS0,S0,XXJDAY,JYEARR,JDAYR,ISPARE(98)
c    X             ,DMOICE,DMLICE,TRSLCR,TRDFSL,TRUFSL
C     INPUT DATA
c     COMMON/RADCM2/
c    F              PLB(LX),        TLB(LX),TLT(LX),TLM(LX),U0GAS(LX,12)
c    G             ,ULGAS(LX,12),SHL(LX)
c    H             ,TAUWC(LX),TAUIC(LX),SIZEWC(LX),SIZEIC(LX),CLDEPS(LX)
c    I          ,POCEAN,PEARTH,POICE,PLICE,AGESN(3),SNOWE,SNOWOI,SNOWLI
c    J             ,TGO,TGE,TGOI,TGLI,TSL,WMAG,WEARTH,      FSPARE(998)
c    K                              ,COSZ,PVT(11),BXA(153),SRBXAL(15,2)
c    L                              ,JLAT,ILON
C     OUTPUT DATA
c    M             ,TRDFLB(LX),TRUFLB(LX),TRNFLB(LX),TRFCRL(LX)
c    N             ,SRDFLB(LX),SRUFLB(LX),SRNFLB(LX),SRFHRL(LX),SRSLHR
c    O             ,SRIVIS,SROVIS,PLAVIS,SRINIR,SRONIR,PLANIR,SRXATM(4)
c    P             ,SRDVIS,SRUVIS,ALBVIS,SRDNIR,SRUNIR,ALBNIR,FSRNFG(4)
c    Q             ,SRTVIS,SRRVIS,SRAVIS,SRTNIR,SRRNIR,SRANIR,FTRUFG(4)
c    R             ,TRDFGW,TRUFGW,TRUFTW,BTEMPW,              DTRUFG(4)
c    S             ,TRSLTS,TRSLTG,TRSLWV,TRSLBS,TTRUFG,LBOTCL,LTOPCL
c    X             ,Z0(12)
c
c-------------------------------------------------------------------------------
c
C
c     parameter (ntau=7,npres=7,ncol=100)
c     real*8 mnptopI,mntauI
c     COMMON /ISCCP/ fqI(ntau,npres),totcldareaI,mnptopI,
c    +               mntauI,bxtauI(ncol),bxptopI(ncol),
c    +               cldbdz(18,18),cldtdz(18,18)


      real*8 ARMT(LM),ARMQ(LM)
C
C 
C             Record Layout  
C 
C             NSTEPSCM            Time Stamp - note --- add date/time stamp also    
C             P                   Column Pressure (mb)
C             T          (LM)     Temperature (K) 
C             Q          (LM)     Specific Humidity (Kg/Kg)
C             TSAVG               Ts  -  Surface Air T (K) 
C             GTEMP(1,4,itarg,jtarg) Tskin  -  Sking temperature (C)
C             CLCVSS     (LM)     Cloud Cover SS (by area) 
C             CLCVMC     (LM)     Cloud Cover MC 
C             CLTHCK     (LM)     Cloud Thickness 
C             WMCOl      (LM)     Cloud Water Content (Kg/Kg) 
C             SVLHXCOL   (LM)     Liquid/Ice Flag (SS) save Latent Heats (j/Kg) 
C             SVLATCOL   (LM)     Liquid/Ice Flag (MC) save Latent Heats (j/Kg)
C             CSIZE      (LM,2)   Particle Size (10**-06m)     1=mc,2=ss 
C             EFFRAD     (LM)     Effective Radius (10**-06m)
C             CUMFLX     (LM)     Cumulus Mass Flux (mb/s) 
C             CUMHET     (LM)     Cumulus Heating  (10**14 W)
C             CUMOST     (LM)     Cumulus Moistening (10**14 W)
C             SRDFLBTOP           INC SW on Top of Atmos (W/m**2) 
C             SRNFLBTOP           SW ABSORBTION BELOW P0 (W/m**2)
C             TRUFLBTOP           UPward LW at P0  (W/m**2)
C             SRDFLBBOT           SW INC on Z0  (W/m**2)
C             SRNFLBBOT           SW ABS on Z0  (W/m**2)
C             TRUFLBBOT           Upward LW at Z0  (W/m**2)
C             TRDFLBBOT           LW INC on Z0  (W/m**2)
C             PRCSS               Precip - Large Scale SS (mm/hour) 
C             PRCMC               Precip - Convective (mm/hour)
C             EVPFLX              Evaporation Flux 
C             SHFLX               Sensible Heat Flux 
C             SOILMS              Soil Moisture 
C             SRFHRLCOL(1-LM)*COSZ1  SW Heating 
C             TRFCRLCOL (1-LM)    LW Heating 
C             TAUSSC     (LM)     Cloud Optical Thickness - SS 
C             TAUMCC     (LM)     Cloud Optical Thickness - MC 
C             SG_P       (LM)     Pressure at sigma layers (mb)
C             SG_T       (LM)     Arm T interpolated to Sig levels (K)
C             SG_Q       (LM)     ARM Q interpolated to Sig levels (g/KG)
C             SG_OMEGA   (LM)     ARM Omega (mb/hour)
C             APREC               ARM precipitation (mm/hour)
C             ALH                 ARM Latent heat upwards (W/m2)
C             ASH                 ARM Sensible heat upwards (W/m2)
C             AMEANPS             ARM Area mean Suface Pressure (mb)
C             ATSAIR              ARM Surface Temperature Air (C)
C             ATSKIN              ARM Surface Skin Temperature (C)
C             ARHSAIR             ARM Surface Relative Humidity (%)
C             SG_U       (LM)     ARM u-wind (m/s)
C             SG_V       (LM)     ARM v-wind (m/s)
C             SG_HOR_TMP_ADV(LM)  ARM Horizontal Temperature Advection(K/s) 
C             SG_VER_S_ADV(LM)    ARM Vertical S Advection (K/s)
C             SG_HOR_Q_ADV(LM)    ARM Horizontal Q Advection (kg/kg/s) 
C             SG_VER_Q_ADV(LM)    ARM Vertiacl Q Advection (kg/Kg/s)
C             CLSAV(LM)           SCM cloud fraction SS (by volume)
C             CLDFLG(LM)          SCM Cloud flag from Radia-outcome of Rand(0,1)
C             DWNFLX(LM)          SCM downdraft cloud mass flux (mb/s)
C             RHC(LM)             SCM Relative Humidity saved after Cloud 
c                                     routines (Qs(over water))
c             ALWP                ARM MWR cloud liquid water path (cm)
c             ADWDT               ARM d(Column_H2O)/dt (mm/hr)
c             ADWADV              ARM Column_H2O_Advection)_(mm/hr) 
c             ATLWUP              ARM TOA LW Up (W/m**2)
c             ATSWDN              ARM TOA SW Dn (W/m**2)
c             ATSWIN              ARM TOA SW IN (W/m**2)
c             SG_ARSCL(LM)        ARM  ARSCL CLD AMOUNT (%)
c             PRESAV(LM)        SCM Large Scale Precip by layer (kg/kg)
c             PREMC(LM)         SCM MSTCNV Precip by layer (kg/kg)
c             LHPSAV(LM)        SCM Phase of Precip for SS
c             LHPMC(LM)         SCM Phase of Precip for MC
c             WCUSCM(LM,2)        SCM Cumulus updraft speed for 2 two plume types
c             WCUDEEP(LM,2)       SCM  Cumulus updraft speed for deep convection - 2 plumes
c             PRCCDEEP(LM,2,LM)   SCM Precipitating convective condensate for Deep conv
c             NPRCCDEEP(LM,2,LM)  SCM Non-Precipitating conv condensate for Deep conv
c             TPALL(LM,2,LM)      SCM Plume Temperature for Deep Conv
c             MCCOND(LM,2,LM)     SCM convective condensate for deep and shallow  
c             PRCCGRP(LM,2,LM)    SCM deep convective condensate graupel
c             PRCCICE(LM,2,LM)    SCM deep convective condensate ICE
c

c             isccp record layout 
c
c             fqI(ntau,npres)     SCM (ISCCP diagnostics) fraction of model grid box
c                                     covered by each of the 49 ISCCP D level cloud 
c                                     types
c             totcldareaI         SCM (ISCCP diagnostics) fraction of model grid box
c                                     columns with cloud somewhere in them
c                                     (sum over all fqI)
c             mnptopI             SCM mean cloud top pressure (mb)
c             mntauI              SCM mean optical thickness
c             bxtauI(ncol)        SCM optical thickness in each column
c             bxptopI(ncol)       SCM cloud top pressure (mb) in each column
c             cldbdz(18,18)       SCM (cld diagnostics) cld base vs cld thickness
c             cldtdz(18,18)       SCM (cld diagnostics) cld top vs cld thickness
c            
C              
      real*8 TPRT(LM), QPRT(LM) ,TSURF, TSKIN, WMCOL(LM)    
      real*8 TDIFF,QDIFF
      real*8 PCOL, SVLHXCOL(LM),SVLATCOL(LM)    


      INTEGER L,LMIN,IC,IU    
      INTEGER IPLUM,IPL,IPLUMSV      

      if (NSTEPSCM.eq.0) then
          call openunit('scm.save.sige',iu,.true.,.false.)
          WRITE(iu) SIGE
          call closeunit(iu)
          call openunit('scm.output',iu_scm_diag,.true.,.false.)
      endif


      PCOL = P(I_TARG,J_TARG)
      TSURF = TSAVG(I_TARG,J_TARG)
      TSKIN = GTEMP(1,4,I_TARG,J_TARG)
      
      do L = 1,LM
         TPRT(L) = T(I_TARG,J_TARG,L)*PK(L,I_TARG,J_TARG) 
         QPRT(L) = Q(I_TARG,J_TARG,L)
         WMCOL(L) = WM(I_TARG,J_TARG,L)
         SVLHXCOL(L) = SVLHX(L,I_TARG,J_TARG)
         SVLATCOL(L) = SVLAT(L,I_TARG,J_TARG)
         CLCVSS(L) = CLDSS(L,I_TARG,J_TARG)
         CLCVMC(L) = CLDMC(L,I_TARG,J_TARG)
         CLSAV(L) = CLDSAV(L,I_TARG,J_TARG)   
         TAUSSC(L) = TAUSS(L,I_TARG,J_TARG)
         TAUMCC(L) = TAUMC(L,I_TARG,J_TARG)
      enddo      

      do L=1,LM
         WCUSCM(L,1) = 0.0
         WCUSCM(L,2) = 0.0
      enddo

      do ic=1,2
         IPLUM = 0
         do LMIN = 1,LM
            do L=1,LM
               if (WCUALL(L,ic,LMIN).ne.0.0d0) then
                   IPLUM = LMIN 
c                  write(iu_scm_prt,24) ic,lmin,iplum
c24                format(1x,
c    *               'WCUPLUM found   ic lmin iplum ',i5,i5,i5)
                   go to 25
               endif
            enddo
         enddo
25       continue
         if (IPLUM.gt.0) then
             do L=1,LM
                WCUSCM(L,ic) = WCUALL(L,ic,IPLUM) 
             enddo
         endif
      enddo
c     do ic=1,2
c        do l=1,LM 
c           write(iu_scm_prt,26) ic,l,wcudeep(l,ic)
c26         format(1x,'ic l wcudeep ',i5,i5,f10.4)
c        enddo
c     enddo

Ccccc before writing out diagnostics convert cumulus diagnostics
c     do L = 1,LM
c        write(iu_scm_prt,80) L,CUMFLX(L),DWNFLX(L),CUMHET(L),CUMOST(L)
c 80     format(1x,'L  flx dwn  het mst ',i5,4(f12.3))    
c        CUMFLX(L) = CUMFLX(L)/DTSRC
c        DWNFLX(L) = DWNFLX(L)/DTSRC
c        CUMHET(L) = CUMHET(L)*10.E-13*SHA*DXYP(J_TARG)/(GRAV*DTSRC)
c        CUMOST(L) = CUMOST(L)*10.E-13*SHA*DXYP(J_TARG)/(GRAV*DTSRC)
c     enddo

c     Use the hourly version of the ARM data to save as a diagnostic
      do L = 1,LM
         ARMT(L) = THR(L,NSTEPSCM+IKT)
         ARMQ(L) = QHR(L,NSTEPSCM+IKT)
      enddo
C
      WRITE(iu_scm_diag) NSTEPSCM,PCOL,TPRT,QPRT,TSURF,TSKIN,CLCVSS,
cc    WRITE(97) NSTEPSCM,PCOL,TPRT,QPRT,TSURF,TSKIN,CLCVSS,
     *           CLCVMC,CLTHCK,WMCOL,SVLHXCOL,SVLATCOL,CSIZE,EFFRAD,
     *           CUMFLX,CUMHET, CUMOST,SRDFLBTOP,SRNFLBTOP,TRUFLBTOP,
     *           SRDFLBBOT,SRNFLBBOT,TRUFLBBOT,TRDFLBBOT,PRCSS,PRCMC,
     *           EVPFLX,SHFLX,SOILMS,SRFHRLCOL,
     *           TRFCRLCOL,TAUSSC,TAUMCC,SG_P,ARMT,ARMQ,
     *           SG_OMEGA,APREC,ALH,ASH,AMEANPS,ATSAIR,ATSKIN,
     *           ARHSAIR,SG_U,SG_V,SG_HOR_TMP_ADV,
     *           SG_VER_S_ADV,SG_HOR_Q_ADV,SG_VER_Q_ADV,CLSAV,
     *           CLDFLG,DWNFLX,RHC,ALWP,ADWDT,ADWADV,ATLWUP,
     *           ATSWDN,ATSWIN,SG_ARSCL,PRESAV,PREMC,LHPSAV,LHPMC,
     *           WCUSCM,WCUDEEP,PRCCDEEP,NPRCCDEEP,TPALL,MCCOND,
     *           PRCCGRP,PRCCICE   

c     WRITE(3) TAU,P,fqI,totcldareaI,mnptopI,mntauI,bxtauI,bxptopI,
c    *         cldbdz,cldtdz
C 
      write(iu_scm_prt,99) NSTEPSCM
 99   format(//1x,'END OF TIME STEP      NSTEPSCM  ',i6)
      write(iu_scm_prt,100) NSTEPSCM,PCOl,TSURF,TSKIN 
100   format(1x,'NSTEP P tsurf tskin  ',i5,3(f10.4))
      write(iu_scm_prt,110) PRCSS,PRCMC
110   format(1x,'PRCSS MC ',2(f10.4))
      write(iu_scm_prt,120) SRNFLBBOT,SRNFLBTOP,SRDFLBBOT,SRDFLBTOP, 
     *              TRUFLBTOP,TRUFLBBOT,TRDFLBBOT
120   format(1x,'rad  sw nbot top dbot top ',4(f10.2),
     *       '    lw utop bot dbot ',3(f10.2)) 

c
      do L=1,LM
      
         write(iu_scm_prt,140) L,SG_P(L),TPRT(L),ARMT(L),
     +                 QPRT(L)*1000.0,ARMQ(L)*1000.0,    
     +                 WMCOL(L)*1000.0,TAUSSC(L),TAUMCC(L),
     +                 CLCVSS(L)*100.,CLCVMC(L)*100.,SG_ARSCL(L)
 140     format(1x,i3,f8.2,' T ',2(f7.2),'  Q ',2(f7.3),'  WM ',
     +          f7.4,'  tauss mc ',2(f7.3),'  camtss mc ',
     +          2(f7.3),' arscl ',f6.2)
      enddo 


      RETURN 

      END SUBROUTINE SCM_DIAG 
