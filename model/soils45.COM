C**** INCLUDE FILE FOR SOILS45
      PARAMETER (NGM=6,NG=NGM+1,IMT=5)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/SOILS/   PR,HTPR,PRS,HTPRS,W(0:NGM,2),HT(0:NGM,2)
     & ,SNOWD(2),WS(0:NGM,2),TP(0:NGM,2),FICE(0:NGM,2),HOUR,COST,SINT
     & ,DZ(NGM),Q(IMT,NGM),QK(IMT,NGM),SL,FV,FB,ALAI,ATRG,ASHG,ALHG
     & ,ABETAD,ABETAV,ABETAT,ABETAP,ABETAB,ABETA,ACNA,ACNC,AEVAPW,AEVAPD
     & ,AEVAPB,ARUNS,ARUNU,AERUNS,AERUNU,ADIFS,AEDIFS,AEPC,AEPB,AEPP
     & ,AFHG,AF0DT,AF1DT,CNC,ZW(2),FD,FW,FM,VH,ALAIE,TBCS,TCS
     & ,SNOWM
     & ,ID,IGCM
      COMMON/SPAR/THETA(0:NG,2),THETS(0:NG,2),F(0:NG,2)
     & ,FH(0:NG,2),ZB(NG),ZC(NG),ZHTB,SB(33),SHW,SHI,FSN,ELH,SDSTNC,SHA
     & ,SHV,SHC(0:NG,2),ALAMW,ALAMI,ALAMA,ALAMSN,ALAMBR,ALAMS(IMT-1)
     & ,XKH(NG,2),XKHM(NG,2),FR(NG),RNF(2),RNFF(NG,2),H(0:NG,2)
     & ,XK(0:NG,2),XINFC(2),SNK(0:NG,2),SNKH(0:NG,2),SHCLC,D(0:NG,2)
     & ,XKU(0:NG,2),IST(33)
      COMMON/GPAR/       DT,TTOT,ZMAX,ZMIN,ALGDEL,XIN,AKAPPA,DTR,SGM
     & ,DTS,SHTPR,DTOUT,TFRZ,BETAD,HW,DR,DRS,RS,PRFR,C1,AL(2),XINC
     & ,PRE(2),BETADL(NGM),SPGSN,N,NMAX,LC
C     COMMON/WEIGHT/A(4,IMT-1),B(4,IMT-1),P(4,IMT-1)
      COMMON/FPAR/DTH,THETM(0:NG,2),THMAX(0:NG,2),HLM(0:64),
     &  XKLM(0:64,IMT-1),DLM(0:64,IMT-1),THM(0:64,IMT-1),ALPH0,NTH
      COMMON/EVPTR/QM1,Q1,QS,QB,QC,EVAP(2),EVAPW,EVAPD,EVAPS,EPB,EPC,
     &  QL,QLDT,SNOWF,SNOWFS
      REAL*8 FS(2),FS1(2,2)
      COMMON/SURFC/PRES,RHO,TS,VSM,CH,CD,SNHT,SRHT,TRHT,ZS,Z1,CDN,
     &     P1,PPBL,FCOR,FMAG,TG,T1,VG,SINA,EDDY,RI1,RIS,FPAD,FS,FS1,
     & T1G,TSG,IZ,IZX
      COMMON /CMHEAT/ THRM(2),SNSH(2),XLTH(2),ALTH(2)
      double precision, dimension(4,imt-1), parameter :: a=reshape(
     &     (/                  ! matric potential coefficients for
     &     .2514,  0.0136, -2.8319,  0.5958,   ! sand
     &     .1481,  1.8726,  0.1025, -3.6416,   ! loam
     &     .2484,  2.4842,  0.4583, -3.9470,   ! clay
     &     .8781, -5.1816, 13.2385,-11.9501/), ! peat
     &     (/4,imt-1/))
      double precision, dimension(4,imt-1), parameter :: b=reshape(
     &     (/                  ! conductivity coefficients for
     &     -0.4910, -9.8945,  9.7976, -3.2211,   ! sand
     &     -0.3238,-12.9013,  3.4247,  4.4929,   ! loam
     &     -0.5187,-13.4246,  2.8899,  5.0642,   ! clay
     &     -3.0848,  9.5497,-26.2868, 16.6930/), ! peat
     &     (/4,imt-1/))
      double precision, dimension(4,imt-1), parameter :: p=reshape(
     &     (/                  ! diffusivity coefficients for
     &     -0.1800, -7.9999,  5.5685, -1.8868,   ! sand
     &     -0.1000,-10.0085,  3.6752,  1.2304,   ! loam
     &     -0.1951, -9.7055,  2.7418,  2.0054,   ! clay
     &     -2.1220,  5.9983,-16.9824,  8.7615/), ! peat
     &     (/4,imt-1/))
      COMMON/SOILS2/SDATA(IM,JM,11*NGM+1),ALA(3,IM,JM),ACS(3,IM,JM)
     &   ,AFB(IM,JM),AFR(IM,JM,NGM),AVH(IM,JM)
      COMMON/SOILS3/GHDATA(IM,JM,4*NGM+5)
