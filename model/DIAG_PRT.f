      MODULE BDJ
!@sum  stores information for outputting composite zonal diagnostics
!@auth M. Kelley

      IMPLICIT NONE

!@param nj_out number of j-format output fields = 10
      integer, parameter :: nj_out=10
!@var units string containing output field units
      CHARACTER(LEN=50), DIMENSION(nj_out) :: UNITS_J_O
!@var lname string describing output field
      CHARACTER(LEN=50), DIMENSION(nj_out) :: LNAME_J_O
!@var sname string referencing output field in self-desc. output file
      CHARACTER(LEN=30), DIMENSION(nj_out) :: NAME_J_O
!@var stitle short title for print out
      CHARACTER(LEN=16), DIMENSION(nj_out) :: STITLE_J_O
!@var INUM_J_O,IDEN_J_O numerator and denominator for calculated J diags
      INTEGER, DIMENSION(nj_out) :: INUM_J_O, IDEN_J_O

      END MODULE BDJ

      SUBROUTINE J_TITLES
!@sum  J_TITLES calculated zonal diagnostics
!@auth M. Kelley/G. Schmidt
!@ver  1.0
      USE DAGCOM, only : j_srincp0,j_srnfp0,j_plavis,j_planir,j_srnfg
     *     ,j_srincg,j_albvis,j_albnir,j_srrvis,j_srrnir,j_sravis
     *     ,j_sranir
      USE BDJ
      IMPLICIT NONE
      INTEGER :: K
C**** These information are for J zonal budget calulated diagnostics
C**** Note that we assume that they are all ratios of two existing
C**** records.
c
      k = 0
c
      k = k + 1
      name_j_o(k) = 'plan_alb'
      lname_j_o(k) = ' TOTAL PLANETARY ALBEDO'
      units_j_o(k) = '%'
      stitle_j_o(k)= ' PLANETARY ALBDO'
      inum_j_o(k)  = J_SRNFP0
      iden_j_o(k)  = J_SRINCP0
c
      k = k + 1
      name_j_o(k) = 'plan_alb_vis'
      lname_j_o(k) = 'PLANETARY ALBEDO IN VISUAL'
      units_j_o(k) = '%'
      stitle_j_o(k)= ' PLAN ALB VISUAL'
      inum_j_o(k)  = J_PLAVIS
      iden_j_o(k)  = J_SRINCP0
c
      k = k + 1
      name_j_o(k) = 'plan_alb_nir'
      lname_j_o(k) = 'PLANETARY ALBEDO IN NEAR IR'
      units_j_o(k) = '%'
      stitle_j_o(k)= ' PLAN ALB NEARIR'
      inum_j_o(k)  = J_PLANIR
      iden_j_o(k)  = J_SRINCP0
c
      k = k + 1
      name_j_o(k) = 'surf_alb'
      lname_j_o(k) = 'GROUND ALBEDO'
      units_j_o(k) = '%'
      stitle_j_o(k)= ' SURFACE G ALBDO'
      inum_j_o(k)  = J_SRNFG
      iden_j_o(k)  = J_SRINCG
c
      k = k + 1
      name_j_o(k) = 'surf_alb_vis'
      lname_j_o(k) = 'GROUND ALBEDO IN VISUAL'
      units_j_o(k) = '%'
      stitle_j_o(k)= ' SURF ALB VISUAL'
      inum_j_o(k)  = J_ALBVIS
      iden_j_o(k)  = J_SRINCP0
c
      k = k + 1
      name_j_o(k) = 'surf_alb_nir'
      lname_j_o(k) = 'GROUND ALBEDO IN NEAR IR'
      units_j_o(k) = '%'
      stitle_j_o(k)= ' SURF ALB NEARIR'
      inum_j_o(k)  = J_ALBNIR
      iden_j_o(k)  = J_SRINCP0
c
      k = k + 1
      name_j_o(k) = 'atm_alb_vis'
      lname_j_o(k) = 'ATMOSPHERIC ALBEDO IN VISUAL'
      units_j_o(k) = '%'
      stitle_j_o(k)= '0ATMO ALB VISUAL'
      inum_j_o(k)  = J_SRRVIS
      iden_j_o(k)  = J_SRINCP0
c
      k = k + 1
      name_j_o(k) = 'atm_alb_nir'
      lname_j_o(k) = 'ATMOSPHERIC ALBEDO IN NEAR IR'
      units_j_o(k) = '%'
      stitle_j_o(k)= ' ATMO ALB NEARIR'
      inum_j_o(k)  = J_SRRNIR
      iden_j_o(k)  = J_SRINCP0
c
      k = k + 1
      name_j_o(k) = 'atm_abs_vis'
      lname_j_o(k) = 'ATMOSPHERIC ABSORPTION IN VISUAL'
      units_j_o(k) = 'W/M**2'
      stitle_j_o(k)= ' ATMO ABS VISUAL'
      inum_j_o(k)  = J_SRAVIS
      iden_j_o(k)  = J_SRINCP0
c
      k = k + 1
      name_j_o(k) = 'atm_abs_nir'
      lname_j_o(k) = 'ATMOSPHERIC ABSORPTION IN NEAR IR'
      units_j_o(k) = 'W/M**2'
      stitle_j_o(k)= ' ATMO ABS NEARIR'
      inum_j_o(k)  = J_SRANIR
      iden_j_o(k)  = J_SRINCP0

      RETURN
      END SUBROUTINE J_TITLES


      SUBROUTINE DIAGJ
C****
C**** THIS ROUTINE PRODUCES AREA WEIGHTED STATISTICS OF
C****
C   K   N
C****
C***1   1  SOLAR RADIATION INCIDENT ON PLANET (W/M**2)
C****
C**1A   2/1  PLANETARY ALBEDO (10**-2)
C**1B  72/1  PLANETARY ALBEDO VISUAL (10**-2)
C**1C  73/1  PLANETARY ALBEDO NEAR IR (10**-2)
C**1D   6/5  GROUND ALBEDO (10**-2)
C**1E  74/1  GROUND ALBEDO VISUAL (10**-2)
C**1F  75/1  GROUND ALBEDO NEAR IR (10**-2)
C**1G  76/1  ATMOSPHERIC ALBEDO VISUAL (10**-2)
C**1H  77/1  ATMOSPHERIC ALBEDO NEAR IR (10**-2)
C**1I  78/1  ATMOSPHERIC ABSORPTION VISUAL (10**-2)
C**1J  79/1  ATMOSPHERIC ABSORPTION NEAR IR (10**-2)
C****
C***2   2  SOLAR RADIATION ABSORBED BY PLANET (W/M**2)
C***3   3  SOLAR RADIATION ABSORBED BELOW PTOP (W/M**2)
C***4   4  SOLAR RADIATION ABSORBED BY ATMOSPHERE (W/M**2)
C***5   5  SOLAR RADIATION INCIDENT ON GROUND (W/M**2)
C***6   6  SOLAR RADIATION ABSORBED BY GROUND (W/M**2)
C***7   7  THERMAL RADIATION EMITTED BY PLANET (W/M**2)
C***8   8  THERMAL RADIATION AT PTOP (W/M**2)
C***9   9  THERMAL RADIATION EMITTED BY GROUND (W/M**2)
C**10  67  THERMAL RADIATION INCIDENT ON GROUND (W/M**2)
C****
C**11  55  BRIGHTNESS TEMPERATURE THROUGH WINDOW REGION (K-TF)
C****  10  NET RADIATION ABSORBED BY PLANET (W/M**2)
C****  11  NET RADIATION ABSORBED BELOW PTOP (W/M**2)
C****  12  NET RADIATION ABSORBED BY GROUND (W/M**2)
C****  81  NET CLEAR SKY RADIATION AT P0 W/M**2
C****  82  NET CLEAR SKY RADIATION AT TROPOPAUSE (WMO) W/M**2
C****  83  NET RADIATION AT TROPOPAUSE (WMO) W/M**2
C****  13  SENSIBLE HEAT FLUX INTO THE GROUND (W/M**2)
C****  14  EVAPORATION HEAT FLUX INTO THE GROUND (W/M**2)
C****  39  PRECIPITATION HEAT FLUX INTO THE GROUND (W/M**2)
C****
C**21  40  HEAT RUNOFF FROM FIRST GROUND LAYER (W/M**2)
C****  44  NET HEATING AT Z0 (W/M**2)
C****  42  CONDUCTION AT -Z1 (W/M**2)
C****  41  HEAT OF WATER OR ICE DUFFUSION AT -Z1 (W/M**2)
C****  16  NET HEATING AT -Z1 (W/M**2)
C****  15  CONDUCTION AT -Z1-Z2 (W/M**2)
C****  43  ENERGY OF ICE MELTING (OR TRANSPORTING) AT -Z1-Z2 (W/M**2)
C****  56  NET HEATING AT -Z1-Z2 (W/M**2)
C****  33  OCEAN TRANSPORT (W/M**2)
C****  48  HEAT RUNOFF THROUGH THE MIXED LAYER DEPTH (W/M**2)
C****
C**31  68  ENERGY DIFFUSION INTO THE THERMOCLINE (W/M**2)
C****  18  MEAN TEMPERATURE OF FIRST GROUND LAYER (.1 K-TF)
C****  17  MEAN TEMPERATURE OF SECOND GROUND LAYER (.1 K-TF)
C****  34  OCEAN TEMPERATURE AT THE MAXIMUM MIXED LAYER DEPTH
C****  23  SURFACE AIR TEMPERATURE (.1 K-TF)
C****  22  FIRST LAYER AIR TEMPERATURE (.1 K-TF)
C****  21  COMPOSITE AIR TEMPERATURE (.1 K-TF)
C****  35  STRATO TEMPERATURE CHANGE PER DEGREE LATITUDE (10**-2 K)
C****  36  TROPO TEMPERATURE CHANGE PER DEGREE LATITUDE (10**-2 K)
C****  24  STRATOSPHERIC STATIC STABILITY (10**-3 K/M)
C****
C**41  25  TROPOSPHERIC STATIC STABILITY (10**-3 K/M)
C****  26  STRATOSPHERIC RICHARDSON NUMBER (1)
C****  27  TROPOSPHERIC RICHARDSON NUMBER (1)
C****  28  STRATOSPHERIC ROSSBY NUMBER (1)
C****  29  TROPOSPHERIC ROSSBY NUMBER (1)
C****  37  L IN THE STRATOSPHERE (10**5 M)
C****  38  L IN THE TROPOSPHERE (10**5 M)
C****  64  GAM  (10**-3 K/M)
C****  65  GAMM  (10**-3 K/M)
C****  66  GAMC  (10**-3 K/M)
C****
C**51  57  INTEGRATED SUPER-SATURATION CLOUD COVER (10**-2)
C****  58  INTEGRATED MOIST CONVECTIVE CLOUD COVER (10**-2)
C****  59  INTEGRATED TOTAL CLOUD COVER (10**-2)
C****  60  MOIST CONVECTIVE CLOUD DEPTH (100 N)
C****  61  SUPER SATURATION PRECIPITATION (KG/M**2/86400 S)
C****  62  MOIST CONVECTIVE PRECIPITATION (KG/M**2/86400 S)
C****  20  PRECIPITATION (KG/M**2/86400 S)
C****  19  EVAPORATION (KG/M**2/86400 S)
C****  63  WATER CONTENT OF ATMOSPHERE (KG/M**2)
C****  54  WATER RUNOFF AT Z0 (KG/M**2/86400 S)
C****
C**61  45  WATER OR ICE DIFFUSION AT -Z1 (KG/M**2/86400 S)
C****  46  ICE MELTING (OR TRANSPORTING) AT -Z1-Z2 (KG/M**2/86400 S)
C****  47  WATER RUNOFF THROUGH MIXED LAYER DEPTH (KG/M**2/86400 S)
C****  49  WATER CONTAINED IN FIRST GROUND LAYER (KG/M**2)
C****  50  ICE CONTAINED IN FIRST GROUND LAYER (KG/M**2)
C****  51  WATER CONTAINED IN SECOND GROUND LAYER (KG/M**2)
C****  52  ICE CONTAINED IN SECOND GROUND LAYER (KG/M**2)
C****  53  SNOW DEPTH (KG/M**2)
C****  31  SNOW COVER (10**-2)
C****  30  OCEAN ICE COVER (10**-2)
C**71  69  SURFACE TYPE FRACTION (1)
C****
c      USE PRTCOM, only :
      USE CONSTANT, only :
     &     grav,rgas,sday,twopi,omega,kapa,bygrav
      USE MODEL_COM, only :
     &     im,jm,lm,fim,flice,
     &     dtsrc,fland,idacc,jhour,jhour0,jdate,jdate0,amon,amon0,
     &     jyear,jyear0,ls1,sige,itime,itime0,nday,xlabel,lrunid,ntype
      USE GEOM, only :
     &     dxyp,lat,LAT_DG
      USE DAGCOM, only :
     &     QCHECK,acc_period,aj,areg,jreg,kdiag,namreg,nreg,kaj,ia_j,
     &     j_ctopp,j_cdldep,j_pcldmc,j_srabs,j_srnfp0,j_srnfg,j_trnfp0,
     &     j_hsurf,j_trhdt,j_trnfp1,j_hatm,j_rnfp0,j_rnfp1,j_srnfp1,
     &     j_rhdt,j_hz1,j_edifs,j_f1dt,j_prcp,j_prcpss,j_prcpmc,j_hz0,
     &     j_shdt,j_evhdt,j_eprcp,j_erun1,j_hz2,j_f2dt,j_erun2,j_type,
     *     scale_j,stitle_j,lname_j,name_j,units_j
      USE BDJ
      IMPLICIT NONE
      SAVE
      DOUBLE PRECISION, DIMENSION(JM) ::
     &     CONTJ,CONTO,CONTL,CONTOI,S1,FLAT
      DOUBLE PRECISION, DIMENSION(NTYPE,JM) :: SPTYPE
      INTEGER, DIMENSION(JM) :: MLAT
      DOUBLE PRECISION, DIMENSION(NREG) :: SAREA
      INTEGER, DIMENSION(2) :: MHEM
      DOUBLE PRECISION, DIMENSION(2) :: FHEM
      INTEGER III
!@param NTYPE_OUT number of output budgets pages
      INTEGER, PARAMETER :: NTYPE_OUT=NTYPE+2  ! to include composites
C**** Expanded version of surfaces (including composites)
!@var TERRAIN name of surface type
      CHARACTER*16, DIMENSION(0:NTYPE_OUT) :: TERRAIN = (/
     *     '    (GLOBAL)','(OPEN OCEAN)',' (OCEAN ICE)','     (OCEAN)',
     *     '      (LAND)','  (LAND ICE)',' (OPEN LAKE)','  (LAKE ICE)',
     *     '     (LAKES)'/)
C**** Arrays needed for full output
      REAL*8, DIMENSION(JM+3,KAJ) :: BUDG
      CHARACTER*16, DIMENSION(KAJ) :: TITLEO
      CHARACTER*50, DIMENSION(KAJ) :: LNAMEO
      CHARACTER*30, DIMENSION(KAJ) :: SNAMEO
      CHARACTER*50, DIMENSION(KAJ) :: UNITSO
C**** weighting functions for surface types
      DOUBLE PRECISION, DIMENSION(0:NTYPE_OUT,NTYPE) ::
     *     WT=RESHAPE(          ! seperate types + composites
     *     (/1.,1.,0.,1.,0.,0.,0.,0.,0., 1.,0.,1.,1.,0.,0.,0.,0.,0.,
     *       1.,0.,0.,0.,1.,0.,0.,0.,0., 1.,0.,0.,0.,0.,1.,0.,0.,0.,
     *       1.,0.,0.,0.,0.,0.,1.,0.,1., 1.,0.,0.,0.,0.,0.,0.,1.,1./),
     *     (/NTYPE_OUT+1,NTYPE/) )
      INTEGER, DIMENSION(71) :: NDEX = (/
     *  1,2,3,4,5,6,7,8,9,67,           55,10,11,12,81,82,83,13,14,39,
     *     40,44,
     *  42,41,16,15,43,56,33,48,68,18,  17,34,23,22,21,35,36,24,25,26,
     *  27,28,29,37,38,64,65,66,57,58,  59,60,61,62,20,19,63,54,45,46,
     *  47,49,50,51,52,53,31,30,69/)

      DOUBLE PRECISION :: A1BYA2,A2BYA1,AMULT,BYA1,BYIACC,
     &     FGLOB,GSUM,GSUM2,GWT,HSUM,HSUM2,HWT,QDEN,QJ,QNUM,DAYS,WTX

      INTEGER :: I,IACC,INC,J,JH,JHEMI,JMHALF,JR,K,KA,KD1M,M,MD,N,ND,NN
     *     ,IT
      INTEGER :: IFIRST = 1
      IF (IFIRST.EQ.1) THEN
      IFIRST=0
C**** INITIALIZE CERTAIN QUANTITIES
      call j_titles
      KD1M=71  ! was 68, now includes J_TYPE,J_CLRTOA,J_CLRTPR,J_TOTTRP
               ! removed J_SWCOR
      INC=1+(JM-1)/24
      JMHALF=JM/2
      SAREA=0.
      DO J=1,JM
        S1(J)=IM
        DO I=1,IM
          JR=JREG(I,J)
          SAREA(JR)=SAREA(JR)+DXYP(J)
        END DO
      END DO
      S1(1)=1.
      S1(JM)=1.
      END IF
C**** OPEN PLOTTABLE OUTPUT FILE IF DESIRED
      IF (QCHECK)  ! the +1 is because types dimensioned 0:ntype_out
     & call open_j(trim(acc_period)//'.j'//XLABEL(1:LRUNID),ntype_out+1)

C**** CALCULATE THE DERIVED QUANTTIES
      BYA1=1./(IDACC(1)+1.D-20)
      A2BYA1=DFLOAT(IDACC(2))/DFLOAT(IDACC(1))
      A1BYA2=IDACC(1)/(IDACC(2)+1.D-20)
      DO JR=1,23 ! only 23 will fit on a green sheet
        AREG(JR,J_TRNFP0)=AREG(JR,J_HSURF)+A2BYA1*AREG(JR,J_TRHDT)/DTSRC
        AREG(JR,J_TRNFP1)=AREG(JR,J_HATM)+A2BYA1*AREG(JR,J_TRHDT)/DTSRC
        AREG(JR,J_SRABS) =AREG(JR,J_SRNFP0)-AREG(JR,J_SRNFG)
        AREG(JR,J_RNFP0) =AREG(JR,J_SRNFP0)+AREG(JR,J_TRNFP0)
        AREG(JR,J_RNFP1) =AREG(JR,J_SRNFP1)+AREG(JR,J_TRNFP1)
        AREG(JR,J_RHDT)  =A1BYA2*AREG(JR,J_SRNFG)*DTSRC+AREG(JR,J_TRHDT)
        AREG(JR,J_PRCP)  =AREG(JR,J_PRCPSS)+AREG(JR,J_PRCPMC)
        AREG(JR,J_HZ1)   =AREG(JR,J_EDIFS)+AREG(JR,J_F1DT)
        AREG(JR,J_HZ0)   =AREG(JR,J_RHDT)+AREG(JR,J_SHDT)+
     *       AREG(JR,J_EVHDT)+AREG(JR,J_EPRCP)-AREG(JR,J_ERUN1)
        AREG(JR,J_HZ2)   =AREG(JR,J_F2DT)+AREG(JR,J_ERUN2)
        AREG(JR,J_CTOPP) =IDACC(2)*SAREA(JR)*AREG(JR,J_CDLDEP)/
     *       (AREG(JR,J_PCLDMC)+1.D-20)
      END DO
      DO J=1,JM
      DO IT=1,NTYPE
        SPTYPE(IT,J) =AJ(J,J_TYPE,IT)*BYA1
        AJ(J,J_TRNFP0,IT)=AJ(J,J_HSURF,IT)+A2BYA1*AJ(J,J_TRHDT,IT)/DTSRC
        AJ(J,J_TRNFP1,IT)=AJ(J,J_HATM,IT) +A2BYA1*AJ(J,J_TRHDT,IT)/DTSRC
c       AJ(J,J_SWCOR,IT) =(1.-SRCOR)*AJ(J,J_SRNFG,IT)   !obsolete
        AJ(J,J_SRABS ,IT)=AJ(J,J_SRNFP0,IT)-AJ(J,J_SRNFG,IT)
        AJ(J,J_RNFP0 ,IT)=AJ(J,J_SRNFP0,IT)+AJ(J,J_TRNFP0,IT)
        AJ(J,J_RNFP1 ,IT)=AJ(J,J_SRNFP1,IT)+AJ(J,J_TRNFP1,IT)
        AJ(J,J_RHDT  ,IT)=A1BYA2*AJ(J,J_SRNFG,IT)*DTSRC+AJ(J,J_TRHDT,IT)
        AJ(J,J_PRCP  ,IT)=AJ(J,J_PRCPSS,IT)+AJ(J,J_PRCPMC,IT)
        AJ(J,J_HZ1   ,IT)=AJ(J,J_EDIFS,IT)+AJ(J,J_F1DT,IT)
        AJ(J,J_HZ0   ,IT)=AJ(J,J_RHDT,IT)+AJ(J,J_SHDT,IT)+
     *       AJ(J,J_EVHDT,IT)+AJ(J,J_EPRCP,IT)-AJ(J,J_ERUN1,IT)
        AJ(J,J_HZ2,   IT)=AJ(J,J_F2DT,IT)+AJ(J,J_ERUN2,IT)
C**** this is incorrect: ratios should not be summed over type
        AJ(J,J_CTOPP,IT) =IDACC(2)*SPTYPE(IT,J)*AJ(J,J_CDLDEP,IT)/
     *       (AJ(J,J_PCLDMC,IT)+1.D-20)
      END DO
      END DO
      DAYS=(Itime-Itime0)/DFLOAT(nday)
C****
C**** LOOP OVER SURFACE TYPES: 1 TO NTYPE
C****
      IF (KDIAG(1).GT.7) GO TO 510
      DO M=0,NTYPE_OUT
      WRITE (6,901) XLABEL
      WRITE (6,902) TERRAIN(M),JYEAR0,AMON0,JDATE0,JHOUR0,
     *  JYEAR,AMON,JDATE,JHOUR,ITIME,DAYS

      WRITE (6,903) (NINT(LAT_DG(J,1)),J=JM,INC,-INC)
      WRITE (6,905)
      DO K=1,KD1M
      N=NDEX(K)
      IACC=IDACC(IA_J(N))
C**** set weighting for denominator (different only for J_TYPE)
      MD=M
      IF (name_j(N).eq.'J_surf_type_frac') MD=0
      GSUM=0.
      GWT=0.
      DO JHEMI=1,2
        HSUM=0.
        HWT=0.
        DO JH=1,JMHALF
          J=(JHEMI-1)*JMHALF+JH
C**** Sum over types
          QJ=0
          WTX=0
          DO IT=1,NTYPE
            QJ =QJ +WT(M ,IT)*AJ(J,N,IT)
            WTX=WTX+WT(MD,IT)*SPTYPE(IT,J)
          END DO
          QJ=QJ*SCALE_J(N)
          WTX=WTX*IACC
          FLAT(J)=QJ/(WTX+1.D-20)
          MLAT(J)=NINT(FLAT(J))
          HSUM=HSUM+QJ*DXYP(J)*(FIM+1.-S1(J))
          HWT=HWT+WTX*DXYP(J)*(FIM+1.-S1(J))
        END DO
        FHEM(JHEMI)=HSUM/(HWT+1.D-20)
        GSUM=GSUM+HSUM
        GWT=GWT+HWT
      END DO
      FGLOB=GSUM/(GWT+1.D-20)
      IF (M.EQ.0) CALL KEYDJ (N,FGLOB,FHEM(2))
C**** Save BUDG for full output
      BUDG(1:JM,K)=FLAT(1:JM)
      BUDG(JM+1,K)=FHEM(1)
      BUDG(JM+2,K)=FHEM(2)
      BUDG(JM+3,K)=FGLOB
      TITLEO(K)=STITLE_J(N)
      LNAMEO(K)=LNAME_J(N)
      SNAMEO(K)=NAME_J(N)
      UNITSO(K)=UNITS_J(N)
C**** select output format depending on field name
      SELECT CASE (name_j(N)(3:len_trim(name_j(N))))
      CASE ('sstab_trop','SWCOR')
        WRITE (6,906) STITLE_J(N),FGLOB,FHEM(2),FHEM(1),
     *       (FLAT(J),J=JM,INC,-INC)
      CASE ('evap','prec','ross_num_strat','ross_num_trop'
     *       ,'ross_radius_strat','ross_radius_trop','ht_runoff_z0'
     *       ,'h2o_difs_z1','ice_thru_z1z2','h2o_runoff_mld'
     *       ,'ht_runoff_mld','wat_runoff_z0','ssprec','mcprec'
     *       ,'lapse_rate','lapse_rate_m','lapse_rate_c'
     *       ,'ht_thermocline')
        WRITE (6,911) STITLE_J(N),FGLOB,FHEM(2),FHEM(1),
     *       (FLAT(J),J=JM,INC,-INC)
      CASE ('ht_cond_z1z2','net_ht_z1','ocn_ht_trans','prec_ht_flx'
     *       ,'ht_wtr_difs_z1','ht_cond_z1','ht_ice_z1z2','nt_ht_z0'
     *       ,'net_ht_z1z2')
        WRITE (6,912) STITLE_J(N),FGLOB,FHEM(2),FHEM(1),
     *       (MLAT(J),J=JM,INC,-INC)
      CASE DEFAULT
        WRITE (6,907) STITLE_J(N),FGLOB,FHEM(2),FHEM(1),
     *       (MLAT(J),J=JM,INC,-INC)
      END SELECT
      IF (N.EQ.1) THEN
C**** CALCULATE AND PRINT ALBEDOS
      DO KA=1,nj_out
        NN=INUM_J_O(KA)
        ND=IDEN_J_O(KA)
        AMULT=1.
        IF (name_j_o(ka).eq.'plan_alb'.or.name_j_o(ka).eq.'surf_alb')
     *       AMULT=-1.
        GSUM=0.
        GSUM2=0.
        DO JHEMI=1,2
          HSUM=0.
          HSUM2=0.
          DO JH=1,JMHALF
            J=(JHEMI-1)*JMHALF+JH
C**** Sum over types
            QNUM=0
            QDEN=0
            DO IT=1,NTYPE
              QNUM=QNUM+WT(M,IT)*AJ(J,NN,IT)
              QDEN=QDEN+WT(M,IT)*AJ(J,ND,IT)
            END DO
            FLAT(J)=AMULT*(100.*QNUM/(QDEN+1.D-20)-50.)+50.
            MLAT(J)=FLAT(J)+.5
            HSUM=HSUM+QNUM*DXYP(J)*(FIM+1.-S1(J))
            HSUM2=HSUM2+QDEN*DXYP(J)*(FIM+1.-S1(J))
          END DO
          FHEM(JHEMI)=50.+AMULT*(100.*HSUM/(HSUM2+1.D-20)-50.)
          GSUM=GSUM+HSUM
          GSUM2=GSUM2+HSUM2
        END DO
        FGLOB=50.+AMULT*(100.*GSUM/(GSUM2+1.D-20)-50.)
        IF (M.EQ.0.AND.name_j_o(ka).eq.'plan_alb') CALL KEYDJA (FGLOB)
C**** Save BUDG for full output
      BUDG(1:JM,KA+KD1M)=FLAT(1:JM)
      BUDG(JM+1,KA+KD1M)=FHEM(1)
      BUDG(JM+2,KA+KD1M)=FHEM(2)
      BUDG(JM+3,KA+KD1M)=FGLOB
      TITLEO(KA+KD1M)=STITLE_J_O(KA)
      LNAMEO(KA+KD1M)=LNAME_J_O(KA)
      SNAMEO(KA+KD1M)=NAME_J_O(KA)
      UNITSO(KA+KD1M)=UNITS_J_O(KA)
C****
      WRITE (6,912) STITLE_J_O(KA),FGLOB,FHEM(2),FHEM(1),
     *     (MLAT(J),J=JM,INC,-INC)
      END DO
      END IF
      END DO
      WRITE (6,903) (NINT(LAT_DG(J,1)),J=JM,INC,-INC)
      WRITE (6,905)
      IF (QCHECK) CALL POUT_J(TITLEO,SNAMEO,LNAMEO,UNITSO,BUDG,KD1M+10
     *     ,TERRAIN(M),M+1) ! the +1 is because M starts at 0
      IF (KDIAG(1).GT.1) RETURN
      END DO
      if(qcheck) call close_j
  510 IF (KDIAG(1).GT.0.AND.KDIAG(1).NE.8) RETURN
C****
C**** PRODUCE REGIONAL STATISTICS
C****
      WRITE (6,901) XLABEL
      WRITE (6,902) '   (REGIONS)    ',JYEAR0,AMON0,JDATE0,JHOUR0,
     *  JYEAR,AMON,JDATE,JHOUR,ITIME,DAYS
      WRITE(6,918)(NAMREG(1,K),K=1,23),(NAMREG(2,K),K=1,23)
      DO K=1,KD1M
      N=NDEX(K)
      BYIACC=1./(IDACC(IA_J(N))+1.D-20)
      DO JR=1,23
        FLAT(JR)=AREG(JR,N)*SCALE_J(N)*BYIACC/SAREA(JR)
        MLAT(JR)=NINT(FLAT(JR))
      END DO
CF       DO 523 J=1,23
CF523    RBUDG(J,K)=FLAT(J)
C**** select output format based on field name
      SELECT CASE (name_j(N)(3:len_trim(name_j(N))))
      CASE ('evap','prec','sstab_trop','ross_num_strat','ross_num_trop',
     *     'ocn_lak_ice_frac','snow_cover','SWCOR','ross_radius_strat'
     *     ,'ross_radius_trop','ht_runoff_z0','h2o_difs_z1'
     *     ,'ice_thru_z1z2','h2o_runoff_mld','ht_runoff_mld','wat_g1'
     *     ,'ice_g1','snowdp','wat_runoff_z0','ssprec','mcprec','atmh2o'
     *     ,'lapse_rate','lapse_rate_m','lapse_rate_c','ht_thermocline')
        WRITE (6,910) STITLE_J(N),(FLAT(JR),JR=1,23)
      CASE DEFAULT
        WRITE (6,909) STITLE_J(N),(MLAT(JR),JR=1,23)
      END SELECT
      IF (N.EQ.1) THEN ! do albedo printout next
C**** CALCULATE AND PRINT ALBEDOS FOR REGIONAL STATISTICS
      DO KA=1,nj_out
        NN=INUM_J_O(KA)
        ND=IDEN_J_O(KA)
        AMULT=1.
        IF (name_j_o(ka).eq.'plan_alb'.or.name_j_o(ka).eq.'surf_alb')
     *       AMULT=-1.
        DO JR=1,23
          FLAT(JR)=AMULT*(100.*AREG(JR,NN)/(AREG(JR,ND)+1.D-20)-50.)+50.
          MLAT(JR)=FLAT(JR)+.5
        END DO
CF       DO 613 J=1,23
CF613    RBUDG(J,KA+KD1M)=FLAT(J)
        WRITE (6,909) STITLE_J_O(KA),(MLAT(JR),JR=1,23)
      END DO
      END IF
      END DO
      WRITE (6,905)
      WRITE(6,918)(NAMREG(1,K),K=1,23),(NAMREG(2,K),K=1,23)
      RETURN
C****
  901 FORMAT ('1',A)
  902 FORMAT ('0** BUDGETS',A16,'**   From:',I6,A6,I2,',  Hr',I3,
     *  6X,'To:',I6,A6,I2,', Hr',I3,'  Model-Time:',I9,5X,
     *  'Dif:',F7.2,' Days')
  903 FORMAT ('0',131('-')/20X,'G     NH    SH  ',24I4)
  904 FORMAT (A16,3I6,2X,24I4)
  905 FORMAT (1X,131('-'))
  906 FORMAT (A16,3F6.1,2X,24F4.1)
  907 FORMAT (A16,3F6.1,2X,24I4)
  909 FORMAT (A16,1X,23I5)
  910 FORMAT (A16,1X,23F5.1)
  911 FORMAT (A16,3F6.2,2X,24F4.1)
  912 FORMAT (A16,3F6.2,2X,24I4)
  918 FORMAT ('0',16X,23(1X,A4)/17X,23(1X,A4)/1X,131('-'))
      END SUBROUTINE DIAGJ


      MODULE BDjkl
!@sum  stores information for outputting lat-sigma/pressure diagnostics
!@auth M. Kelley
      IMPLICIT NONE

!@param njkl_out number of derived jk/jl-format output fields
      integer, parameter :: njkl_out=50

!@param names of derived jk/jl output fields
      INTEGER :: jkl_rad_cool,jkl_nt_qgpv,jkl_dudt_econv,
     *  jkl_psi_cp,jkl_dudt_epdiv,jkl_stdev_dp,
     *  jkl_dtempdt_econv,jkl_phi_amp_wave1,jkl_phi_phase_wave1,
     *  jkl_epflx_div,jkl_vt_dse_e,jkl_vt_lh_eddy,jkl_vt_se_eddy,
     *  jkl_tot_vt_se,jkl_psi_tem,jkl_epflx_v,
     *  jkl_nt_eqgpv,jkl_dyn_conv_eddy_geop,jkl_nt_sheat_e,
     *  jkl_dyn_conv_dse,jkl_seke,jkl_eke,
     *  jkl_nt_dse_se,jkl_nt_dse_e,jkl_tot_nt_dse,
     *  jkl_nt_lh_e,jkl_nt_see,jkl_tot_nt_se,
     *  jkl_nt_am_stand_eddy,jkl_nt_am_eddy,jkl_tot_nt_am,
     *  jkl_we_flx_nor,jkl_we_flx_div,jkl_refr_ind_wave1,
     *  jkl_del_qgpv,jkl_nt_lh_se,
     *  jkl_mcdrgpm10,jkl_mcdrgpm40,jkl_mcdrgpm20,
     *  jkl_sumdrg

!@var units string containing output field units
      CHARACTER(LEN=50), DIMENSION(njkl_out) :: UNITS_jkl
!@var lname string describing output field
      CHARACTER(LEN=50), DIMENSION(njkl_out) :: LNAME_jkl
!@var sname string referencing output field in self-desc. output file
      CHARACTER(LEN=30), DIMENSION(njkl_out) :: SNAME_jkl

      END MODULE BDjkl


      SUBROUTINE jkl_TITLES
      USE BDjkl
      IMPLICIT NONE
      INTEGER :: K
c
      k = 0
c
c derived JL-arrays
c
      k = k + 1 ! 01
      jkl_rad_cool = k
      sname_jkl(k) = 'rad_cool'
      lname_jkl(k) = 'TOTAL RADIATION COOLING RATE'
      units_jkl(k) = '10**13 WATTS/UNIT SIGMA'
      k = k + 1 ! 02
      jkl_nt_qgpv = k
      sname_jkl(k) = 'nt_qgpv'
      lname_jkl(k) = 'NORTHWARD TRANSPORT OF Q-G POT. VORTICITY'
      units_jkl(k) = '10**18 JOULES/DSIG'
      k = k + 1 ! 03
      jkl_dudt_econv = k
      sname_jkl(k) = 'dudt_eddy_conv'
      lname_jkl(k) = 'DU/DT BY EDDY CONVERGENCE (CP)'
      units_jkl(k) = '10**-6 M/S/S'
      k = k + 1 ! 04
      jkl_psi_cp = k
      sname_jkl(k) = 'psi_cp'
      lname_jkl(k) = 'STREAM FUNCTION (CP)'
      units_jkl(k) = '10**9 KILOGRAMS/SECOND'
      k = k + 1 ! 05
      jkl_dudt_epdiv = k
      sname_jkl(k) = 'dudt_epdiv'
      lname_jkl(k) = 'DU/DT BY ELIASSEN-PALM DIVERGENCE (CP)'
      units_jkl(k) = '10**-6 M/S/S'
      k = k + 1 ! 06
      jkl_stdev_dp = k
      sname_jkl(k) = 'stdev_dp'
      lname_jkl(k) = 'STANDARD DEVIATION OF PRESSURE DIFFERENCES'
      units_jkl(k) = 'MB'
      k = k + 1 ! 07
      jkl_dtempdt_econv = k
      sname_jkl(k) = 'dtempdt_eddy_conv'
      lname_jkl(k) = 'DTEMP/DT BY EDDY CONVERGENCE (CP)'
      units_jkl(k) = '10**-1 DEG-K/DAY'
      k = k + 1 ! 08
      jkl_phi_amp_wave1 = k
      sname_jkl(k) = 'phi_amp_wave1'
      lname_jkl(k) ='AMPLITUDE OF GEOPOTENTIAL HEIGHT FOR WAVE NUMBER 1'
      units_jkl(k) = 'METERS'
      k = k + 1 ! 09
      sname_jkl(k) = 'phi_amp_wave2'
      lname_jkl(k) ='AMPLITUDE OF GEOPOTENTIAL HEIGHT FOR WAVE NUMBER 2'
      units_jkl(k) = 'METERS'
      k = k + 1 ! 10
      sname_jkl(k) = 'phi_amp_wave3'
      lname_jkl(k) ='AMPLITUDE OF GEOPOTENTIAL HEIGHT FOR WAVE NUMBER 3'
      units_jkl(k) = 'METERS'
      k = k + 1 ! 11
      sname_jkl(k) = 'phi_amp_wave4'
      lname_jkl(k) ='AMPLITUDE OF GEOPOTENTIAL HEIGHT FOR WAVE NUMBER 4'
      units_jkl(k) = 'METERS'
      k = k + 1 ! 12
      jkl_phi_phase_wave1 = k
      sname_jkl(k) = 'phi_phase_wave1'
      lname_jkl(k) = 'PHASE OF GEOPOTENTIAL HEIGHT FOR WAVE NUMBER 1'
      units_jkl(k) = 'DEG WEST LONG'
      k = k + 1 ! 13
      sname_jkl(k) = 'phi_phase_wave2'
      lname_jkl(k) = 'PHASE OF GEOPOTENTIAL HEIGHT FOR WAVE NUMBER 2'
      units_jkl(k) = 'DEG WEST LONG'
      k = k + 1 ! 14
      sname_jkl(k) = 'phi_phase_wave3'
      lname_jkl(k) = 'PHASE OF GEOPOTENTIAL HEIGHT FOR WAVE NUMBER 3'
      units_jkl(k) = 'DEG WEST LONG'
      k = k + 1 ! 15
      sname_jkl(k) = 'phi_phase_wave4'
      lname_jkl(k) = 'PHASE OF GEOPOTENTIAL HEIGHT FOR WAVE NUMBER 4'
      units_jkl(k) = 'DEG WEST LONG'
      k = k + 1 ! 16
      jkl_epflx_div = k
      sname_jkl(k) = 'epflx_div'
      lname_jkl(k) = 'DIVERGENCE OF THE ELIASSEN-PALM FLUX'
      units_jkl(k) = '10**17 JOULES/UNIT SIGMA'
      k = k + 1 ! 17
      jkl_vt_dse_e = k
      sname_jkl(k) = 'vt_dse_e'
      lname_jkl(k) = 'VERT. TRANS. OF DRY STATIC ENERGY BY EDDIES (CP)'
      units_jkl(k) = '10**12 WATTS'
      k = k + 1 ! 18
      jkl_vt_lh_eddy = k
      sname_jkl(k) = 'vt_lh_eddy'
      lname_jkl(k) = 'VERTICAL TRANSPORT OF LATENT HEAT BY EDDIES (CP)'
      units_jkl(k) = '10**12 WATTS'
      k = k + 1 ! 19
      jkl_vt_se_eddy = k
      sname_jkl(k) = 'vt_se_eddy'
      lname_jkl(k) ='VERTICAL TRANSPORT OF STATIC ENERGY BY EDDIES (CP)'
      units_jkl(k) = '10**13 WATTS'
      k = k + 1 ! 20
      jkl_tot_vt_se = k
      sname_jkl(k) = 'tot_vt_se'
      lname_jkl(k) = 
     &    'TOTAL LARGE SCALE VERT. TRANS. OF STATIC ENERGY (CP)'
      units_jkl(k) = '10**14 WATTS'
      k = k + 1 ! 21
      jkl_psi_tem = k
      sname_jkl(k) = 'psi_tem'
      lname_jkl(k) = 'TRANSFORMED STREAM FUNCTION (CP)'
      units_jkl(k) = '10**9 KG/SEC'
      k = k + 1 ! 22
      jkl_epflx_v = k
      sname_jkl(k) = 'epflx_vert_cp'
      lname_jkl(k) = 'VERTICAL ELIASSEN-PALM FLUX (CP)'
      units_jkl(k) = '10**11 JOULES/METER'
c
c derived JK-arrays
c
      k = k + 1 ! 23
      jkl_nt_eqgpv = k
      sname_jkl(k) = 'nt_eddy_qgpv'
      lname_jkl(k) = 'NORTH. TRANS. OF EDDY Q-G POT. VORTICITY'
      units_jkl(k) = '10**-6 M/S**2'
      k = k + 1 ! 24
      jkl_dyn_conv_eddy_geop = k
      sname_jkl(k) = 'dyn_conv_eddy_geop'
      lname_jkl(k) = 'DYNAMIC CONVERGENCE OF EDDY GEOPOTENTIAL'
      units_jkl(k) = '.1 WATTS/M**2/DSIGMA'
      k = k + 1 ! 25
      jkl_nt_sheat_e = k
      sname_jkl(k) = 'nt_sheat_eddy'
      lname_jkl(k) = 'NORTH. TRANS. OF SENSIBLE HEAT BY EDDIES'
      units_jkl(k) = '10**14 WATTS/DSIGMA'
      k = k + 1 ! 26
      jkl_dyn_conv_dse = k
      sname_jkl(k) = 'dyn_conv_dse'
      lname_jkl(k) = 'DYNAMIC CONVERGENCE OF DRY STATIC ENERGY'
      units_jkl(k) = '10 WATTS/M**2/DSIGMA'
      k = k + 1 ! 27
      jkl_seke = k
      sname_jkl(k) = 'stand_eddy_ke'
      lname_jkl(k) = 'STANDING EDDY KINETIC ENERGY'
      units_jkl(k) = '10**4 JOULES/M**2/UNIT SIGMA'
      k = k + 1 ! 28
      jkl_eke = k
      sname_jkl(k) = 'eddy_ke'
      lname_jkl(k) = 'EDDY KINETIC ENERGY'
      units_jkl(k) = '10**4 JOULES/M**2/UNIT SIGMA'
      k = k + 1 ! 29
      jkl_nt_dse_se = k
      sname_jkl(k) = 'nt_dse_stand_eddy'
      lname_jkl(k) = 'NOR. TRANS. OF DRY STAT. ENERGY BY STAND. EDDIES'
      units_jkl(k) = '10**14 W/DSIG'
      k = k + 1 ! 30
      jkl_nt_dse_e = k
      sname_jkl(k) = 'nt_dse_eddy'
      lname_jkl(k) = 'NORTH. TRANS. OF DRY STATIC ENERGY BY EDDIES'
      units_jkl(k) = '10**14 WATTS/DSIG'
      k = k + 1 ! 31
      jkl_tot_nt_dse = k
      sname_jkl(k) = 'tot_nt_dse'
      lname_jkl(k) = 'TOTAL NORTH. TRANSPORT OF DRY STATIC ENERGY'
      units_jkl(k) = '10**15 WATTS/DSIG'
      k = k + 1 ! 32
      jkl_nt_lh_e = k
      sname_jkl(k) = 'nt_lh_e'
      lname_jkl(k) = 'NORTHWARD TRANSPORT OF LATENT HEAT BY EDDIES'
      units_jkl(k) = '10**13 WATTS/DSIG'
      k = k + 1 ! 33
      jkl_nt_see = k
      sname_jkl(k) = 'nt_se_eddy'
      lname_jkl(k) = 'NORTH.TRANSPORT OF STATIC ENERGY BY EDDIES'
      units_jkl(k) = '10**14 WATTS/DSIGMA'
      k = k + 1 ! 34
      jkl_tot_nt_se = k
      sname_jkl(k) = 'tot_nt_se'
      lname_jkl(k) = 'TOTAL NORTHWARD TRANSPORT OF STATIC ENERGY'
      units_jkl(k) = '10**15 WATTS/DSIGMA'
      k = k + 1 ! 35
      jkl_nt_am_stand_eddy = k
      sname_jkl(k) = 'nt_am_stand_eddy'
      lname_jkl(k) = 'NORTH. TRANS. OF ANG. MOMENTUM BY STAND. EDDIES'
      units_jkl(k) = '10**18 J/DSIG'
      k = k + 1 ! 36
      jkl_nt_am_eddy = k
      sname_jkl(k) = 'nt_am_eddy'
      lname_jkl(k) = 'NORTH. TRANS. OF ANG. MOMENTUM BY EDDIES'
      units_jkl(k) = '10**18 JOULES/DSIGMA'
      k = k + 1 ! 37
      jkl_tot_nt_am = k
      sname_jkl(k) = 'tot_nt_am'
      lname_jkl(k) = 'TOTAL NORTHWARD TRANSPORT OF ANG. MOMENTUM'
      units_jkl(k) = '10**19 JOULES/DSIG'
      k = k + 1 ! 38
      jkl_we_flx_nor = k
      sname_jkl(k) = 'we_flx_nor'
      lname_jkl(k) = 'NORTHWARD WAVE ENERGY FLUX'
      units_jkl(k) = '10**11 JOULES/METER/UNIT SIGMA'
      k = k + 1 ! 39
      jkl_we_flx_div = k
      sname_jkl(k) = 'we_flx_div'
      lname_jkl(k) = 'DIVERGENCE OF THE WAVE ENERGY FLUX'
      units_jkl(k) = '10**-5 M/S**2'
      k = k + 1 ! 40
      jkl_refr_ind_wave1 = k  !!!!! Refraction Inicies must be in order
      sname_jkl(k) = 'refr_ind_wave1'
      lname_jkl(k) = 'REFRACTION INDEX FOR WAVE NUMBER 1'
      units_jkl(k) = '10**-8 PER METER**2'
      k = k + 1 ! 41
      sname_jkl(k) = 'refr_ind_wave2'
      lname_jkl(k) = 'REFRACTION INDEX FOR WAVE NUMBER 2'
      units_jkl(k) = '10**-8 PER METER**2'
      k = k + 1 ! 42
      sname_jkl(k) = 'refr_ind_wave3'
      lname_jkl(k) = 'REFRACTION INDEX FOR WAVE NUMBER 3'
      units_jkl(k) = '10**-8 PER METER**2'
      k = k + 1 ! 43
      sname_jkl(k) = 'refr_ind_wave6'
      lname_jkl(k) = 'REFRACTION INDEX FOR WAVE NUMBER 6'
      units_jkl(k) = '10**-8 PER METER**2'
      k = k + 1 ! 44
      sname_jkl(k) = 'refr_ind_wave9'
      lname_jkl(k) = 'REFRACTION INDEX FOR WAVE NUMBER 9'
      units_jkl(k) = '10**-8 PER METER**2'
      k = k + 1 ! 45
      jkl_del_qgpv = k
      sname_jkl(k) = 'del_qgpv'
      lname_jkl(k) = 'Q-G POT. VORTICITY CHANGE OVER LATITUDES'
      units_jkl(k) = '10**-12 1/(SEC-M)'
      k = k + 1 ! 46
      jkl_nt_lh_se = k
      sname_jkl(k) = 'nt_lh_stand_eddy'
      lname_jkl(k) = 'N. TRANSPORT OF LATENT HEAT BY STAND. EDDIES'
      units_jkl(k) = '10**13 WATTS/DSIG'
c
c More derived arrays
c
      k=k+1      !47
      jkl_mcdrgpm10 = k
      sname_jkl(k) = 'dudt_mcdrgpm10' !ajl22+23
      lname_jkl(k) = 'DU/DT BY STRAT. MC DRAG  C=+/-10R'
      units_jkl(k) = '10**-6 M/S/S'
      k=k+1      !48
      jkl_mcdrgpm40 = k
      sname_jkl(k) = 'dudt_mcdrgpm40' !AJL24+25
      lname_jkl(k) = 'DU/DT BY STRAT. MC DRAG  C=+/-40R'
      units_jkl(k) = '10**-6 M/S/S'
      k=k+1      !49
      jkl_mcdrgpm20 = k
      sname_jkl(k) = 'dudt_mcdrgpm20' !AJL26+27
      lname_jkl(k) = 'DU/DT BY STRAT MC DRAG C=+/-20R'
      units_jkl(k) = '10**-6 M/S/S'
      k=k+1      !50
      jkl_sumdrg = k
      sname_jkl(k) = 'dudt_sumdrg' !AJL(18+20-27)
      lname_jkl(k) = 'ZONAL WIND CHANGE BY MTN+DEFORM+SHR+MC DRAG'
      units_jkl(k) = '10**-6 M/S/S'

      if (k .gt. njkl_out) then
        write (6,*) 'Increase njkl_out=',njkl_out,' to at least ',k
        stop 'jkl_TITLES: njkl_out too small'
      end if

      RETURN
      END SUBROUTINE jkl_TITLES


      SUBROUTINE DIAGJK
c      USE PRTCOM, only :
      USE CONSTANT, only :
     &     grav,rgas,kapa,sday,lhe,twopi,omega,sha,bygrav,tf
      USE MODEL_COM, only :
     &     im,jm,lm,fim, xlabel,lrunid,DO_GWDRAG,
     &     BYIM,DSIG,BYDSIG,DT,DTsrc,IDACC,IMH,LS1,NDAA,nidyn,
     &     PTOP,PMTOP,PSFMPT,SIG,SIGE,JHOUR
      USE GEOM, only :
     &     AREAG,BYDXYP,COSP,COSV,DLON,DXV,DXYP,DXYV,DYP,FCOR,RADIUS,WTJ
      USE DAGPCOM, only :
     &     PLM,PLE,P1000K,linect,jmby2
      USE DAGCOM
      USE BDjkl
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(JM) ::
     &     BYP,BYPV,BYDAPO,BYPDA,DXCOSV,ONESPO
     &    ,DXYPPO,DACOSV
      DOUBLE PRECISION, DIMENSION(JM+LM) :: ONES
      DOUBLE PRECISION, DIMENSION(JM,LM) :: AX,BX,CX,DX,VX,EX
      DOUBLE PRECISION, DIMENSION(JM,LM_REQ) :: ARQX
      DOUBLE PRECISION, DIMENSION(LM_REQ) :: BYDPS,BYPKS
      DOUBLE PRECISION, DIMENSION(2,IMH+1) :: CN
      DOUBLE PRECISION, DIMENSION(0:IMH) :: AN,BN

      DOUBLE PRECISION, DIMENSION(JM,LM,2) :: DSJK
      DOUBLE PRECISION, DIMENSION(2,LM,2) :: DSHEM
      DOUBLE PRECISION, DIMENSION(LM,2) :: DSGLOB
      COMMON/WORK5/DSJK,DSHEM,DSGLOB

      DOUBLE PRECISION, DIMENSION(LM) :: PM,PKM,PME
      DOUBLE PRECISION, DIMENSION(JM,2) :: PJ
      DOUBLE PRECISION, DIMENSION(JM,kgz+1,4) :: AMPLTD,PHASE
      INTEGER, PARAMETER, DIMENSION(5) :: MW=(/1,2,3,6,9/)
      DOUBLE PRECISION, PARAMETER :: ONE=1.

      INTEGER ::
     &     I,J,J0,J1,JH,JHEMI,K,K1,KDN,KM,KUP,L,LR,M,N,NS,NT,IX

      DOUBLE PRECISION ::
     &     BY100G,BYDP2,BYDPK,BYFSQ,BYIACN,BYIADA,BYIARD,BYIMDA,BYN,
     &     BYRCOS,BDN,BUP,
     &     DAM4,DXCVN,DXCVS,DALPHA,DE4TI,DP,DPG,DPH,DPTI,DTHDP,DTHETA,
     &     EL4QI,ELOFIM,   GBYRSQ,GSQ,
     &     PDN,PIG,PIH,PMK,PUP,PUTI,PVTI,SCALET,SCALES,SCALEV,SDDP,
     &     SKEI,SMALL,SN,SNAMI,SNDEGI,SNELGI,SQM,SQN,SZNDEG,SZNELG,
     &     THDN,THETA,THUP,TX,UDXN,UDXS,UX,WTKP1,WTN,WTNK,
     &     WTS,WTSK,XWON


C**** OPEN PLOTTABLE OUTPUT FILE IF DESIRED
      IF(QCHECK) call open_jl(trim(acc_period)//'.jk'//XLABEL(1:LRUNID))

C**** INITIALIZE CERTAIN QUANTITIES
      call jkl_titles

      IF (KDIAG(2).GE.8) GO TO 120
      XWON=TWOPI/(DLON*FIM)
      KM=LM
      BY100G=.01*BYGRAV
      DO 30 L=1,LM
      PKM(L)=PLM(L)**KAPA
      PME(L)=PSFMPT*SIGE(L)+PTOP
   30 PM(L)=PSFMPT*SIGE(L+1)+PTOP
      BYDPS(1)=1./(.5*PMTOP)
      BYDPS(2)=1./(.3*PMTOP)
      BYDPS(3)=1./(.2*PMTOP)
      BYPKS(1)=1./(.75*PMTOP)**KAPA
      BYPKS(2)=1./(.35*PMTOP)**KAPA
      BYPKS(3)=1./(.1*PMTOP)**KAPA
      ONES=1.
      DO 40 J=1,JM
      DXYPPO(J)=DXYP(J)
      ONESPO(J)=1.
c      BYDXYP(J)=1./DXYP(J)
      BYDAPO(J)=BYDXYP(J)
   40 CONTINUE
      DXYPPO(JM)=DXYP(JM)*FIM
      DXYPPO(1)=DXYP(1)*FIM
      ONESPO(1)=FIM
      ONESPO(JM)=FIM
      BYDAPO(1)=BYDAPO(1)*FIM
      BYDAPO(JM)=BYDAPO(JM)*FIM
      DO 50 J=2,JM
      DACOSV(J)=DXYV(J)*COSV(J)
      DXCOSV(J)=DXV(J)*COSV(J)
   50 CONTINUE
      LINECT=65
      WRITE (6,901)
      BYIACN=1./(IDACC(1)+1.D-20)
      BYIARD=1./(IDACC(2)+1.D-20)
      BYIADA=1./(IDACC(4)+1.D-20)
      BYIMDA=BYIADA*BYIM
      DO 52 J=1,JM
      PJ(J,1)=0
      PJ(J,2)=0
      DO 51 K=1,KM
      PJ(J,1)=PJ(J,1)+AJK(J,K,JK_DPA)
   51 PJ(J,2)=PJ(J,2)+AJK(J,K,JK_DPB)
      BYP(J)=1./(PJ(J,1)+1.D-20)
      BYPDA(J)=BYP(J)*BYDXYP(J)
   52 BYPV(J)=1./(PJ(J,2)+1.D-20)
C****
C**** INITIALIZE DELTA SIGMA IN PRESSURE COORDINATES
C****
      DO 60 J1=1,2
      J0=J1-1
      DO 60 K=1,KM
      DPG=0.
      PIG=0.
      DO 58 JHEMI=1,2
      DPH=0.
      PIH=0.
      DO 55 JH=1,jmby2
      J=(JHEMI-1)*(jmby2-J0)+JH+J0
      DSJK(J,K,J1)=AJK(J,K,J1)/(PJ(J,J1)+1.D-20)
      DPH=DPH+AJK(J,K,J1)*WTJ(J,2,J1)
   55 PIH=PIH+PJ(J,J1)*WTJ(J,2,J1)
      DSHEM(JHEMI,K,J1)=DPH/(PIH+1.D-20)
      DPG=DPG+DPH
   58 PIG=PIG+PIH
   60 DSGLOB(K,J1)=DPG/(PIG+1.D-20)
C****
C**** V-V* IS D/DP(V'TH'/DTH/DP) , DX=4*V'TH'/DTH/DP AT INTERFACES
C****
      DO 70 J=2,JM
      KDN=1
      DO 70 K=1,KM
      CX(J,K)=FIM*IDACC(4)
      AX(J,K)=AJK(J,K,JK_SHETH)/(AJK(J,K,JK_DPB)+1.D-20)
      KUP=K+1
      IF (K.EQ.KM) KUP=KM
      VX(J,K)=0.
      IF (AJK(J,K,JK_DPB).EQ.0.) GO TO 70
      IF (AJK(J,KDN,JK_DPB).EQ.0.) KDN=KDN+1
      VX(J,K)=AJK(J,K,JK_DPB)*(AJK(J,KUP,JK_SHETH)/AJK(J,KUP,JK_DPB)-
     *   AJK(J,KDN,JK_SHETH)/AJK(J,KDN,JK_DPB))/(PM(KUP)-PM(KDN)+.5*
     *   (AJK(J,KUP,JK_DPB)/AJK(J,KUP,JK_NPTSAVG)-
     &    AJK(J,KDN,JK_DPB)/AJK(J,KDN,JK_NPTSAVG)))
   70 KDN=K
C     DX(J,K)=AX*CX=4*(TRANSFORMED STREAM FUNCTION-STREAM FUNCTION)
   90 DO 95 J=2,JM
      DX(J,KM)=AX(J,KM)*CX(J,KM)
      DO 95 K=1,KM-1
      WTKP1=AJK(J,K,JK_DPB)/(AJK(J,K+1,JK_DPB)+AJK(J,K,JK_DPB)+1.D-20)
   95 DX(J,K)=(AX(J,K)*(1.-WTKP1)+AX(J,K+1)*WTKP1)*CX(J,K)
C****
C**** PROGNOSTIC QUANTITIES AT CONSTANT PRESSURE
C****
C**** # OF GRIDPOINTS, DELTA P, S.D. OF DELTA P
      n = JK_NPTSAVG
      SCALET=scale_jk(n)/idacc(ia_jk(n))
      CALL JLMAP(LNAME_JK(n),SNAME_JK(n),UNITS_JK(n),
     &     PLM,AJK(1,1,n),SCALET,ONES,ONES,LS1-1,1,JGRID_JK(n))
      n = JK_DPB
      SCALET=scale_jk(n)/idacc(ia_jk(n))
      CALL JLMAP(LNAME_JK(n),SNAME_JK(n),UNITS_JK(n),
     &     PLM,AJK(1,1,n),SCALET,ONES,ONES,LS1-1,2,JGRID_JK(n))
      DO 98 J=2,JM
      DO 98 K=1,LS1-1
      BYN=1./(AJK(J,K,JK_NPTSAVG)+1.D-10)
      AX(J,K)=0.
      SDDP=(AJK(J,K,JK_DPSQR)-AJK(J,K,JK_DPB)*AJK(J,K,JK_DPB)*BYN)*BYN
   98 IF (SDDP.GT.0.) AX(J,K)=SQRT(SDDP)
      n = jkl_stdev_dp
      CALL JLMAP(LNAME_jkl(n),SNAME_jkl(n),UNITS_jkl(n),
     &     PLM,AX,ONE,ONES,ONES,LS1-1,2,2)
C**** TEMPERATURE, HEIGHT, SPECIFIC AND RELATIVE HUMIDITY
      n = JK_TEMP
      CALL JKMAPS(LNAME_JK(n),SNAME_JK(n),UNITS_JK(n),
     &     PLM,AJK(1,1,n),ONE,BYP,ONES,KM,2,JGRID_JK(n),
     *     ASJL,BYIMDA,ONESPO,ONES)
      SCALES=BYIMDA*BY100G
      n = JK_HGHT
      CALL JKMAPS(LNAME_JK(n),SNAME_JK(n),UNITS_JK(n),
     &  PLM,AJK(1,1,n),BY100G,BYP,ONES,KM,2,JGRID_JK(n),
     *  ASJL(1,1,2),SCALES,ONESPO,ONES)
      SCALET=1.D5
      n = JK_Q
      CALL JKMAP(LNAME_JK(n),SNAME_JK(n),UNITS_JK(n),
     &    PLM,AJK(1,1,n),SCALET,BYP,ONES,KM,2,JGRID_JK(n))
      SCALET=100.
      n = JK_RH
      CALL JKMAP(LNAME_JK(n),SNAME_JK(n),UNITS_JK(n),
     &    PLM,AJK(1,1,n),SCALET,BYP,ONES,KM,2,JGRID_JK(n))
C**** LIQUID WATER CONTENT
      SCALET=1.D6
      n = JK_CLDH2O
      CALL JKMAP(LNAME_JK(n),SNAME_JK(n),UNITS_JK(n),
     &     PLM,AJK(1,1,n),SCALET,BYP,ONES,KM,2,JGRID_JK(n))
C**** U AND V WINDS, STREAM FUNCTION
      SCALET=10.
      n = JK_U
      CALL JKMAP(LNAME_JK(n),SNAME_JK(n),UNITS_JK(n),
     &    PLM,AJK(1,1,n),SCALET,BYPV,ONES,KM,2,JGRID_JK(n))
      SCALET=100.
      n = JK_V
      CALL JKMAP(LNAME_JK(n),SNAME_JK(n),UNITS_JK(n),
     &    PLM,AJK(1,1,n),SCALET,BYPV,ONES,KM,2,JGRID_JK(n))
      DO 100 J=2,JM
      AX(J,1)=AJK(J,1,JK_V)
      DO 100 K=2,KM
  100 AX(J,K)=AX(J,K-1)+AJK(J,K,JK_V)
      SCALET=100.D-9*XWON*BYIADA*BYGRAV
      n = jkl_psi_cp
      CALL JLMAP(LNAME_jkl(n),SNAME_jkl(n),UNITS_jkl(n),
     &     PM,AX,SCALET,DXV,ONES,KM,2,2)
      DO 110 K=1,KM
      DO 110 J=2,JM
  110 BX(J,K)=AX(J,K)+.25*DX(J,K)
      n = jkl_psi_tem
      CALL JLMAP(LNAME_jkl(n),SNAME_jkl(n),UNITS_jkl(n),
     &     PM,BX,SCALET,DXV,ONES,KM,2,2)
C**** VERTICAL WINDS
      SCALET=-1.D5*BYIMDA
      n = JK_VVEL
      CALL JLMAP(LNAME_JK(n),SNAME_JK(n),UNITS_JK(n),
     &     PME,AJK(1,1,n),SCALET,BYDXYP,ONES,KM,2,JGRID_JK(n))
C****
C**** CALCULATIONS FOR STANDING EDDIES
C****
  120 DO J=2,JM
      DO K=1,KM
        EX(J,K)=0.
        AX(J,K)=0.
        BX(J,K)=0.
        CX(J,K)=0.
      END DO
      END DO
      DO 170 J=2,JM
      DO 170 K=1,KM
      DPTI=0.
      PUTI=0.
      PVTI=0.
      DE4TI=0.
      EL4QI=0.
      SKEI=0.
      SNDEGI=0.
      SNELGI=0.
      SNAMI=0.
      DO 160 I=1,IM
      IF (AIJK(I,J,K,IJK_DP).EQ.0.) GO TO 160
      DPTI=DPTI+AIJK(I,J,K,IJK_DP)
      BYDPK=1./(AIJK(I,J,K,IJK_DP)+1.D-20)
      PUTI=PUTI+AIJK(I,J,K,IJK_U)
      PVTI=PVTI+AIJK(I,J,K,IJK_V)
      DE4TI=DE4TI+AIJK(I,J,K,IJK_DSE)
      EL4QI=EL4QI+AIJK(I,J,K,IJK_Q)
      SKEI=SKEI+(AIJK(I,J,K,IJK_U)*AIJK(I,J,K,IJK_U)
     *            +AIJK(I,J,K,IJK_V)*AIJK(I,J,K,IJK_V))*BYDPK
      SNDEGI=SNDEGI+(AIJK(I,J,K,IJK_DSE)*AIJK(I,J,K,IJK_V)*BYDPK)
      SNELGI=SNELGI+(AIJK(I,J,K,IJK_Q)*AIJK(I,J,K,IJK_V)*BYDPK)
      SNAMI=SNAMI+AIJK(I,J,K,IJK_U)*AIJK(I,J,K,IJK_V)*BYDPK
  160 CONTINUE
      AX(J,K)=SKEI-(PUTI*PUTI+PVTI*PVTI)/(DPTI+1.D-20)
      SZNDEG=DE4TI*PVTI/(DPTI+1.D-20)
      SZNELG=EL4QI*PVTI/(DPTI+1.D-20)
      EX(J,K)=SNELGI-SZNELG
      BX(J,K)=SNDEGI-SZNDEG
  170 CX(J,K)=SNAMI-PUTI*PVTI/(DPTI+1.D-20)
      IF (KDIAG(2).GE.8) RETURN
C**** STANDING EDDY, EDDY AND TOTAL KINETIC ENERGY
      SCALET=50.D-4*BYIMDA*BYGRAV
      n = jkl_seke
      CALL JKMAP(LNAME_jkl(n),SNAME_jkl(n),UNITS_jkl(n),
     &     PLM,AX,SCALET,ONES,ONES,KM,2,2)
      DO 200 K=1,KM
      DO 200 J=2,JM
  200 AX(J,K)=AJK(J,K,JK_TOTKE)-AJK(J,K,JK_ZMFKE)
      n = jkl_eke
      CALL JKMAP(LNAME_jkl(n),SNAME_jkl(n),UNITS_jkl(n),
     &     PLM,AX,SCALET,ONES,ONES,KM,2,2)
      n = jk_totke
      CALL JKMAP(LNAME_JK(n),SNAME_JK(n),UNITS_JK(n),
     &     PLM,AJK(1,1,n),SCALET,ONES,ONES,KM,2,JGRID_JK(n))
C**** POTENTIAL TEMPERATURE, POTENTIAL VORTICITY
      DO 205 LR=1,LM_REQ
      DO 205 J=1,JM
  205 ARQX(J,LR)=ASJL(J,LR,1)*BYIMDA*ONESPO(J)+TF
      N = JK_THETA
      CALL JKMAPS(LNAME_JK(N),SNAME_JK(N),UNITS_JK(N),
     &     PLM,AJK(1,1,N),P1000K,BYP,ONES,KM,2,JGRID_JK(N), 
     &     ARQX,P1000K,ONES,BYPKS)
      SCALET=1.D6*BYIMDA*P1000K
      N = JK_POTVORT
      CALL JLMAP(LNAME_JK(n),SNAME_JK(n),UNITS_JK(n),
     &     PLM,AJK(1,1,N),SCALET,BYDXYP,ONES,KM,2,JGRID_JK(n))
C****
C**** NORTHWARD TRANSPORTS AT CONSTANT PRESSURE
C****
C**** NORTHWARD TRANSPORT OF SENSIBLE HEAT BY EDDIES
      SCALET=25.D-14*XWON*SHA*BYIADA*BYGRAV
      DO 210 K=1,KM
      DO 210 J=2,JM
  210 AX(J,K)=AJK(J,K,JK_TOTNTSH)-AJK(J,K,JK_ZMFNTSH)
      N = jkl_nt_sheat_e
      CALL JKMAP(LNAME_jkl(n),SNAME_jkl(n),UNITS_jkl(n),
     &    PLM,AX,SCALET,DXV,ONES,KM,2,2)
C**** NORTHWARD TRANSPORT OF DRY STATIC ENERGY BY STANDING EDDIES,
C****   EDDIES, AND TOTAL
C**** Individual wave transports commented out. (gas - 05/2001)
      SCALET=25.D-14*XWON*BYIADA*BYGRAV
      N = jkl_nt_dse_se
      CALL JKMAP(LNAME_jkl(n),SNAME_jkl(n),UNITS_jkl(n),
     &     PLM,BX,SCALET,DXV,ONES,KM,2,2)
      DO 230 K=1,KM
      DO 230 J=2,JM
      AX(J,K)=SHA*(AJK(J,K,JK_TOTNTSH)-AJK(J,K,JK_ZMFNTSH))+
     &            (AJK(J,K,JK_TOTNTGEO)-AJK(J,K,JK_ZMFNTGEO))
  230 BX(J,K)=SHA*AJK(J,K,JK_TOTNTSH)+AJK(J,K,JK_TOTNTGEO)
      n = jkl_nt_dse_e
      CALL JKMAP(LNAME_jkl(n),SNAME_jkl(n),UNITS_jkl(n),
     &     PLM,AX,SCALET,DXV,ONES,KM,2,2)
      SCALET=SCALET*.1
      n = jkl_tot_nt_dse
      CALL JKMAP(LNAME_jkl(n),SNAME_jkl(n),UNITS_jkl(n),
     &     PLM,BX,SCALET,DXV,ONES,KM,2,2)
C**** NORTHWARD TRANSPORT OF LATENT HEAT BY STAND. EDDY, EDDIES AND TOTA
      DO 240 K=1,KM
      DO 240 J=2,JM
      DX(J,K)=AJK(J,K,JK_TOTNTLH)-AJK(J,K,JK_ZMFNTLH)
  240 AX(J,K)=AX(J,K)+LHE*DX(J,K)
      SCALET=25.D-13*XWON*LHE*BYIADA*BYGRAV
      n = jkl_nt_lh_se
      CALL JKMAP(LNAME_jkl(n),SNAME_jkl(n),UNITS_jkl(n),
     &     PLM,EX,SCALET,DXV,ONES,KM,2,2)
      n = jkl_nt_lh_e
      CALL JKMAP(LNAME_jkl(n),SNAME_jkl(n),UNITS_jkl(n),
     &     PLM,DX,SCALET,DXV,ONES,KM,2,2)
      SCALET=SCALET*.1
      n = jk_totntlh
      CALL JKMAP(LNAME_JK(n),SNAME_JK(n),UNITS_JK(n),
     &     PLM,AJK(1,1,n),SCALET,DXV,ONES,KM,2,JGRID_JK(n))
C**** NORTHWARD TRANSPORT OF STATIC ENERGY BY EDDIES AND TOTAL
      DO 245 K=1,KM
      DO 245 J=2,JM
  245 DX(J,K)=BX(J,K)+LHE*AJK(J,K,JK_TOTNTLH)
      SCALET=25.D-14*XWON*BYIADA*BYGRAV
      n = jkl_nt_see
      CALL JKMAP(LNAME_jkl(n),SNAME_jkl(n),UNITS_jkl(n),
     &     PLM,AX,SCALET,DXV,ONES,KM,2,2)
      SCALET=SCALET*.1
      n = jkl_tot_nt_se
      CALL JKMAP(LNAME_jkl(n),SNAME_jkl(n),UNITS_jkl(n),
     &     PLM,DX,SCALET,DXV,ONES,KM,2,2)
C**** NORTHWARD TRANSPORT OF KINETIC ENERGY
      SCALET=50.D-12*XWON*BYIADA*BYGRAV
      n = jk_totntke
      CALL JKMAP(LNAME_JK(n),SNAME_JK(n),UNITS_JK(n),
     &     PLM,AJK(1,1,n),SCALET,DXV,ONES,KM,2,JGRID_JK(n))
C**** NOR. TRANS. OF ANG. MOMENTUM BY STANDING EDDIES, EDDIES AND TOTAL
      SCALET=100.D-18*XWON*RADIUS*BYIADA*BYGRAV
      n = jkl_nt_am_stand_eddy
      CALL JKMAP(LNAME_jkl(n),SNAME_jkl(n),UNITS_jkl(n),
     &     PLM,CX,SCALET,DXCOSV,ONES,KM,2,2)
      DO 260 K=1,KM
      DO 260 J=2,JM
      CX(J,K)=AJK(J,K,JK_TOTNTMOM)-AJK(J,K,JK_ZMFNTMOM)
  260 DX(J,K)=AJK(J,K,JK_TOTNTMOM)+RADIUS*OMEGA*COSV(J)*AJK(J,K,JK_V)
      n = jkl_nt_am_eddy
      CALL JKMAP(LNAME_jkl(n),SNAME_jkl(n),UNITS_jkl(n),
     &     PLM,CX,SCALET,DXCOSV,ONES,KM,2,2)
      SCALET=.1*SCALET
      n = jkl_tot_nt_am
      CALL JKMAP(LNAME_jkl(n),SNAME_jkl(n),UNITS_jkl(n),
     &     PLM,DX,SCALET,DXCOSV,ONES,KM,2,2)
C****
C**** DYNAMIC CONVERGENCE OF ENERGY
C****
      SCALET=25.D-1*BYIMDA*BYGRAV
      DO 370 K=1,KM
C     CX(1,K)=-BX(2,K)*SCALET*DXV(2)
      CX(1,K)=0.
      CX(JM,K)=BX(JM,K)*SCALET*DXV(JM)
C     DX(1,K)=-(AJK(2,K,JK_TOTNTGEO)-AJK(2,K,JK_ZMFNTGEO))*SCALET*DXV(2)
      DX(1,K)=0.
      DX(JM,K)=(AJK(JM,K,JK_TOTNTGEO)-
     &     AJK(JM,K,JK_ZMFNTGEO))*SCALET*DXV(JM)
      DO 370 J=2,JM-1
      CX(J,K)=(BX(J,K)*DXV(J)-BX(J+1,K)*DXV(J+1))*SCALET
      DX(J,K)=((AJK(J,K,JK_TOTNTGEO)-AJK(J,K,JK_ZMFNTGEO))*DXV(J) -
     *  (AJK(J+1,K,JK_TOTNTGEO)-AJK(J+1,K,JK_ZMFNTGEO))*DXV(J+1))*SCALET
  370 CONTINUE
      SCALET=-100.D-1*BYIMDA*BYGRAV
      DO 380 K=1,KM-1
      DO 380 J=1,JM
      CX(J,K)=CX(J,K)-AJK(J,K,JK_TOTVTDSE)*SCALET
      CX(J,K+1)=CX(J,K+1)+AJK(J,K,JK_TOTVTDSE)*SCALET
      DX(J,K)=DX(J,K)-AJK(J,K,JK_VTGEOEDDY)*SCALET
  380 DX(J,K+1)=DX(J,K+1)+AJK(J,K,JK_VTGEOEDDY)*SCALET
      n = jkl_dyn_conv_dse
      CALL JKMAP(LNAME_jkl(n),SNAME_jkl(n),UNITS_jkl(n),
     &     PLM,CX,ONE,BYDXYP,ONES,KM,2,1)
      SCALET=100.
      n = jkl_dyn_conv_eddy_geop
      CALL JKMAP(LNAME_jkl(n),SNAME_jkl(n),UNITS_jkl(n),
     &     PLM,DX,SCALET,BYDXYP,ONES,KM,2,1)
C**** BAROCLINIC EKE GENERATION, P-K BY EDDY PRESSURE GRADIENT FORCE
      SCALET=50.E1*BYIMDA*RGAS*BYGRAV
      n = jk_barekegen
      CALL JKMAP(LNAME_JK(n),SNAME_JK(n),UNITS_JK(n),
     &     PLM,AJK(1,1,n),SCALET,BYDXYP,ONES,KM,2,JGRID_JK(n))
      SCALET=100.E1*BYIMDA/(GRAV*2.*DT)
      n = jk_p2kedpgf
      CALL JKMAP(LNAME_JK(n),SNAME_JK(n),UNITS_JK(n),
     &     PLM,AJK(1,1,n),SCALET,ONES,ONES,KM,2,JGRID_JK(n))
C****
C**** VERTICAL TRANSPORTS
C****
C**** VERTICAL TRANSPORT OF GEOPOTENTIAL ENERGY BY EDDIES
      SCALET=-100.D-12*XWON*BYIADA*BYGRAV
      n = jk_vtgeoeddy
      CALL JLMAP(LNAME_JK(n),SNAME_JK(n),UNITS_JK(n),
     &     PM,AJK(1,1,n),SCALET,ONES,ONES,KM-1,1,JGRID_JK(n))
C**** VERTICAL TRANSPORT OF DRY STATIC ENERGY BY EDDIES AND TOTAL
      DO 390 K=1,KM-1
      DO 390 J=1,JM
      AX(J,K)=AJK(J,K,JK_TOTVTDSE)-AJK(J,K,JK_ZMFVTDSE)
  390 BX(J,K)=AJK(J,K,JK_TOTVTLH)-AJK(J,K,JK_ZMFVTLH)
      SCALET=-100.D-12*XWON*BYIADA*BYGRAV
      n = jkl_vt_dse_e
      CALL JLMAP(LNAME_jkl(n),SNAME_jkl(n),UNITS_jkl(n),
     &     PM,AX,SCALET,ONES,ONES,KM-1,1,1)
      SCALET=SCALET*.01
      n = jk_totvtdse
      CALL JLMAP(LNAME_JK(n),SNAME_JK(n),UNITS_JK(n),
     &     PM,AJK(1,1,n),SCALET,ONES,ONES,KM-1,1,JGRID_JK(n))
C**** VERTICAL TRANSPORT OF LATENT HEAT BY EDDIES AND TOTAL
      SCALET=-100.D-12*XWON*LHE*BYIADA*BYGRAV
      n = jkl_vt_lh_eddy
      CALL JLMAP(LNAME_jkl(n),SNAME_jkl(n),UNITS_jkl(n),
     &     PM,BX,SCALET,ONES,ONES,KM-1,1,1)
      SCALET=SCALET*.1
      n = jk_totvtlh
      CALL JLMAP(LNAME_JK(n),SNAME_JK(n),UNITS_JK(n),
     &     PM,AJK(1,1,n),SCALET,ONES,ONES,KM-1,1,JGRID_JK(n))
C**** VERTICAL TRANSPORT OF STATIC ENERGY BY EDDIES AND TOTAL
      DO 420 K=1,KM-1
      DO 420 J=1,JM
      AX(J,K)=AX(J,K)+LHE*BX(J,K)
  420 BX(J,K)=AJK(J,K,JK_TOTVTDSE)+LHE*AJK(J,K,JK_TOTVTLH)
      SCALET=-100.D-13*XWON*BYIADA*BYGRAV
      n = jkl_vt_se_eddy
      CALL JLMAP(LNAME_jkl(n),SNAME_jkl(n),UNITS_jkl(n),
     &     PM,AX,SCALET,ONES,ONES,KM-1,1,1)
      SCALET=SCALET*.1
      n = jkl_tot_vt_se
      CALL JLMAP(LNAME_jkl(n),SNAME_jkl(n),UNITS_jkl(n),
     &     PM,BX,SCALET,ONES,ONES,KM-1,1,1)
C**** VERTICAL TRANSPORT OF KINETIC ENERGY
      SCALET=-12.5E-11*XWON*BYIADA*BYGRAV
      n = jk_totvtke
      CALL JLMAP(LNAME_JK(n),SNAME_JK(n),UNITS_JK(n),
     &     PM,AJK(1,1,n),SCALET,ONES,ONES,KM-1,1,JGRID_JK(n))
C**** VERTICAL TRANSPORT OF ANGULAR MOMENTUM BY LARE SCALET MOTIONS
      SCALET=-25.D-16*XWON*RADIUS*BYIADA*BYGRAV
      n = jk_vtameddy
      CALL JLMAP(LNAME_JK(n),SNAME_JK(n),UNITS_JK(n),
     &     PM,AJK(1,1,n),SCALET,COSV,ONES,KM-1,1,JGRID_JK(n))
      SCALET=1.D-2*SCALET
      n = jk_totvtam
      CALL JLMAP(LNAME_JK(n),SNAME_JK(n),UNITS_JK(n),
     &     PM,AJK(1,1,n),SCALET,COSV,ONES,KM-1,1,JGRID_JK(n))
C**** VERTICAL TRANSPORT OF POTENTIAL VORTICITY TOTAL AND BY EDDIES
      SCALET=-25.D-4*XWON*P1000K*BYIADA*BYGRAV
      n = jk_vtpv
      CALL JLMAP(LNAME_JK(n),SNAME_JK(n),UNITS_JK(n),
     &     PM,AJK(1,1,n),SCALET,BYDXYP,ONES,KM-1,1,JGRID_JK(n))
      n = jk_vtpveddy
      CALL JLMAP(LNAME_JK(n),SNAME_JK(n),UNITS_JK(n),PM,AJK(1,1,N),
     &     SCALET,BYDXYP,ONES,KM-1,1,JGRID_JK(n))
C**** NOR. TRANSPORT OF QUASI-GEOSTROPHIC POT. VORTICITY BY EDDIES
      DO 490 K=1,KM
      AX(1,K)=0.
      AX(JM,K)=0.
      DX(1,K)=0.
      DX(JM,K)=0.
      DO 490 J=2,JM-1
      AX(J,K)=((AJK(J,K,JK_TOTNTMOM)-AJK(J,K,JK_ZMFNTMOM))*
     &          DXCOSV(J)/(AJK(J,K,JK_DPB)+1.D-20)-
     *  (AJK(J+1,K,JK_TOTNTMOM)-AJK(J+1,K,JK_ZMFNTMOM))*DXCOSV(J+1)/
     *  (AJK(J+1,K,JK_DPB)+1.D-20))/COSP(J)
      DX(J,K)=FCOR(J)*(VX(J,K)+VX(J+1,K))/
     &     (AJK(J,K,JK_DPB)+AJK(J+1,K,JK_DPB)+1.D-20)
  490 CONTINUE
      DO 500 K=1,KM
      DO 500 J=2,JM-1
  500 AX(J,K)=AJK(J,K,JK_DPA)*(AX(J,K)+.25*DX(J,K))
      SCALET=1.D6
      n = jkl_nt_eqgpv
      CALL JKMAP(LNAME_jkl(n),SNAME_jkl(n),UNITS_jkl(n),
     &     PLM,AX,SCALET,BYPDA,ONES,KM,2,1)
C****
C**** ELIASSEN PALM FLUX:  NORTHWARD, VERTICAL, DIVERGENCE
C****
      SCALET=25.D-11*XWON*BYIADA*BYGRAV
      SMALL=1.D-20
      DO 510 K=1,KM
      AX(1,K)=0.
      DO 510 J=2,JM
      UX=AJK(J,K,JK_U)/(AJK(J,K,JK_DPB)+1.D-20)
      IF (ABS(UX).GE.SMALL) GO TO 510
      SN=+1.
      IF (UX.LT.0.) SN=-1.
      UX=SN*SMALL
  510 AX(J,K)=(AJK(J,K,JK_TOTNTGEO)-AJK(J,K,JK_ZMFNTGEO))/UX*DXV(J)
      n = jkl_we_flx_nor
      CALL JKMAP(LNAME_jkl(n),SNAME_jkl(n),UNITS_jkl(n),
     &     PLM,AX,SCALET,ONES,ONES,KM,2,2)
      SCALET=-100.D-11*XWON*BYIADA*BYGRAV
      DO 520 K=1,KM-1
      BX(1,K)=0.
      BX(JM,K)=0.
      DO 520 J=1,JM
      IF (J.NE.1.AND.J.NE.JM) GO TO 516  ! corrected 4-25-2000
      IF (J.EQ.1) UX=.5*(AJK(J+1,K,JK_U)+AJK(J+1,K+1,JK_U))/
     *     (AJK(J+1,K,JK_DPB)+AJK(J+1,K+1,JK_DPB)+1.D-20)
      IF (J.EQ.JM) UX=.5*(AJK(J,K,JK_U)+AJK(J,K+1,JK_U))/
     *     (AJK(J,K,JK_DPB)+AJK(J,K+1,JK_DPB)+1.D-20)
      GO TO 518
  516 UX=(AJK(J,K,JK_U)+AJK(J+1,K,JK_U)
     &   +AJK(J,K+1,JK_U)+AJK(J+1,K+1,JK_U))/
     *   (AJK(J,K,JK_DPB)+AJK(J+1,K,JK_DPB)+
     &    AJK(J,K+1,JK_DPB)+AJK(J+1,K+1,JK_DPB)+1.D-20)
  518 IF (ABS(UX).GE.SMALL) GO TO 520
      SN=+1.
      IF (UX.LT.0.) SN=-1.
      UX=SN*SMALL
  520 BX(J,K)=AJK(J,K,JK_VTGEOEDDY)/UX
      n = jkl_epflx_v
      CALL JLMAP(LNAME_jkl(n),SNAME_jkl(n),UNITS_jkl(n),
     &     PM,BX,SCALET,ONES,ONES,KM-1,1,1)
      DO 530 K=1,KM
      CX(1,K)=0.
      CX(JM,K)=0.
      DO 530 J=2,JM-1
  530 CX(J,K)=.25*(AX(J+1,K)-AX(J,K))
      DO 540 K=1,KM-1
      DO 540 J=1,JM
      CX(J,K)=CX(J,K)-BX(J,K)
  540 CX(J,K+1)=CX(J,K+1)+BX(J,K)
      SCALET=1.D5
      n = jkl_we_flx_div
      CALL JKMAP(LNAME_jkl(n),SNAME_jkl(n),UNITS_jkl(n),
     &     PLM,CX,SCALET,BYPDA,ONES,KM,2,1)
C****
C**** D/DY OF Q-G POTENTIAL VORTICITY AND REFRACTION INDICES
C****
C**** PRELIMINARIES:  VERTICAL DERIVATIVES AND N**2
      GSQ=GRAV*GRAV
      GBYRSQ=GRAV*GRAV/(RGAS*RGAS)
      IF (AJK(2,KM,JK_DPA).LT.1.D-20) GO TO 670  ! ISTART=4,...
      DO 600 J=1,JM
      K1=1
  580 IF (AJK(J,K1,JK_DPA).GT.1.D-20) GO TO 590
      AX(J,K1)=0.
      BX(J,K1)=0.
      DX(J,K1)=0.
      K1=K1+1
      GO TO 580
  590 KDN=K1
      PDN=PM(KDN)+.5*AJK(J,KDN,JK_DPA)/(AJK(J,KDN,JK_NPTSAVG1)+1.D-20)
      DO 600 K=K1,KM
      DP=AJK(J,K,JK_DPA)
      PMK=PM(K)+.5*AJK(J,K,JK_DPA)/(AJK(J,K,JK_NPTSAVG1)+1.D-20)
      KUP=K+1
      IF (K.EQ.KM) KUP=KM
      PUP=PM(KUP)+.5*AJK(J,KUP,JK_DPA)/(AJK(J,KUP,JK_NPTSAVG1)+1.D-20)
      DALPHA=(AJK(J,KUP,JK_TEMP)/(AJK(J,KUP,JK_DPA)+1.D-20)+TF)/PUP-
     *  (AJK(J,KDN,JK_TEMP)/(AJK(J,KDN,JK_DPA)+1.D-20)+TF)/PDN
      DTHETA=AJK(J,KUP,JK_THETA)/(AJK(J,KUP,JK_DPA)+1.D-20)-
     *  AJK(J,KDN,JK_THETA)/(AJK(J,KDN,JK_DPA)+1.D-20)
      THETA=AJK(J,K,JK_THETA)/(AJK(J,K,JK_DPA)+1.D-20)
      TX=AJK(J,K,JK_TEMP)/(AJK(J,K,JK_DPA)+1.D-20)+TF
      IF (ABS(DTHETA).GE.1.D-20) GO TO 595
      SN=+1.
      IF (DTHETA.LT.0.) SN=-1.
      DTHETA=SN*1.D-20
  595 DX(J,K)=DP*FCOR(J)*PMK*THETA*(PUP-PDN)/(TX*DTHETA*DXYP(J))
      AX(J,K)=DALPHA/(PUP-PDN-1.D-20)
C**** CALCULATE N**2 AT PRESSURE LATITUDES
      BX(J,K)=-DP*GSQ*PMK*DTHETA/(RGAS*TX*THETA*(PUP-PDN-1.D-20))
      KDN=K
  600 PDN=PMK
C**** CALCULATE  Q12 = (D(UDX) + F*DA)/DA
      DO 620 K=1,KM
      UDXS=0.
      DO 610 J=1,JM-1
      UDXN=AJK(J+1,K,JK_U)/(AJK(J+1,K,JK_DPB)+1.D-20)*DXV(J+1)
      CX(J,K)=(UDXS-UDXN+FCOR(J))/DXYP(J)
  610 UDXS=UDXN
      CX(JM,K)=(UDXS+FCOR(JM))/DXYP(JM)
C**** FIND DQ/DY = (Q12(J)-Q12(J-1)+Q3(J)-Q3(J-1))/DY
      DO 620 J=JM,2,-1
      DP=AJK(J,K,JK_DPB)
      AX(J,K)=DP*(CX(J,K)-CX(J-1,K) + (AX(J,K)-AX(J-1,K))*
     *  (DX(J,K)+DX(J-1,K))/
     &     (AJK(J,K,JK_DPA)+AJK(J-1,K,JK_DPA)+1.D-20))/DYP(3)
  620 CONTINUE
      SCALET=1.D12
      n = jkl_del_qgpv
      CALL JKMAP(LNAME_jkl(n),SNAME_jkl(n),UNITS_jkl(n),
     &     PLM,AX,SCALET,BYPV,ONES,KM,2,2)
C**** TERMS FOR THE REFRACTION INDEX EXPRESSION
      DO 640 J=2,JM
      BYFSQ=2.*DXYV(J)*DXYV(J)/(FCOR(J-1)*FCOR(J-1)+FCOR(J)*FCOR(J))
      DO 640 K=1,KM
      BYDP2=1./(AJK(J-1,K,JK_DPA)+AJK(J,K,JK_DPA)+1.D-20)
      TX=BYDP2*(AJK(J-1,K,JK_TEMP)+AJK(J,K,JK_TEMP))+TF
      DX(J,K)=GBYRSQ/(TX*TX)
      SQN=BYDP2*(BX(J-1,K)+BX(J,K))
      CX(J,K)=SQN*BYFSQ
      UX=AJK(J,K,JK_U)
      IF (ABS(UX).GE.1.D-20) GO TO 635
      SN=+1.
      IF (UX.LT.0.) SN=-1.
      UX=SN*1.D-20
  635 AX(J,K)=AX(J,K)/UX
  640 CONTINUE
C**** COMBINE TERMS, PRINT OUT REFRACTION INDICES
      SCALET=1.D8
      IX = jkl_refr_ind_wave1-1
      DO 660 M=1,5
      SQM=MW(M)*MW(M)
      DO 650 J=2,JM
      BYRCOS=1./(RADIUS*RADIUS*COSV(J)*COSV(J))
      DO 650 K=1,KM
      DP=AJK(J,K,JK_DPB)
  650 BX(J,K)=DP*(CX(J,K)*(AX(J,K)-SQM*BYRCOS)-.25*DX(J,K))
  660 CALL JKMAP(LNAME_jkl(M+IX),SNAME_jkl(M+IX),UNITS_jkl(M+IX),
     &    PLM,BX,SCALET,BYPV,ONES,KM,2,2)
  670 CONTINUE
C**** SKIP REMAINING MAPS IF DATA NOT AVAILABLE
      IF (AJK(1,1,JK_EPFLXNCP).NE.0.) GO TO 799
C****
C**** CHANGE OF THE MEAN FIELDS OF WIND AND TEMPERATURE
C****
C**** WIND: RATE OF CHANGE, ADVECTION, EDDY CONVERGENCE
      IF (IDACC(4).LE.1) GO TO 730
      SCALET=1.D6/((IDACC(4)-1)*(DTsrc*NDAA+DT+DT))
      n = JK_TOTDUDT
      CALL JLMAP(LNAME_JK(n),SNAME_JK(n),UNITS_JK(n),
     &     PLM,AJK(1,1,n),SCALET,ONES,ONES,KM,2,JGRID_JK(n))
  730 CONTINUE
C**** Depending on whether EP fluxes have been specially calcualted
C**** output full or approximate version
      IF (KEP.gt.0) THEN
        CALL EPFLXP
      ELSE ! these are not very good
      AX=0.
      BX=0.
      DO 720 K=2,KM-1
      DO 720 J=2,JM-1
      AX(J,K)=((AJK(J,K,JK_TOTNTMOM)-
     &  AJK(J,K,JK_ZMFNTMOM))*DXV(J)-(AJK(J+1,K,JK_TOTNTMOM)-
     *  AJK(J+1,K,JK_ZMFNTMOM))*DXV(J+1))/
     &        (AJK(J,K,JK_DPA)*DXYP(J)+1.D-20)+
     *  .125*((AJK(J,K,JK_VTAMEDDY)-AJK(J,K-1,JK_VTAMEDDY))/
     &        (AJK(J,K,JK_DPB)*DXYV(J)+1.D-20)+
     * (AJK(J+1,K,JK_VTAMEDDY)-AJK(J+1,K-1,JK_VTAMEDDY))/
     &        (AJK(J+1,K,JK_DPB)*DXYV(J+1)+1.D-20))
      BX(J,K)=(AJK(J+1,K,JK_EPFLXNCP)*DXCOSV(J+1)-
     &         AJK(J,K,JK_EPFLXNCP)*DXCOSV(J))/
     *  (DXYP(J)*COSP(J))+.5*(AJK(J,K-1,JK_EPFLXVCP)-
     &                        AJK(J,K,JK_EPFLXVCP)+
     *  AJK(J+1,K-1,JK_EPFLXVCP)-AJK(J+1,K,JK_EPFLXVCP))/(PM(K-1)-PM(K))
  720 CONTINUE
        SCALET=1.D6*BYIADA
        n = jk_dudtmadv
        CALL JLMAP(LNAME_JK(n),SNAME_JK(n),UNITS_JK(n),
     &       PLM,AJK(1,1,n),SCALET,ONES,ONES,KM,2,JGRID_JK(n))
        SCALET=1.D6
        n = jkl_dudt_econv
        CALL JLMAP(LNAME_jkl(n),SNAME_jkl(n),UNITS_jkl(n),
     &       PLM,AX,SCALET,ONES,ONES,KM,2,1)
C**** WIND: TRANSFORMED ADVECTION, LAGRANGIAN CONVERGENCE (DEL.F)
        SCALET=1.D6*BYIADA
        n = jk_dudttem
        CALL JLMAP(LNAME_JK(n),SNAME_JK(n),UNITS_JK(n),
     &       PLM,AJK(1,1,n),SCALET,ONES,ONES,KM,2,JGRID_JK(n))
        n = jkl_dudt_epdiv
        CALL JLMAP(LNAME_jkl(n),SNAME_jkl(n),UNITS_jkl(n),
     &       PLM,BX,SCALET,ONES,ONES,KM,2,1)
      END IF
C**** WIND: DU/DT BY STRAT. DRAG -  MTN, DEFORM., SHEAR ...
      SCALET=scale_jl(jl_dudfmdrg)/idacc(ia_jl(jl_dudfmdrg))
      if (DO_GWDRAG) then
      n = jl_dumtndrg
      CALL JLMAP (LNAME_JL(n),SNAME_JL(n),UNITS_JL(n),
     *     PLM,AJL(1,1,n),SCALET,ONES,ONES,LM,2,JGRID_JL(n))
      n = jl_dudfmdrg
      CALL JLMAP (LNAME_JL(n),SNAME_JL(n),UNITS_JL(n),
     *     PLM,AJL(1,1,n),SCALET,ONES,ONES,LM,2,JGRID_JL(n))
      n = jl_dushrdrg
      CALL JLMAP (LNAME_JL(n),SNAME_JL(n),UNITS_JL(n),
     *     PLM,AJL(1,1,n),SCALET,ONES,ONES,LM,2,JGRID_JL(n))
      DX=0.
      DO 740 L=1,LM
      DO 740 J=1,JM
      AX(J,L)=AJL(J,L,jl_dumcdrgm10)+AJL(J,L,jl_dumcdrgp10)
      BX(J,L)=AJL(J,L,jl_dumcdrgm40)+AJL(J,L,jl_dumcdrgp40)
  740 DX(J,L)=AJL(J,L,jl_dumcdrgm20)+AJL(J,L,jl_dumcdrgp20)
      CALL JLMAP (LNAME_jkl(jkl_mcdrgpm10),SNAME_jkl(jkl_mcdrgpm10),
     &     UNITS_jkl(jkl_mcdrgpm10),PLM,AX,SCALET,ONES,ONES,LM,2,2)
      CALL JLMAP (LNAME_jkl(jkl_mcdrgpm40),SNAME_jkl(jkl_mcdrgpm40),
     &     UNITS_jkl(jkl_mcdrgpm40),PLM,BX,SCALET,ONES,ONES,LM,2,2)
      CALL JLMAP (LNAME_jkl(jkl_mcdrgpm20),SNAME_jkl(jkl_mcdrgpm20),
     &     UNITS_jkl(jkl_mcdrgpm20),PLM,DX,SCALET,ONES,ONES,LM,2,2)
C**** DU/DT BY STRAT. DRAG - TOTAL
      DO 745 L=1,LM
      DO 745 J=1,JM
  745 AX(J,L)=AJL(J,L,jl_dumtndrg)+AJL(J,L,jl_dushrdrg)+
     *   (AX(J,L)+BX(J,L)+DX(J,L)) + AJL(J,L,jl_dudfmdrg)
      n = jkl_sumdrg
      CALL JLMAP (LNAME_jkl(n),SNAME_jkl(n),UNITS_jkl(n),
     *     PLM,AX,SCALET,ONES,ONES,LM,2,2)
      end if
C**** CHANGE IN EAST WIND BY STRATOSPHERIC DIFFUSION
!     CALL JLMAP (28,PLM,AJL(1,1,32),SCALE,ONES,ONES,LM,2,2)
C**** DU/DT BY SDRAG
      n = jl_dudtsdrg
      SCALET=scale_jl(n)/idacc(ia_jl(n))
      CALL JLMAP(LNAME_JL(n),SNAME_JL(n),UNITS_JL(n),
     &     PLM,AJL(1,1,n),SCALET,ONES,ONES,LM,2,JGRID_JL(n))
C**** TEMPERATURE: RATE OF CHANGE, ADVECTION, EDDY CONVERGENCE
      IF (IDACC(4).GT.1) then
      SCALET=1.D1*SDAY/((IDACC(4)-1)*(DTsrc*NDAA+DT+DT))
      n = JK_TOTDTDT
      CALL JLMAP(LNAME_JK(n),SNAME_JK(n),UNITS_JK(n),
     &     PLM,AJK(1,1,n),SCALET,ONES,PKM,KM,2,JGRID_JK(n))
      end if
      SCALET=1.D1*SDAY*BYIADA
      n = JK_DTDTMADV
      CALL JLMAP(LNAME_JK(n),SNAME_JK(n),UNITS_JK(n),
     &     PLM,AJK(1,1,n),SCALET,ONES,PKM,KM,2,JGRID_JK(n))
      cx = 0.
      do k=2,km-1
      do j=2,jm-1
        CX(J,K)=.25*(
     &     (AJK(J,  K,JK_TOTNTSH)-AJK(J,  K,JK_ZMFNTSH))*DXV(J)-
     &     (AJK(J+1,K,JK_TOTNTSH)-AJK(J+1,K,JK_ZMFNTSH))*DXV(J+1))/
     &     (AJK(J,K,JK_DPA)*DXYP(J)+1.D-20)+
     *     (AJK(J,K,JK_EDDVTPT)  -AJK(J,K-1,JK_EDDVTPT))/
     &     (PM(K-1)-PM(K))*BYIADA*PKM(K)
      end do
      end do
      n = jkl_dtempdt_econv
      SCALET=1.D1*SDAY
      CALL JLMAP(LNAME_jkl(n),SNAME_jkl(n),UNITS_jkl(n),
     &      PLM,CX,SCALET,ONES,ONES,KM,2,1)
C**** TEMPERATURE: TRANSFORMED ADVECTION
      n = JK_DTDTTEM
      SCALET=1.D1*SDAY*BYIADA
      CALL JLMAP(LNAME_JK(n),SNAME_JK(n),UNITS_JK(n),
     &     PLM,AJK(1,1,n),SCALET,ONES,PKM,KM,2,JGRID_JK(n))
C**** CHANGE IN TEMPERATURE BY STRATOSPHERIC DRAG
      n = jl_dtdtsdrg
      SCALET=scale_jl(n)/idacc(ia_jl(n))
      CALL JLMAP(LNAME_JL(n),SNAME_JL(n),UNITS_JL(n),
     &     PLM,AJL(1,1,n),SCALET,ONES,PKM,KM,2,JGRID_JL(n))
C**** CHANGE IN TEMPERATURE BY DYNAMICS
      n = JL_DTDYN
      SCALET=scale_jl(n)/idacc(ia_jl(n))
      CALL JLMAP(LNAME_JL(n),SNAME_JL(n),UNITS_JL(n),
     &     PLM,AJL(1,1,n),SCALET,ONESPO,PKM,KM,2,JGRID_JL(n))
  799 CONTINUE

C****
C**** Transplanted from DIAGJL
C****
      LINECT=65
C**** MASS FLUX MOIST CONVECTION
      n = JL_MCMFLX
      SCALET=scale_jl(n)/idacc(ia_jl(n))
      CALL JLMAP(LNAME_JL(n),SNAME_JL(n),UNITS_JL(n),
     &     PLE,AJL(1,1,n),SCALET,DXYPPO,ONES,LM-1,1,JGRID_JL(n))

C****
C**** RADIATION, CONDENSATION AND CONVECTION
C****
C**** SOLAR AND THERMAL RADIATION HEATING
      DO J=1,JM ! temporary redefinition of byp
        BYP(J)=IDACC(4)/(APJ(J,1)+1.D-20)
      ENDDO

      ns = JL_SRHR
      nt = JL_TRCR
      SCALET=scale_jl(ns)/idacc(ia_jl(ns))
      SCALES=scale_sjl(3)/idacc(ia_sjl(3))
      CALL JLMAPS(LNAME_JL(ns),SNAME_JL(ns),UNITS_JL(ns),
     &     PLM,AJL(1,1,ns),SCALET,BYP,BYDSIG,LM,2,JGRID_JL(ns),
     *     ASJL(1,1,3),SCALES,ONESPO,BYDPS)
      SCALET=scale_jl(nt)/idacc(ia_jl(nt))
      SCALES=scale_sjl(4)/idacc(ia_sjl(4))
      CALL JLMAPS(LNAME_JL(nt),SNAME_JL(nt),UNITS_JL(nt),
     &     PLM,AJL(1,1,nt),SCALET,BYP,BYDSIG,LM,2,JGRID_JL(nt),
     *     ASJL(1,1,4),SCALES,ONESPO,BYDPS)
      DO J=1,JM
        DO LR=1,LM_REQ
          ARQX(J,LR)=ASJL(J,LR,3)+ASJL(J,LR,4)
        ENDDO
        DO L=1,LM
          AX(J,L)=AJL(J,L,ns)+AJL(J,L,nt)
        ENDDO
      ENDDO
      SCALET=-1.D-13*XWON*BYIARD
      SCALES=SCALET*PSFMPT
      n = jkl_rad_cool
      CALL JLMAPS(LNAME_jkl(n),SNAME_jkl(n),UNITS_jkl(n),
     &     PLM,AX,SCALET,DXYPPO,BYDSIG,LM,1,1,  
     &     ARQX,SCALES,DXYPPO,BYDPS)

C**** TOTAL, SUPER SATURATION, AND CONVECTIVE CLOUD COVER
      n = JL_TOTCLD
      SCALET=scale_jl(n)/idacc(ia_jl(n))
      CALL JLMAP(LNAME_JL(n),SNAME_JL(n),UNITS_JL(n),
     &     PLM,AJL(1,1,n),SCALET,ONESPO,ONES,LM,2,JGRID_JL(n))
      n = JL_SSCLD
      SCALET=scale_jl(n)/idacc(ia_jl(n))
      CALL JLMAP(LNAME_JL(n),SNAME_JL(n),UNITS_JL(n),
     &     PLM,AJL(1,1,n),SCALET,ONESPO,ONES,LM,2,JGRID_JL(n))
      n = JL_MCCLD
      SCALET=scale_jl(n)/idacc(ia_jl(n))
      CALL JLMAP(LNAME_JL(n),SNAME_JL(n),UNITS_JL(n),
     &     PLM,AJL(1,1,n),SCALET,ONESPO,ONES,LM,2,JGRID_JL(n))
C**** TURBULENT KINETIC ENERGY
      n = JL_TRBKE
      SCALET=scale_jl(JL_TRBKE)/idacc(ia_jl(JL_TRBKE))
      CALL JLMAP(LNAME_JL(JL_TRBKE),SNAME_JL(JL_TRBKE),
     &     UNITS_JL(JL_TRBKE),PLM,AJL(1,1,JL_TRBKE),
     &     SCALET,ONES,ONES,LM,2,JGRID_JL(JL_TRBKE))
C**** HEATING BY LARGE SCALET COND., MOIST CONVECTION AND TURBULENCE
      n = JL_SSHR
      SCALET=scale_jl(n)/idacc(ia_jl(n))
      CALL JLMAP(LNAME_JL(n),SNAME_JL(n),UNITS_JL(n),
     &     PLM,AJL(1,1,n),SCALET,DXYPPO,BYDSIG,LM,1,JGRID_JL(n))
      n = JL_TRBHR
      SCALET=scale_jl(n)/idacc(ia_jl(n))
      CALL JLMAP(LNAME_JL(n),SNAME_JL(n),UNITS_JL(n),
     &     PLM,AJL(1,1,n),SCALET,DXYPPO,ONES,LM,1,JGRID_JL(n))
      n = JL_TRBDLHT
      SCALET=scale_jl(n)/idacc(ia_jl(n))
      CALL JLMAP(LNAME_JL(n),SNAME_JL(n),UNITS_JL(n),
     &     PLM,AJL(1,1,n),SCALET,DXYPPO,BYDSIG,LM,1,JGRID_JL(n))
      n = jl_mcdlht
      SCALET=scale_jl(n)/idacc(ia_jl(n))
      CALL JLMAP(LNAME_JL(n),SNAME_JL(n),UNITS_JL(n),
     &     PLM,AJL(1,1,n),SCALET,DXYPPO,BYDSIG,LM,1,JGRID_JL(n))
      n = JL_MCHEAT
      SCALET=scale_jl(n)/idacc(ia_jl(n))
      CALL JLMAP(LNAME_JL(n),SNAME_JL(n),UNITS_JL(n),
     &     PLM,AJL(1,1,n),SCALET,DXYPPO,ONES,LM,1,JGRID_JL(n))
      n = JL_MCDRY
      SCALET=scale_jl(n)/idacc(ia_jl(n))
      CALL JLMAP(LNAME_JL(n),SNAME_JL(n),UNITS_JL(n),
     &     PLM,AJL(1,1,n),SCALET,DXYPPO,ONES,LM,1,JGRID_JL(n))
C****
C**** ENERGY
C****
C**** AVAILABLE POTENTIAL ENERGY
      n = JL_APE
      SCALET=scale_jl(n)/idacc(ia_jl(n))
      CALL JLMAP(LNAME_JL(n),SNAME_JL(n),UNITS_JL(n),
     &      PLM,AJL(1,1,n),SCALET,ONES,ONES,LM,2,JGRID_JL(n))
C****
C**** NORTHWARD TRANSPORTS
C****
C**** NOR. TRANSPORT OF QUASI-GEOSTROPHIC POT. VORTICITY BY EDDIES
      DO 366 L=1,LM
      CX(1,L)=0.
      CX(2,L)=DXCOSV(2)*(AJL(2,L,JL_TOTNTMOM)-AJL(2,L,JL_ZMFNTMOM))+
     &     .25*FIM*FCOR(2)*COSP(2)*(AJL(2,L,JL_47)+AJL(3,L,JL_47))
      DO 364 J=3,JM-1
      DAM4=DXCOSV(J)*(AJL(J,L,JL_TOTNTMOM)-AJL(J,L,JL_ZMFNTMOM))
      CX(J,L)=DAM4+.25*FIM*FCOR(J)*COSP(J)*
     &     (AJL(J,L,JL_47)+AJL(J-1,L,JL_47))
      CX(J-1,L)=CX(J-1,L)-DAM4
  364 CONTINUE
      CX(JM-1,L)=CX(JM-1,L)-DXCOSV(JM)*(AJL(JM,L,JL_TOTNTMOM)-
     &     AJL(JM,L,JL_ZMFNTMOM))
      CX(JM,L)=0.
  366 CONTINUE
      SCALET=25.D-18*XWON*BYIADA*RADIUS*BYGRAV * 0. ! AJL47 not done
      n = jkl_nt_qgpv
      CALL JLMAP(LNAME_jkl(n),SNAME_jkl(n),UNITS_jkl(n),
     &     PLM,CX,SCALET,ONES,ONES,LM,1,1)
C****
C**** VERTICAL TRANSPORTS
C****
C**** VERTICAL TRANSPORT OF ANGULAR MOMENTUM BY SMALL SCALE MOTIONS
      n = JL_DAMDC
      SCALET=scale_jl(n)/idacc(ia_jl(n))
      CALL JLMAP(LNAME_JL(n),SNAME_JL(n),UNITS_JL(n),
     &     PLM,AJL(1,1,n),SCALET,DACOSV,ONES,LM,1,JGRID_JL(n))
      n = JL_DAMMC
      SCALET=scale_jl(n)/idacc(ia_jl(n))
      CALL JLMAP(LNAME_JL(n),SNAME_JL(n),UNITS_JL(n),
     &     PLM,AJL(1,1,n),SCALET,DACOSV,BYDSIG,LM,1,JGRID_JL(n))
C****
C**** MERIDIONAL LUNES
C****
C**** U, V AND W VELOCITY FOR EAST PACIFIC
      n = JL_UEPAC
      SCALET=scale_jl(n)/idacc(ia_jl(n))
      CALL JLMAP(LNAME_JL(n),SNAME_JL(n),UNITS_JL(n),
     &     PLM,AJL(1,1,n),SCALET,ONES,ONES,LM,2,JGRID_JL(n))
      n = JL_VEPAC
      SCALET=scale_jl(n)/idacc(ia_jl(n))
      CALL JLMAP(LNAME_JL(n),SNAME_JL(n),UNITS_JL(n),
     &     PLM,AJL(1,1,n),SCALET,ONES,ONES,LM,2,JGRID_JL(n))
      n = JL_WEPAC
      SCALET=scale_jl(n)/idacc(ia_jl(n))
      CALL JLMAP(LNAME_JL(n),SNAME_JL(n),UNITS_JL(n),
     &     PLE,AJL(1,1,n),SCALET,BYDXYP,ONES,LM-1,2,JGRID_JL(n))
C**** U, V AND W VELOCITY FOR WEST PACIFIC
      n = JL_UWPAC
      SCALET=scale_jl(n)/idacc(ia_jl(n))
      CALL JLMAP(LNAME_JL(n),SNAME_JL(n),UNITS_JL(n),
     &     PLM,AJL(1,1,n),SCALET,ONES,ONES,LM,2,JGRID_JL(n))
      n = JL_VWPAC
      SCALET=scale_jl(n)/idacc(ia_jl(n))
      CALL JLMAP(LNAME_JL(n),SNAME_JL(n),UNITS_JL(n),
     &     PLM,AJL(1,1,n),SCALET,ONES,ONES,LM,2,JGRID_JL(n))
      n = JL_WWPAC
      SCALET=scale_jl(n)/idacc(ia_jl(n))
      CALL JLMAP(LNAME_JL(n),SNAME_JL(n),UNITS_JL(n),
     &     PLE,AJL(1,1,n),SCALET,BYDXYP,ONES,LM-1,2,JGRID_JL(n))
C****
C**** ELIASSEN-PALM FLUX : NORTHWARD, VERTICAL, DIVERGENCE
C****
      n = JL_EPFLXN
      SCALET=scale_jl(n)/idacc(ia_jl(n))
      CALL JLMAP(LNAME_JL(n),SNAME_JL(n),UNITS_JL(n),
     &     PLM,AJL(1,1,n),SCALET,DXCOSV,ONES,LM,1,JGRID_JL(n))
      n = JL_EPFLXV
      SCALEV=scale_jl(n)/idacc(ia_jl(n))
      CALL JLMAP(LNAME_JL(n),SNAME_JL(n),UNITS_JL(n),
     &     PLE,AJL(1,1,n),SCALEV,COSP,ONES,LM-1,1,JGRID_JL(n))
      DXCVS=DXCOSV(2)
      DO J=2,JM-1
      BDN=0.
      DXCVN=DXCOSV(J+1)
      DO L=1,LM
      BUP=AJL(J,L,JL_EPFLXV)*COSP(J)
      AX(J,L)=AJL(J+1,L,JL_EPFLXN)*DXCVN-AJL(J,L,JL_EPFLXN)*DXCVS+
     *   .125*(BUP-BDN)/DSIG(L)
      BDN=BUP
      ENDDO
      DXCVS=DXCVN
      ENDDO
      DO 550 L=1,LM
      AX(1,L)=0.
  550 AX(JM,L)=0.
      n = jkl_epflx_div
      CALL JLMAP(LNAME_jkl(n),SNAME_jkl(n),UNITS_jkl(n),
     &     PLM,AX,SCALET,ONES,ONES,LM,1,1)

C****
C**** FOURIER ANALYSIS OF GEOPOTENTIAL HEIGHTS FOR WAVE NUMBERS 1 TO 4,
C****   AMPLITUDE AND PHASE
C****
            LINECT=63
      ELOFIM=.5*TWOPI-TWOPI/FIM

      DO K=1,kgz_max
      DO N=1,4
      AMPLTD(1,K,N)=0.
      AMPLTD(JM,K,N)=0.
      PHASE(1,K,N)=0.
      PHASE(JM,K,N)=0.
      ENDDO
      DO J=2,JM-1
      CALL FFT (AIJ(1,J,IJ_PHI1K-1+K),AN,BN)
      DO N=1,4
      AMPLTD(J,K,N)=SQRT(AN(N)*AN(N)+BN(N)*BN(N))
      PHASE(J,K,N)=(ATAN2(BN(N),AN(N))-TWOPI)/N+ELOFIM
      IF (PHASE(J,K,N).LE.-.5*TWOPI) PHASE(J,K,N)=PHASE(J,K,N)+TWOPI
      PHASE(J,K,N)=-PHASE(J,K,N)
      ENDDO
      ENDDO
      ENDDO
      SCALET=BYIADA*BYGRAV
      IX = jkl_phi_amp_wave1-1
      DO N=1,4
      CALL JLMAP(LNAME_jkl(N+ix),SNAME_jkl(N+ix),UNITS_jkl(N+ix),
     &      PMB,AMPLTD(1,1,N),SCALET,ONES,ONES,kgz_max,2,1)
      ENDDO
      SCALET=360./TWOPI
      IX = jkl_phi_phase_wave1-1
      DO N=1,4
      CALL JLMAP(LNAME_jkl(N+ix),SNAME_jkl(N+ix),UNITS_jkl(N+ix),
     &      PMB,PHASE(1,1,N),SCALET,ONES,ONES,kgz_max,2,1)
      ENDDO

      if(qcheck) call close_jl

      RETURN
  901 FORMAT (
     *  '010**14 WATTS = .2067 * 10**19 CALORIES/DAY'/
     *  ' 10**18 JOULES = .864 * 10**30 GM*CM**2/SEC/DAY')
      END SUBROUTINE DIAGJK


      SUBROUTINE JKMAP(LNAME,SNAME,UNITS,
     &     PM,AX,SCALET,SCALEJ,SCALEK,KMAX,JWT,J1)
      USE DAGCOM, only : QCHECK,acc_period,iu_jl,lm_req
      USE DAGPCOM, only : inc=>incj,linect,jmby2
      USE MODEL_COM, only :
     &     jm,lm,JDATE,JDATE0,JMON0,JMON,AMON0,AMON,JYEAR,JYEAR0,XLABEL
      USE GEOM, only :
     &     LAT_DG,WTJ
      IMPLICIT NONE

!@var units string containing output field units
      CHARACTER(LEN=50) :: UNITS
!@var lname string describing output field
      CHARACTER(LEN=50) :: LNAME
!@var sname string referencing output field
      CHARACTER(LEN=30) :: SNAME
!@var title string, formed as concatentation of lname//units
      CHARACTER(LEN=64) :: TITLE

      INTEGER, DIMENSION(JM) :: MLAT
      DOUBLE PRECISION, DIMENSION(JM) :: FLAT,ASUM
      DOUBLE PRECISION, DIMENSION(2) :: AHEM

      DOUBLE PRECISION, DIMENSION(JM,LM,2) :: DSJK
      DOUBLE PRECISION, DIMENSION(2,LM,2) :: DSHEM
      DOUBLE PRECISION, DIMENSION(LM,2) :: DSGLOB
      COMMON/WORK5/DSJK,DSHEM,DSGLOB


      INTEGER :: J1,JWT,KMAX
      DOUBLE PRECISION :: SCALET,SCALER
      DOUBLE PRECISION, DIMENSION(JM,LM) :: AX
      DOUBLE PRECISION, DIMENSION(JM,LM_REQ) :: ARQX
      DOUBLE PRECISION, DIMENSION(JM) :: SCALEJ,SCALJR
      DOUBLE PRECISION, DIMENSION(LM) :: SCALEK
      DOUBLE PRECISION, DIMENSION(LM_REQ) :: SCALLR
      DOUBLE PRECISION, DIMENSION(LM+LM_REQ) :: PM

      DOUBLE PRECISION, DIMENSION(JM,LM) :: CX

      CHARACTER*4 DASH,WORD(4)
      DATA DASH/'----'/,WORD/'SUM','MEAN',' ','.1*'/

      INTEGER :: IWORD,J,J0,JH,JHEMI,K,L ,ksx,klmax
      DOUBLE PRECISION :: AGLOB,FGLOB,FLATJ,G1,H1,H2,SUMFAC

      REAL*8, DIMENSION(JM+3,LM+LM_REQ+1) :: XJL ! for binary output
      CHARACTER XLB*16,CLAT*16,CPRES*16,CBLANK*16,TITLEO*80
      DATA CLAT/'LATITUDE'/,CPRES/'PRESSURE (MB)'/,CBLANK/' '/

C form title string
      title = trim(lname)//' ('//trim(units)//')'
C****
C**** PRODUCE A LATITUDE BY LAYER TABLE OF THE ARRAY A
C****
   10 LINECT=LINECT+KMAX+7
      IF (LINECT.LE.60) GO TO 20
      WRITE (6,907) XLABEL(1:105),JDATE0,AMON0,JYEAR0,JDATE,AMON,JYEAR
      LINECT=KMAX+8
   20 WRITE (6,901) TITLE,(DASH,J=J1,JM,INC)
      WRITE (6,904) WORD(JWT),(NINT(LAT_DG(J,J1)),J=JM,J1,-INC)
      WRITE (6,905) (DASH,J=J1,JM,INC)
      J0=J1-1
         DO 40 L=1,LM+LM_REQ+1
         DO 40 J=1,JM+3
   40    XJL(J,L) = -1.E30
         KSX = 0            ! KSX = LAYERS GENERATED AT ENTRY
  100 DO 110 J=J1,JM
      DO 110 K=1,KMAX
  110 CX(J,K)=AX(J,K)*SCALET*SCALEJ(J)*SCALEK(K)
         KLMAX = KMAX+KSX
C**** HORIZONTAL SUMS AND TABLE ENTRIES
      DO 140 K=KMAX,1,-1
      AGLOB=0.
      DO 130 JHEMI=1,2
      AHEM(JHEMI)=0.
      DO 120 JH=1,jmby2
      J=(JHEMI-1)*(jmby2-J0)+JH+J0
      FLAT(J)=CX(J,K)/(DSJK(J,K,J1)+1.D-20)
         XJL(J,K) = FLAT(J)
         IF (DSJK(J,K,J1).EQ.0.) XJL(J,K) = -1.E30
      MLAT(J)=NINT(FLAT(J))
  120 AHEM(JHEMI)=AHEM(JHEMI)+CX(J,K)*WTJ(J,JWT,J1)
  130 AGLOB=AGLOB+AHEM(JHEMI)/JWT
      H1=AHEM(1)/(DSHEM(1,K,J1)+1.D-20)
      H2=AHEM(2)/(DSHEM(2,K,J1)+1.D-20)
      G1=AGLOB/(DSGLOB(K,J1)+1.D-20)
         XJL(JM+3,K)=H1   ! SOUTHERN HEM
         XJL(JM+2,K)=H2   ! NORTHERN HEM
         XJL(JM+1,K)=G1   ! GLOBAL
      WRITE (6,902) PM(K),G1,H2,H1,(MLAT(J),J=JM,J1,-INC)
         CALL KEYNRL (SNAME,K,FLAT)
  140 CONTINUE
C**** VERTICAL SUMS
      WRITE (6,905) (DASH,J=J1,JM,INC)
      SUMFAC=1.
      IWORD=3
      IF ( SNAME.EQ.'temp' .OR. ! make sumfac an argument to avoid this
     &     SNAME.EQ.'v' .OR.
     &     SNAME.EQ.'tot_nt_dse' .OR.
     &     SNAME.EQ.'tot_nt_lh' .OR.
     &     SNAME.EQ.'tot_nt_se' .OR.
     &     SNAME.EQ.'tot_nt_am') THEN
         SUMFAC=10.
         IWORD=4
      ENDIF
      DO 180 J=J1,JM
      ASUM(J)=0.
      DO 170 K=1,KMAX
  170 ASUM(J)=ASUM(J)+CX(J,K)
         XJL(J,LM+LM_REQ+1)=ASUM(J)
  180 MLAT(J)=NINT(ASUM(J)*SUMFAC)
      AGLOB=0.
      DO 200 JHEMI=1,2
      AHEM(JHEMI)=0.
      DO 190 JH=1,jmby2
      J=(JHEMI-1)*(jmby2-J0)+JH+J0
      AHEM(JHEMI)=AHEM(JHEMI)+ASUM(J)*WTJ(J,JWT,J1)*SUMFAC
  190 CONTINUE
  200 AGLOB=AGLOB+AHEM(JHEMI)/JWT
         XJL(JM+3,LM+LM_REQ+1)=AHEM(1)   ! SOUTHERN HEM
         XJL(JM+2,LM+LM_REQ+1)=AHEM(2)   ! NORTHERN HEM
         XJL(JM+1,LM+LM_REQ+1)=AGLOB     ! GLOBAL
         XLB=' '//acc_period(1:3)//' '//acc_period(4:12)//'  '
         TITLEO=TITLE//XLB
         IF(QCHECK) CALL POUT_JL(TITLEO,LNAME,SNAME,UNITS,
     *        J1,KLMAX,XJL,PM,CLAT,CPRES)
      WRITE (6,903) WORD(IWORD),AGLOB,AHEM(2),AHEM(1),
     *  (MLAT(J),J=JM,J1,-INC)
         CALL KEYVSUMS(SNAME,AGLOB,AHEM,ASUM,SUMFAC)
      RETURN
C****
      ENTRY JKMAPS(LNAME,SNAME,UNITS,
     &     PM,AX,SCALET,SCALEJ,SCALEK,KMAX,JWT,J1,
     *  ARQX,SCALER,SCALJR,SCALLR)
C form title string
      title = trim(lname)//' ('//trim(units)//')'
         KSX = 3
         DO 205 L=1,LM+LM_REQ+1
         DO 205 J=1,JM+3
  205    XJL(J,L) = -1.E30
      LINECT=LINECT+KMAX+10
      IF (LINECT.LE.60) GO TO 230
      WRITE (6,907) XLABEL(1:105),JDATE0,AMON0,JYEAR0,JDATE,AMON,JYEAR
      LINECT=KMAX+11
  230 J0=J1-1
C**** PRODUCE UPPER STRATOSPHERE NUMBERS FIRST
      WRITE (6,901) TITLE,(DASH,J=J1,JM,INC)
      WRITE (6,904) WORD(JWT),(NINT(LAT_DG(J,J1)),J=JM,J1,-INC)
      WRITE (6,905) (DASH,J=J1,JM,INC)
      DO 260 L=LM_REQ,1,-1
      FGLOB=0.
      DO 250 JHEMI=1,2
      AHEM(JHEMI)=0.
      DO 240 JH=1,jmby2
      J=(JHEMI-1)*(jmby2-J0)+JH-J0
      FLATJ=ARQX(J,L)*SCALER*SCALJR(J)*SCALLR(L)
         XJL(J,L+KMAX) = FLATJ
      MLAT(J)=NINT(FLATJ)
  240 AHEM(JHEMI)=AHEM(JHEMI)+FLATJ*WTJ(J,JWT,J1)
  250 FGLOB=FGLOB+AHEM(JHEMI)/JWT
         XJL(JM+3,L+KMAX)=AHEM(1)   ! SOUTHERN HEM
         XJL(JM+2,L+KMAX)=AHEM(2)   ! NORTHERN HEM
         XJL(JM+1,L+KMAX)=FGLOB     ! GLOBAL
  260 WRITE (6,902) PM(L+LM),FGLOB,AHEM(2),AHEM(1),
     *  (MLAT(J),J=JM,J1,-INC)
      GO TO 100
  901 FORMAT ('0',30X,A64,'  CP'/1X,32('-'),24A4)
  902 FORMAT (1X,F7.3,3F8.1,1X,24I4)
  903 FORMAT (A6,2X,3F8.1,1X,24I4)
  904 FORMAT (' P(MB)   ',A4,' G      NH      SH  ',24I4)
  905 FORMAT (1X,32('-'),24A4)
  907 FORMAT ('1',A,I3,1X,A3,I5,' - ',I3,1X,A3,I5)
      END SUBROUTINE JKMAP


      SUBROUTINE JLMAP(LNAME,SNAME,UNITS,
     &     PL,AX,SCALET,SCALEJ,SCALEL,LMAX,JWT,J1)
C****
C**** THIS ROUTINE PRODUCES LAYER BY LATITUDE TABLES ON THE LINE
C**** PRINTER.  THE INTERIOR NUMBERS OF THE TABLE ARE CALCULATED AS
C****               AX * SCALET * SCALEJ * SCALEL.
C**** WHEN JWT=1, THE INSIDE NUMBERS ARE NOT AREA WEIGHTED AND THE
C****    HEMISPHERIC AND GLOBAL NUMBERS ARE SUMMATIONS.
C**** WHEN JWT=2, ALL NUMBERS ARE PER UNIT AREA.
C**** J1 INDICATES PRIMARY OR SECONDARY GRID.
C**** THE BOTTOM LINE IS CALCULATED AS THE SUMMATION OF DSIG TIMES THE
C**** NUMBERS ABOVE (POSSIBLY MULTIPLIED BY A FACTOR OF 10)
C****
      USE DAGCOM, only : QCHECK,acc_period,iu_jl,LM_REQ
      USE DAGPCOM, only : inc=>incj,linect,jmby2
      USE MODEL_COM, only :
     &     jm,lm,DSIG,JDATE,JDATE0,AMON,AMON0,JYEAR,JYEAR0,SIGE,XLABEL
      USE GEOM, only :
     &     LAT_DG,WTJ
      IMPLICIT NONE

!@var units string containing output field units
      CHARACTER(LEN=50) :: UNITS
!@var lname string describing output field
      CHARACTER(LEN=50) :: LNAME
!@var sname string referencing output field
      CHARACTER(LEN=30) :: SNAME
!@var title string, formed as concatentation of lname//units
      CHARACTER(LEN=64) :: TITLE

      INTEGER, DIMENSION(JM) :: MLAT
      DOUBLE PRECISION, DIMENSION(JM) :: FLAT,ASUM
      DOUBLE PRECISION, DIMENSION(2) :: FHEM,HSUM


      INTEGER :: J1,JWT,LMAX
      DOUBLE PRECISION :: SCALET,SCALER
      DOUBLE PRECISION, DIMENSION(JM,LM) :: AX
      DOUBLE PRECISION, DIMENSION(JM,LM_REQ) :: ARQX
      DOUBLE PRECISION, DIMENSION(JM) :: SCALEJ,SCALJR
      DOUBLE PRECISION, DIMENSION(LM) :: SCALEL
      DOUBLE PRECISION, DIMENSION(LM_REQ) :: SCALLR
      DOUBLE PRECISION, DIMENSION(LM+LM_REQ) :: PL

      CHARACTER*4 DASH,WORD(4)
      DATA DASH/'----'/,WORD/'SUM','MEAN',' ','.1*'/

      INTEGER :: IWORD,J,J0,JH,JHEMI,K,L  ,ksx,klmax
      DOUBLE PRECISION :: FGLOB,FLATJ,GSUM,SDSIG,SUMFAC

      REAL*8, DIMENSION(JM+3,LM+LM_REQ+1) :: XJL ! for binary output
      CHARACTER XLB*16,CLAT*16,CPRES*16,CBLANK*16,TITLEO*80
      DATA CLAT/'LATITUDE'/,CPRES/'PRESSURE (MB)'/,CBLANK/' '/

C form title string
      title = trim(lname)//' ('//trim(units)//')'
C****
C**** PRODUCE A LATITUDE BY LAYER TABLE OF THE ARRAY A
C****
   10 LINECT=LINECT+LMAX+7
      IF (LINECT.LE.60) GO TO 20
      WRITE (6,907) XLABEL(1:105),JDATE0,AMON0,JYEAR0,JDATE,AMON,JYEAR
      LINECT=LMAX+8
   20 WRITE (6,901) TITLE,(DASH,J=J1,JM,INC)
      WRITE (6,904) WORD(JWT),(NINT(LAT_DG(J,J1)),J=JM,J1,-INC)
      WRITE (6,905) (DASH,J=J1,JM,INC)
      J0=J1-1
         DO 40 L=1,LM+LM_REQ+1
         DO 40 J=1,JM+3
   40    XJL(J,L) = -1.E30
         KSX = 0            ! KSX = LAYERS GENERATED AT ENTRY
  100 SDSIG=1.-SIGE(LMAX+1)
         KLMAX = LMAX+KSX
      DO 110 J=1,JM
  110 ASUM(J)=0.
      HSUM(1)=0.
      HSUM(2)=0.
      GSUM=0.
      SUMFAC=1.
      IWORD=3
      if(sname.eq.'dudt_mtndrg') then ! make sumfac an argument ... ???
         SUMFAC=10.                   ! ... to avoid this if-block  ???
         IWORD=4
      endif
      DO 140 L=LMAX,1,-1
      FGLOB=0.
      DO 130 JHEMI=1,2
      FHEM(JHEMI)=0.
      DO 120 JH=1,jmby2
      J=(JHEMI-1)*(jmby2-J0)+JH+J0
      FLAT(J)=AX(J,L)*SCALET*SCALEJ(J)*SCALEL(L)
         XJL(J,L) = FLAT(J)
      MLAT(J)=NINT(FLAT(J))
  115 ASUM(J)=ASUM(J)+FLAT(J)*DSIG(L)/SDSIG
  120 FHEM(JHEMI)=FHEM(JHEMI)+FLAT(J)*WTJ(J,JWT,J1)
  130 FGLOB=FGLOB+FHEM(JHEMI)/JWT
         XJL(JM+3,L)=FHEM(1)   ! SOUTHERN HEM
         XJL(JM+2,L)=FHEM(2)   ! NORTHERN HEM
         XJL(JM+1,L)=FGLOB     ! GLOBAL
      WRITE (6,902) PL(L),FGLOB,FHEM(2),FHEM(1),(MLAT(J),J=JM,J1,-INC)
         CALL KEYNRL (SNAME,L,FLAT)
  136 HSUM(1)=HSUM(1)+FHEM(1)*SUMFAC*DSIG(L)/SDSIG
      HSUM(2)=HSUM(2)+FHEM(2)*SUMFAC*DSIG(L)/SDSIG
  140 GSUM=GSUM+FGLOB*SUMFAC*DSIG(L)/SDSIG
      WRITE (6,905) (DASH,J=J1,JM,INC)
      ASUM(jmby2+1)=ASUM(jmby2+1)/J1
      DO 150 J=J1,JM
  150 MLAT(J)=NINT(ASUM(J)*SUMFAC)
         DO 180 J=J1,JM
  180    XJL(J   ,LM+LM_REQ+1)=ASUM(J)
         XJL(JM+3,LM+LM_REQ+1)=HSUM(1)   ! SOUTHERN HEM
         XJL(JM+2,LM+LM_REQ+1)=HSUM(2)   ! NORTHERN HEM
         XJL(JM+1,LM+LM_REQ+1)=GSUM      ! GLOBAL
         XLB=' '//acc_period(1:3)//' '//acc_period(4:12)//'  '
         TITLEO=TITLE//XLB
         IF(QCHECK) CALL POUT_JL(TITLEO,LNAME,SNAME,UNITS,
     *        J1,KLMAX,XJL,PL,CLAT,CPRES)
      if(  sname.eq.'phi_amp_wave1' .or.
     &     sname.eq.'phi_amp_wave2' .or.
     &     sname.eq.'phi_amp_wave3' .or.
     &     sname.eq.'phi_amp_wave4' .or.
     &     sname.eq.'phi_phase_wave1' .or.
     &     sname.eq.'phi_phase_wave2' .or.
     &     sname.eq.'phi_phase_wave3' .or.
     &     sname.eq.'phi_phase_wave4' ) return
      WRITE (6,903) WORD(IWORD),GSUM,HSUM(2),HSUM(1),
     *  (MLAT(J),J=JM,J1,-INC)
      RETURN
C****
      ENTRY JLMAPS(LNAME,SNAME,UNITS,
     &     PL,AX,SCALET,SCALEJ,SCALEL,LMAX,JWT,J1,
     *  ARQX,SCALER,SCALJR,SCALLR)

C form title string
      title = trim(lname)//' ('//trim(units)//')'

         KSX = 3
         DO 205 L=1,LM+LM_REQ
         DO 205 J=1,JM
  205    XJL(J,L) = -1.E30
      LINECT=LINECT+LMAX+10
      IF (LINECT.LE.60) GO TO 200
      WRITE (6,907) XLABEL(1:105),JDATE0,AMON0,JYEAR0,JDATE,AMON,JYEAR
      LINECT=LMAX+11
  200 J0=J1-1
C**** PRODUCE UPPER STRATOSPHERE NUMBERS FIRST
      WRITE (6,901) TITLE,(DASH,J=J1,JM,INC)
      WRITE (6,904) WORD(JWT),(NINT(LAT_DG(J,J1)),J=JM,J1,-INC)
      WRITE (6,905) (DASH,J=J1,JM,INC)
      DO 230 L=LM_REQ,1,-1
      FGLOB=0.
      DO 220 JHEMI=1,2
      FHEM(JHEMI)=0.
      DO 210 JH=1,jmby2
      J=(JHEMI-1)*(jmby2-J0)+JH-J0
      FLATJ=ARQX(J,L)*SCALER*SCALJR(J)*SCALLR(L)
         XJL(J,L+LMAX) = FLATJ
      MLAT(J)=NINT(FLATJ)
  210 FHEM(JHEMI)=FHEM(JHEMI)+FLATJ*WTJ(J,JWT,J1)
  220 FGLOB=FGLOB+FHEM(JHEMI)/JWT
         XJL(JM+3,L+LMAX)=FHEM(1)   ! SOUTHERN HEM
         XJL(JM+2,L+LMAX)=FHEM(2)   ! NORTHERN HEM
         XJL(JM+1,L+LMAX)=FGLOB     ! GLOBAL
  230 WRITE (6,902) PL(L+LM),FGLOB,FHEM(2),FHEM(1),
     *  (MLAT(J),J=JM,J1,-INC)
      GO TO 100
  901 FORMAT ('0',30X,A64/2X,32('-'),24A4)
  902 FORMAT (1X,F8.3,3F8.1,1X,24I4)
  903 FORMAT (1X,A6,2X,3F8.1,1X,24I4)
  904 FORMAT ('  P(MB)   ',A4,' G      NH      SH  ',24I4)
  905 FORMAT (2X,32('-'),24A4)
  907 FORMAT ('1',A,I3,1X,A3,I5,' - ',I3,1X,A3,I5)
      END SUBROUTINE JLMAP


      SUBROUTINE DIAGIL
!@sum  DIAGIL prints out longitude/height diagnostics
!@auth Original Development Team
!@ver  1.0
      USE MODEL_COM, only : im,lm,bydsig,idacc,xlabel,lrunid
      USE DAGPCOM, only : plm,ple,linect
      USE DAGCOM, only : ail,lm_req,acc_period, qcheck,lname_il,name_il
     *     ,units_il,scale_il,ia_il,kail
      IMPLICIT NONE
      CHARACTER sname*20,unit*20,lname*80
      DOUBLE PRECISION, DIMENSION(LM) :: ONES
      DOUBLE PRECISION, DIMENSION(IM,LM) :: XIL
      INTEGER :: K

C**** OPEN PLOTTABLE OUTPUT FILE IF DESIRED
      IF(QCHECK) call open_il(trim(acc_period)//'.il'//XLABEL(1:LRUNID))

C**** INITIALIZE CERTAIN QUANTITIES
      ONES(1:LM)=1.

      linect = 65

      DO K=1,KAIL
        XIL=AIL(:,:,K)*SCALE_IL(K)/IDACC(IA_IL(K))
        sname=name_il(k)
        lname=lname_il(k)
        unit=units_il(k)
        SELECT CASE (sname)
! Centered in L; secondary grid; hor. mean; vert. sum
        CASE ('u_equator','v_equator','u_70N','u_50N')
          CALL ILMAP(sname,lname,unit,PLM,XIL,ONES,LM,2,2)
! Vertical edges; primary grid; hor. mean; vert. sum
        CASE ('vvel_equator','vvel_50N','vvel_70N')
          CALL ILMAP(sname,lname,unit,PLE,XIL,ONES,LM-1,2,1)
! Centered in L; primary grid; hor. mean; vert. sum
        CASE ('temp_equator','rh_equator','temp_50N','temp_70N')
          CALL ILMAP(sname,lname,unit,PLM,XIL,ONES,LM,2,1)
! Centered in L; primary grid; hor. sum; vert. sum
        CASE ('mcheat_equator')
          CALL ILMAP(sname,lname,unit,PLM,XIL,ONES,LM,1,1)
! Centered in L; primary grid; hor. sum; vert. mean
        CASE ('rad_cool_equator') ! also 'rad_cool_50N','rad_cool_70N'
          CALL ILMAP(sname,lname,unit,PLM,XIL,BYDSIG,LM,1,1)
        END SELECT
      END DO
      if(qcheck) call close_il
      RETURN
      END SUBROUTINE DIAGIL


      SUBROUTINE ILMAP (sname,lname,unit,PL,AX,SCALEL,LMAX,JWT
     *     ,ISHIFT)
      USE DAGCOM, only : qcheck,acc_period,iu_il
      USE DAGPCOM, only : inc=>inci,linect
      USE CONSTANT, only : twopi
      USE MODEL_COM, only : im,jm,lm,dsig,jdate,jdate0,amon,amon0,jyear
     *     ,jyear0,sige,xlabel
      USE GEOM, only : dlon,lon_dg
      IMPLICIT NONE
      CHARACTER XLB*80,CWORD*8
      character(len=20), intent(in) :: sname,unit
      character(len=80), intent(in) :: lname
      CHARACTER*64 :: TITLE
      CHARACTER*4, PARAMETER :: DASH='----'
      CHARACTER*4, DIMENSION(2), PARAMETER :: WORD=(/'SUM ','MEAN'/)
      CHARACTER*16, PARAMETER :: CBLANK=' ', CLAT='LONGITUDE',
     *     CPRES='PRESSURE (MB)'
      DOUBLE PRECISION, DIMENSION(LM), INTENT(IN) :: PL,SCALEL
      DOUBLE PRECISION, DIMENSION(IM,LM), INTENT(IN) :: AX
      REAL*8 :: XIL(IM,LM),ZONAL(LM) ! used for post-proc
      DOUBLE PRECISION :: FGLOB,FLON,GSUM,SDSIG
      DOUBLE PRECISION, DIMENSION(IM) :: ASUM
      INTEGER, DIMENSION(IM) :: MLON
      INTEGER, INTENT(IN) :: JWT,ISHIFT
      INTEGER :: I,K,L,LMAX
C****
C**** PRODUCE A LONGITUDE BY LAYER TABLE OF THE ARRAY A
C****
!@var ISHIFT: When=2, print longitude indices off center (U-grid)
      LINECT=LINECT+LMAX+7
      IF (LINECT.GT.60) THEN
        WRITE (6,907) XLABEL(1:105),JDATE0,AMON0,JYEAR0,JDATE,AMON,JYEAR
        LINECT=LMAX+8
      END IF
      SDSIG=1.-SIGE(LMAX+1)

      TITLE=trim(lname)//" ("//trim(unit)//")"
      WRITE (6,901) TITLE,(DASH,I=1,IM,INC)
      IF (ISHIFT.EQ.1) WRITE (6,904) WORD(JWT),(I,I=1,IM,INC)
      IF (ISHIFT.EQ.2) WRITE (6,906) WORD(JWT),(I,I=1,IM,INC)
      WRITE (6,905) (DASH,I=1,IM,INC)

      ASUM(:)=0. ; GSUM=0.
      DO L=LMAX,1,-1
        FGLOB=0.
        DO I=1,IM
          FLON=AX(I,L)*SCALEL(L)
          XIL(I,L)=FLON
          MLON(I)=NINT(FLON)
          ASUM(I)=ASUM(I)+FLON*DSIG(L)/SDSIG
          FGLOB=FGLOB+FLON
        END DO
        FGLOB=FGLOB/IM
        IF (JWT.EQ.1) FGLOB=FGLOB*TWOPI/DLON
        ZONAL(L)=FGLOB
        WRITE (6,902) PL(L),FGLOB,(MLON(I),I=1,IM,INC)
        GSUM=GSUM+FGLOB*DSIG(L)/SDSIG
      END DO
      MLON(:)=NINT(ASUM(:))

      WRITE (6,905) (DASH,I=1,IM,INC)
      WRITE (6,903) GSUM,(MLON(I),I=1,IM,INC)
C**** Output for post-processing
      CWORD=WORD(JWT)           ! pads out to 8 characters
      XLB=TITLE
      XLB(65:80)=' '//acc_period(1:3)//' '//acc_period(4:12)//'  '
      IF(QCHECK) CALL POUT_IL(XLB,sname,lname,unit,ISHIFT,LMAX,XIL
     *     ,PL,CLAT,CPRES,ASUM,GSUM,ZONAL)
      RETURN
C****
  901 FORMAT ('0',30X,A64,2('-'),/1X,16('-'),36A3)
  902 FORMAT (F8.3,F8.1,1X,36I3)
  903 FORMAT (F16.1,1X,36I3)
  904 FORMAT (' P(MB)',6X,A4,1X,36I3)  ! U-grid (i.e., centers)
  905 FORMAT (1X,16('-'),36A3)
  906 FORMAT (' P(MB)',6X,A4,36I3)     ! V-grid (i.e., edges)
  907 FORMAT ('1',A,I3,1X,A3,I5,' - ',I3,1X,A3,I5)
C****
      END SUBROUTINE ILMAP


      BLOCK DATA BDWP
C****
C**** TITLES FOR SUBROUTINE DIAG7
C****
      COMMON/D7COM/TITLE
      CHARACTER*64 TITLE(12)
      DATA TITLE/
     1'WAVE POWER FOR U NEAR 850 MB AND EQUATOR (DAY*(M/S)**2)       ',
     2'WAVE POWER FOR V NEAR 850 MB AND EQUATOR (DAY*(M/S)**2)       ',
     3'WAVE POWER FOR U NEAR 300 MB AND EQUATOR (10 DAY*(M/S)**2)    ',
     4'WAVE POWER FOR V NEAR 300 MB AND EQUATOR (DAY*(M/S)**2)       ',
     5'WAVE POWER FOR U NEAR 50 MB AND EQUATOR (10 DAY*(M/S)**2)     ',
     6'WAVE POWER FOR V NEAR 50 MB AND EQUATOR (DAY*(M/S)**2)        ',
     7'WAVE POWER FOR PHI AT 922 MB AND 50 DEG NORTH (10**3 DAY*M**2)',
     8'WAVE POWER FOR PHI AT 700 MB AND 50 DEG NORTH (10**3 DAY*M**2)',
     9'WAVE POWER FOR PHI AT 500 MB AND 50 DEG NORTH (10**3 DAY*M**2)',
     A'WAVE POWER FOR PHI AT 300 MB AND 50 DEG NORTH (10**3 DAY*M**2)',
     B'WAVE POWER FOR PHI AT 100 MB AND 50 DEG NORTH (10**4 DAY*M**2)',
     C'WAVE POWER FOR PHI AT 10 MB AND 50 DEG NORTH (10**4 DAY*M**2) '/
      END BLOCK DATA BDWP


      SUBROUTINE DIAG7P
C****
C**** THIS ENTRY PRINTS THE TABLES
C****
c      USE PRTCOM, only :
      USE DAGCOM, only :
     &     nwav_dag,wave,Max12HR_sequ,Min12HR_sequ
      USE MODEL_COM, only :
     &     im,IDACC,JDATE,JDATE0,AMON,AMON0,JYEAR,JYEAR0,XLABEL
      IMPLICIT NONE

      INTEGER, DIMENSION(44) :: IPOWER
      DOUBLE PRECISION, DIMENSION(120) :: POWER
      DOUBLE PRECISION, DIMENSION(43) :: XPOWER
      DOUBLE PRECISION, DIMENSION(13) :: FPE
      DOUBLE PRECISION, DIMENSION(13,3,4) :: FFPE
      DOUBLE PRECISION, DIMENSION(43,9,3,4) :: FPOWER

      INTEGER, PARAMETER :: MMAX=12,NUAMAX=120,NUBMAX=15

      COMMON/D7COM/TITLE
      CHARACTER*64 TITLE(12)

      DOUBLE PRECISION, DIMENSION(12) :: SCALET
      DATA SCALET/1.,1., .1,1., .1,1., 4*1.D-3,1.D-4,1.D-5/

      DOUBLE PRECISION ::
     &     BYIA12,PNU,POWX,VAR

      INTEGER ::
     &     IDACC9,K,KPAGE,KQ,KTABLE,
     &     M,MMAXP1,N,NMAX,NS,NUA,NX

      NMAX=NWAV_DAG
      IDACC9=IDACC(9)
      IF (IDACC9.LE.MMAX) RETURN
C**** PATCH NEEDED IF SEVERAL RESTART FILES WERE ACCUMULATED
      IF (IDACC(12).LE.1) GO TO 320
      IDACC9=Min12HR_sequ           ! in case a February was included
      BYIA12=1./IDACC(12)
      WAVE(:,:,:,:)=WAVE(:,:,:,:)*BYIA12
  320 CONTINUE
      IF (IDACC9.GT.Max12HR_sequ) IDACC9=Max12HR_sequ
C****
C**** OUTPUT WAVE POWER AT THE EQUATOR
C****
      MMAXP1=MMAX+1
      DO 400 KPAGE=1,2
      WRITE (6,907) XLABEL(1:105),JDATE0,AMON0,JYEAR0,JDATE,AMON,JYEAR
      DO 390 KTABLE=1,3
      KQ=3*(KPAGE-1)+KTABLE
      WRITE (6,901) TITLE(KQ)
      DO 380 NX=1,NMAX
      N=NMAX+1-NX
      CALL MEM (WAVE(1,1,N,KQ),IDACC9,MMAX,NUAMAX,NUBMAX,POWER,FPE,
     *  VAR,PNU)
      POWX=.5*POWER(1)
      DO 330 NUA=2,27
  330 POWX=POWX+POWER(NUA)
      XPOWER(1)=SCALET(KQ)*POWX/26.5
      POWX=0.
      DO 340 NUA=28,34
  340 POWX=POWX+POWER(NUA)
      XPOWER(2)=SCALET(KQ)*POWX/7.
      XPOWER(3)=SCALET(KQ)*(POWER(35)+POWER(36)+POWER(37)+POWER(38))/4.
      XPOWER(4)=SCALET(KQ)*(POWER(39)+POWER(40))/2.
      DO 350 NUA=41,76
  350 XPOWER(NUA-36)=SCALET(KQ)*POWER(NUA)
      POWX=.5*POWER(1)
      DO 360 NUA=77,120
  360 POWX=POWX+POWER(NUA)
      XPOWER(41)=SCALET(KQ)*POWX/44.5
      XPOWER(42)=10.*SCALET(KQ)*VAR
      XPOWER(43)=1000.*SCALET(KQ)*(VAR-PNU)
      DO 370 NS=1,43
         FPOWER(NS,NX,KTABLE,KPAGE)=XPOWER(NS)
      IPOWER(NS)=XPOWER(NS)+.5
  370 CONTINUE
  380 WRITE (6,902) N,(IPOWER(NS),NS=1,43)
         DO 385 M=1,MMAXP1
  385    FFPE(M,KTABLE,KPAGE)=FPE(M)
  390 WRITE (6,903) (FPE(M),M=1,MMAXP1)
  400 CONTINUE
C****
C**** OUTPUT WAVE POWER AT 50 DEG NORTH
C****
      DO 500 KPAGE=3,4
      WRITE (6,907) XLABEL(1:105),JDATE0,AMON0,JYEAR0,JDATE,AMON,JYEAR
      DO 490 KTABLE=1,3
      KQ=3*(KPAGE-1)+KTABLE
  410 WRITE (6,911) TITLE(KQ)
      DO 480 NX=1,NMAX
      N=NMAX+1-NX
      CALL MEM (WAVE(1,1,N,KQ),IDACC9,MMAX,NUAMAX,NUBMAX,POWER,FPE,
     *  VAR,PNU)
      DO 420 M=1,MMAXP1
  420 FPE(M)=1000.*SCALET(KQ)*FPE(M)
      POWX=.5*POWER(1)
      DO 430 NUA=2,45
  430 POWX=POWX+POWER(NUA)
      XPOWER(1)=SCALET(KQ)*POWX/44.5
      DO 440 NUA=46,81
  440 XPOWER(NUA-44)=SCALET(KQ)*POWER(NUA)
      XPOWER(38)=SCALET(KQ)*(POWER(82)+POWER(83))/2.
      XPOWER(39)=SCALET(KQ)*(POWER(84)+POWER(85)+POWER(86)+POWER(87))/4.
      POWX=0.
      DO 450 NUA=88,94
  450 POWX=POWX+POWER(NUA)
      XPOWER(40)=SCALET(KQ)*POWX/7.
      POWX=.5*POWER(1)
      DO 460 NUA=95,120
  460 POWX=POWX+POWER(NUA)
      XPOWER(41)=SCALET(KQ)*POWX/26.5
      XPOWER(42)=10.*SCALET(KQ)*VAR
      XPOWER(43)=1000.*SCALET(KQ)*(VAR-PNU)
      DO 470 NS=1,43
         FPOWER(NS,NX,KTABLE,KPAGE)=XPOWER(NS)
      IPOWER(NS)=XPOWER(NS)+.5
  470 CONTINUE
  480 WRITE (6,902) N,(IPOWER(NS),NS=1,43)
         DO 485 M=1,MMAXP1
  485    FFPE(M,KTABLE,KPAGE)=FPE(M)
  490 WRITE (6,903) (FPE(M),M=1,MMAXP1)
  500 CONTINUE
      RETURN
C****
  901 FORMAT ('0',30X,A64,8X,'*1/60 (1/DAY)'/'   PERIOD EASTWARD--',
     *   35('---')/' N    -2      *-3   -3.3      -4       -5    -6   -7
     *.5  -10-12-15-20-30-60    60 30 20 15 12 10    7.5    6     5
     *   4*   VAR ERR'/'   --',40('---'))
  902 FORMAT (I2,41I3,I4,I4)
  903 FORMAT ('   --',40('---')/(1X,13F10.4))
  907 FORMAT ('1',A,I3,1X,A3,I5,' - ',I3,1X,A3,I5)
  911 FORMAT ('0',30X,A64,8X,'*1/60 (1/DAY)'/'   PERIOD EASTWARD--',
     *  35('---')/               ' N   *-4       -5    -6   -7.5  -10-12
     *-15-20-30-60    60 30 20 15 12 10    7.5    6     5        4
     * 3.3    3*       2    VAR ERR'/'   --',40('---'))
      END SUBROUTINE DIAG7P


      SUBROUTINE MEM (SERIES,ITM,MMAX,NUAMAX,NUBMAX,POWER,FPE,VAR,PNU)
      USE CONSTANT, only : pi
      IMPLICIT NONE
      DIMENSION C(1800),S(1800),B1(62),B2(62),A(12),AA(11),P(13)
      DIMENSION SERIES(*),POWER(*),FPE(*)
      DOUBLE PRECISION ARG,PP,POWERX,P,C,S,POWER,FPE
      COMPLEX*16 CI,CSUM,SS,A,AA,B1,B2,ANOM,ADEN
      COMPLEX*16 SERIES
      DOUBLE PRECISION :: PNU,VAR
      INTEGER ::
     &     I,ITM,L,M,MM1,MMAX,MMAXP1,NU,NUA,
     &     NUAMAX,NUB,NUBMAX,NUMAX,NUTM
      CI=DCMPLX(0.D0,1.D0)
      MMAXP1=MMAX+1
C**COSINE AND SINE FUNCTION
      NUMAX=NUAMAX*NUBMAX
      DO 20 NU=1,NUMAX
      ARG=2.0*PI*DFLOAT(NU)/DFLOAT(NUMAX)
      C(NU)=DCOS(ARG)
   20 S(NU)=DSIN(ARG)
   50 PP=0.0
      DO 60 I=1,ITM
   60 PP=PP+SERIES(I)*CONJG(SERIES(I))
      P(1)=PP/DFLOAT(ITM)
      VAR=P(1)
      M=1
      B1(1)=SERIES(1)
      B2(ITM-1)=SERIES(ITM)
      DO 70 I=2,ITM-1
      B1(I)=SERIES(I)
   70 B2(I-1)=SERIES(I)
      GO TO 80
  100 DO 110 I=1,M
  110 AA(I)=A(I)
      M=M+1
      DO 120 I=1,ITM-M
      B1(I)=B1(I)-CONJG(AA(M-1))*B2(I)
  120 B2(I)=B2(I+1)-AA(M-1)*B1(I+1)
   80 ANOM=DCMPLX(0.D0,0.D0)
      ADEN=DCMPLX(0.D0,0.D0)
      DO 90 I=1,ITM-M
      ANOM=ANOM+CONJG(B1(I))*B2(I)
   90 ADEN=ADEN+B1(I)*CONJG(B1(I))+B2(I)*CONJG(B2(I))
      A(M)=(ANOM+ANOM)/ADEN
      P(M+1)=P(M)*(1.0-CONJG(A(M))*A(M))
      IF (M.EQ.1) GO TO 100
  130 CONTINUE
      DO 140 I=1,M-1
  140 A(I)=AA(I)-A(M)*CONJG(AA(M-I))
      IF (M.LT.MMAX) GO TO 100
C**FINAL PREDICTION ERROR
      DO 150 M=1,MMAXP1
  150 FPE(M)=P(M)*DFLOAT(ITM+M-1)/DFLOAT(ITM-M+1)
      DO 180 NUA=1,NUAMAX
      POWERX=0.
C**FREQUENCY BAND AVERAGE
      DO 170 NUB=1,NUBMAX
      NU=NUB+NUA*NUBMAX+(NUMAX-3*NUBMAX-1)/2
      CSUM=1.
      DO 160 M=1,MMAX
      NUTM=MOD(NU*M-1,NUMAX)+1
  160 CSUM=CSUM-A(M)*(C(NUTM)-CI*S(NUTM))
  170 POWERX=POWERX+P(MMAXP1)/(CSUM*CONJG(CSUM))
      POWER(NUA)=.5*POWERX/DFLOAT(NUBMAX)
  180 CONTINUE
      PNU=0.0
      DO 210 L=1,NUAMAX
  210 PNU=PNU+POWER(L)
      PNU=PNU/(.5*NUAMAX)
      RETURN
      END SUBROUTINE MEM


      MODULE BDIJ
!@sum  stores information for outputting lon-lat diagnostics
!@auth M. Kelley

      use MODEL_COM, only : IM,JM
      use DAGCOM, only : kaij,kaijx,lname_ij,name_ij,units_ij
     *  ,kgz_max,pmb,ght

      IMPLICIT NONE

!@param nij_o total number of diagnostic ij-format fields
      integer nij_o

!@var ij_xxx non single-aij diagnostic names
      INTEGER :: ij_topo, ij_jet, ij_wsmn, ij_jetdir, ij_wsdir, ij_grow,
     *  ij_netrdp, ij_albp, ij_albg, ij_albv, ij_ntdsese, ij_fland,
     *  ij_ntdsete, ij_dzt1 

!@param LEGEND "contour levels" for ij-maps
      CHARACTER(LEN=40), DIMENSION(24), PARAMETER :: LEGEND=(/ !
     1  '0=0,1=5...9=45,A=50...K=100             ', ! ir_pct    fac=.2
     2  '0=0...9=90,A=100...I=180...R=270        ', ! ir_angl       .1
     3  '1=.5...9=4.5,A=5...Z=17.5,+=MORE        ', ! ir_0_18        2
     4  '1=.1...9=.9,A=1...Z=3.5,+=MORE          ', ! ir_0_4        10
     5  '1=2...9=18,A=20...Z=70,+=MORE           ', ! ir_0_71       .5
     6  '1=50...9=450,A=500...Z=1750,+=MORE      ', ! ir_0_1775     .02
     7  '1=100...9=900,A=1000...Z=3500,+=MORE    ', ! ir_0_3550     .01
     8  '1=20...9=180,A=200...Z=700,+=MORE       ', ! ir_0_710      .05
     9  'A=1...Z=26,3=30...9=90,+=100-150,*=MORE ', ! ir_0_26_150    1
     O  '0=0,A=.1...Z=2.6,3=3...9=9,+=10-15      ', ! ir_0_3_15     10
     1  '-=LESS,Z=-78...0=0...9=27,+=MORE        ', ! ir_m80_28     .33
     2  '-=LESS,Z=-260...0=0...9=90,+=MORE       ', ! ir_m265_95    .1
     3  '-=LESS,Z=-520...0=0...9=180,+=MORE      ', ! ir_m530_190   .05
     4  '-=LESS,Z=-1300...0=0...9=450,+=MORE     ', ! ir_m1325_475  .02
     5  '-=LESS,Z=-2600...0=0...9=900,+=MORE     ', ! ir_m2650_950  .01
     6  '-=LESS,Z=-3900...0=0...9=1350,+=MORE    ', ! ir_m3975_1425 .007
     7  '-=LESS,Z=-5200...0=0...9=1800,+=MORE    ', ! ir_m5300_1900 .005
     8  '-=LESS,9=-.9...0=0,A=.1...Z=2.6,+=MORE  ', ! ir_m1_3       10
     9  '-=LESS,9=-45...0=0,A=5...I=45...+=MORE  ', ! ir_m45_130    .2
     O  '-=LESS,9=-90...0=0,A=10...Z=260,+=MORE  ', ! ir_m95_265    .1
     1  '-=LESS,9=-180...A=20...Z=520,+=MORE     ', ! ir_m190_530   .05
     2  '-=LESS,9=-9...0=0,A=1...Z=26,+=MORE     ', ! ir_m9_26       1
     3  '-=LESS,9=-36...0=0,A=4...Z=104,+=MORE   ', ! ir_m38_106    .25
     4  '1=5...9=45,A=50...Z=175,+=MORE          ' /) ! ir_0_180    .2
!@var ir_xxxx names for indices to LEGEND indicating the (rounded) range
      integer, parameter :: ir_pct=1, ir_angl=2, ir_0_18=3, ir_0_4=4,
     * ir_0_71=5, ir_0_1775=6, ir_0_3550=7, ir_0_710=8, ir_0_26_150=9,
     * ir_0_3_15=10, ir_m80_28=11, ir_m265_95=12, ir_m530_190=13,
     * ir_m1325_475=14, ir_m2650_950=15, ir_m3975_1425=16,
     * ir_m5300_1900=17, ir_m1_3=18, ir_m45_130=19, ir_m95_265=20,
     * ir_m190_530=21, ir_m9_26=22, ir_m38_106=23, ir_0_180=24
!@var fac_legnd = 1/(range_of_1_colorbox)
      real*8, dimension(24) :: fac_legnd=(/
     1      1d0/5,  1d0/10,    2.d0,   10.d0,   1d0/2,
     6     1d0/50, 1d0/100,  1d0/20,    1.d0,   10.d0,
     1      1d0/3,  1d0/10,  1d0/20,  1d0/50, 1d0/100,
     6    1d0/150, 1d0/200,   10.d0,   1d0/5,  1d0/10,
     1     1d0/20,    1.d0,   1d0/4,   1d0/5           /)

!@param CBAR "color bars" for ij-maps
      CHARACTER(LEN=38), PARAMETER, DIMENSION(5) :: CBAR=(/
     &     ' 0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ+',  ! ib_pos
     &     ' 0123456789ABCDEFGHIJKX               ',  ! ib_pct
     &     '-9876543210ABCDEFGHIJKLMNOPQRSTUVWXYZ+',  ! ib_npp
     &     ' 0ABCDEFGHIJKLMNOPQRSTUVWXYZ3456789+* ',  ! ib_hyb
     &     '-ZYXWVUTSRQPONMLKJIHGFEDCBA0123456789+'/) ! ib_nnp
!@var ib_xxx indices for color bars
      integer, parameter :: ib_pos=1,ib_pct=2,ib_npp=3,ib_hyb=4,ib_nnp=5

!@var SENTDSE stand.eddy northw. transport of dry static energy * 16
!@var TENTDSE trans.eddy northw. transport of dry static energy * 16
      DOUBLE PRECISION, DIMENSION(IM,JM) :: SENTDSE,TENTDSE


      contains

      function mark (val,ibar,undef)
!@sum  mark selects a character (color) based on value and color bar
!@auth R. Ruedy
!@ver  1.0
      real*8 val,undef
      integer ibar,n
      character*1 mark

      select case (ibar)
      case (ib_pct)                                ! 0.....100 %
        n = 2.5 + val
        if (val .ge. 20.) n=23
        if (val .le.  0.) n= 1
        mark = cbar(ib_pct)(n:n)
      case (ib_pos)                                ! 0++++++++++
        n = 2.5 + val
c          non-unif scaling: (currently not used)
c          if (n .gt. 13) n = (n+123)/10
        if (n .gt. 38) n=38
        if (n .lt. 1 .or. val .le. 0.) n= 1
        mark = cbar(ib_pos)(n:n)
      case (ib_npp)                                ! ---0+++++++
        n = 11.5 + val
        if (n .gt. 38.) n=38
        if (n .lt.  1 ) n= 1
        mark = cbar(ib_npp)(n:n)
      case (ib_nnp)                                ! -------0+++
        n = 28.5 + val
        if (n .gt. 38.) n=38
        if (n .lt.  1 ) n= 1
        mark = cbar(ib_nnp)(n:n)
      case (ib_hyb)                      ! hybrid: multiple scales
        n = 2.5 + val
        if (n .gt. 28) n=(n+263)/10
        if (n .gt. 35) n=(n+180)/6
        if (n .gt. 37) n=37
        if (val .le.  0.) n=1
        mark = cbar(ib_hyb)(n:n)
      end select
      if (val .eq. undef) mark=' '

      return
      end function mark

      function ib_of_legnd (leg)
!@sum  ib_of_legnd finds the 'colorbar' for the given legend
!@auth R. Ruedy
!@ver  1.0
      integer ib_of_legnd, leg

      ib_of_legnd = ib_pos
      if (legend(leg)(7:8) .eq. ',Z') ib_of_legnd = ib_nnp
      if (legend(leg)(7:8) .eq. ',9') ib_of_legnd = ib_npp
      if (index(legend(leg)(21:40),'-') .gt. 0) ib_of_legnd = ib_hyb
      if (index(legend(leg),'100 ') .gt. 0) ib_of_legnd = ib_pct

      return
      end function ib_of_legnd

      END MODULE BDIJ


      SUBROUTINE IJ_TITLEX
!@sum  IJ_TITLEX titles etc for composite ij output
!@auth G. Schmidt/M. Kelley
!@ver  1.0
      USE BDIJ
      IMPLICIT NONE
      INTEGER :: k,k1
c
      k = kaij
c
      k = k + 1

      ij_topo = k
      name_ij(k) = 'topog'
      lname_ij(k) = 'TOPOGRAPHY'
      units_ij(k) = 'METERS'

      k = k + 1
      ij_fland = k
      name_ij(k) = 'frac_land'
      lname_ij(k) = 'LAND COVERAGE'
      units_ij(k) = '%'

      k = k + 1
      ij_jet = k
      name_ij(k) = 'jet_speed'
      lname_ij(k) = 'JET SPEED'
      units_ij(k) = 'METERS/SEC'

      k = k + 1
      ij_wsmn = k
      name_ij(k) = 'rt_usmn2_vsmn2'
      lname_ij(k) = 'SURF WIND SPEED FROM Uav,Vav'
      units_ij(k) = 'M/S'

      k = k + 1
      ij_jetdir = k
      name_ij(k) = 'jet_dir'
      lname_ij(k) = 'JET DIRECTION'
      units_ij(k) = 'CW NOR'

      k = k + 1
      ij_wsdir = k
      name_ij(k) = 'srf_wind_dir'
      lname_ij(k) = 'SURFACE WIND DIRECTION'
      units_ij(k) = 'CW NOR'

      k = k + 1
      ij_netrdp = k
      name_ij(k) = 'net_rad_planet'
      lname_ij(k) = 'NET RAD. OF PLANET'
      units_ij(k) = 'WATTS/M**2'

      k = k + 1
      ij_albp = k
      name_ij(k) = 'plan_alb'
      lname_ij(k) = 'PLANETARY ALBEDO'
      units_ij(k) = '%'

      k = k + 1
      ij_albg = k
      name_ij(k) = 'grnd_alb'
      lname_ij(k) = 'GROUND ALBEDO'
      units_ij(k) = '%'

      k = k + 1
      ij_albv = k
      name_ij(k) = 'vis_alb'
      lname_ij(k) = 'VISUAL ALBEDO'
      units_ij(k) = '%'

      k = k + 1
      ij_ntdsese = k
      name_ij(k) = 'stand_eddy_nt_dse'
      lname_ij(k) = 'NT DRY STAT ENR BY ST ED' ! NORTHWD TRANSP
      units_ij(k) = 'E14 WT'

      k = k + 1
      ij_ntdsete = k
      name_ij(k) = 'trans_eddy_nt_dse'
      lname_ij(k) = 'NT DRY STAT ENR BY TR ED' ! NORTHWD TRANSP
      units_ij(k) = 'E14 WT'

      ij_dzt1 = k+1
      do k1 = 1,kgz_max-1
        name_ij(k+k1) = 'dztemp_1000-850'
        if(k1.gt.1) write(name_ij(k+k1)(8:15),
     *    '(i3.3,a1,i3.3,a1)') nint(pmb(k1)),'-',nint(pmb(k1+1)),' '
        lname_ij(k+k1) = 'THICKNESS TEMP 1000-850'
        if(k1.gt.1) write(lname_ij(k+k1)(16:23),
     *     '(i3.3,a1,i3.3,a1)') nint(pmb(k1)),'-',nint(pmb(k1+1)),' '
        units_ij(k+k1) = 'C'
      end do
      k = k + kgz_max -1

      k = k + 1
      ij_grow = k
      name_ij(k) = 'grow_seas'
      lname_ij(k) = 'GROWING SEASON'
      units_ij(k) = 'days'

c Check the count
      if (k .gt. kaijx) then
        write (6,*) 'Increase kaijx=',kaijx,' to at least ',k
        stop 'IJ_TITLES: kaijx too small'
      end if

      END SUBROUTINE IJ_TITLEX


      subroutine IJ_MAPk (k,smap,smapj,gm,jgrid,irange,name,lname,units)
!@sum IJ_MAPk returns the map data and related terms for the k-th field
      USE DAGCOM
      USE CONSTANT, only :
     &     grav,rgas,sday,twopi,sha,kapa,bygrav,tf,undef,teeny
      USE MODEL_COM, only :
     &     im,jm,fim,jeq,byim,DTsrc,ptop,
     &     FLAND,FLICE,FEARTH,FOCEAN,   IDACC,
     &     JHOUR,JHOUR0,JDATE,JDATE0,AMON,AMON0,JYEAR,JYEAR0,
     &     NDAY,Itime,Itime0,XLABEL,LRUNID,ZATMO
      USE DAGPCOM, only :
     *     P1000K
      USE BDIJ

      IMPLICIT NONE

      REAL*8, DIMENSION(IM,JM) :: anum,adenom,smap
      REAL*8, DIMENSION(JM) :: smapj
      integer i,j,k1,k2,k,iwt,jgrid,irange,n1,n2
      character(len=30) name,units
      character(len=80) lname
      real*8 :: gm,nh,sh, off, byiacc, scalek

c**** Find & scale the numerators and find the appropriate denominators
c****
      adenom = 1.                                             ! default
      if (k.eq.ij_dsev) adenom(1,1) = undef    ! zonal mean -> zonal sum

c**** the standard cases: aij(.,.,k) or aij(.,.,k)/aij(.,.,k1)
      if (k .le. kaij) then
        name = name_ij(k) ; lname = lname_ij(k) ; units = units_ij(k)
        iwt = iw_ij(k) ; jgrid = jgrid_ij(k) ; irange = ir_ij(k)
c**** offsets ("  + " or "  - " in lname_ij, i.e. 2 blanks,+|-,1 blank)
        off = 0.
        k1 = index(lname_ij(k),'  - ')
        if (k1 .le. 0) k1 = index(lname_ij(k),'  + ')
        if (k1 .gt. 0) then
          if (index(lname_ij(k),'  + TF ') .gt. 0) then
            off = TF                       ! should do in accum-phase ??
          else if (index(lname_ij(k),'  - PTOP') .gt. 0)  then
            off = -ptop
          end if
          lname_ij(k)(k1:80) = ' ' ; lname = lname_ij(k)
        end if
        byiacc = 1./(idacc(ia_ij(k))+teeny)
        do j=1,jm
        do i=1,im
          anum(i,j) = aij(i,j,k)*(scale_ij(k)*byiacc) - off
        end do
        end do

c**** ratios (the denominators)
        k1 = index(lname_ij(k),' x ')
        if (k1 .gt. 0) then
          if (index(lname_ij(k),' x POCEAN') .gt. 0) then
            do j=1,jm
            do i=1,im
              adenom(i,j) = 1-fland(i,j) - aij(i,j,ij_rsoi)*byiacc
            end do
            end do
          else if (index(lname_ij(k),' x POICE') .gt. 0) then
            do j=1,jm
            do i=1,im
              adenom(i,j)=aij(i,j,ij_rsoi) * byiacc
            end do
            end do
          else if (index(lname_ij(k),' x PLICE') .gt. 0) then
            adenom = flice
          else if (index(lname_ij(k),' x PSOIL') .gt. 0) then
            adenom = fearth
          else if (index(lname_ij(k),' x TOTAL CLOUD') .gt. 0) then
            do j=1,jm
            do i=1,im
              adenom(i,j)=aij(i,j,ij_cldcv) * byiacc
            end do
            end do
          end if
          lname_ij(k)(k1:80) = ' '    ; lname = lname_ij(k)
        end if
        go to 100
      end if

c**** compound quantities defined with their attributes (k > kaij)
c****
      iwt = iw_all ; jgrid = 1 ; irange = ir_pct    ! defaults
      name  = name_ij(k)
      lname = lname_ij(k) ; units = units_ij(k)

c**** time independent arrays
      if      (k.eq.ij_topo)  then
        anum = zatmo*bygrav    ; irange = ir_0_3550

      else if (k.eq.ij_fland) then
        anum = 100.*fland

c**** vectors: magnitude
      else if (k.eq.ij_jet.or.k.eq.ij_wsmn) then
          jgrid = 2 ;  n1 = ij_ujet ; n2 = ij_vjet ; irange = ir_0_71
        if (k.eq.ij_wsmn) then
          jgrid = 1 ;  n1 = ij_us   ; n2 = ij_vs   ; irange = ir_0_18
        end if
        byiacc=1./(idacc(ia_ij(n1))+teeny)
        do j=1,jm
        do i=1,im
          anum(i,j)=sqrt(aij(i,j,n1)**2+aij(i,j,n2)**2)*
     *      (scale_ij(n1)*byiacc)
        end do
        end do

c**** vectors: direction clockwise north (-180 -> 180)
      else if (k.eq.ij_jetdir .or. k.eq.ij_wsdir) then
        irange = ir_angl
        if (k.eq.ij_jetdir) then
          jgrid = 2 ; n1 = ij_ujet ; n2 = ij_vjet
        else if (k.eq.ij_wsdir)  then
          jgrid = 1 ; n1 = ij_us   ; n2 = ij_vs
        end if
        do j=2,jm
        do i=1,im
          anum(i,j)=360.*atan2(aij(i,j,n1)+teeny,aij(i,j,n2)+teeny)
     *         /twopi
        end do
        end do

c**** linear combinations (sums, differences, etc)
      else if (k.eq.ij_netrdp) then
        irange = ir_m530_190
        byiacc = 1./(idacc(ia_ij(ij_trnfp0))+teeny)
        do j=1,jm
        do i=1,im
          anum(i,j)=(aij(i,j,ij_trnfp0)+aij(i,j,ij_srnfp0))*byiacc
        end do
        end do

c**** ratios of lin. comb.: albedos from net radiation
      else if (k.eq.ij_albp .or. k.eq.ij_albg) then
        n2=ij_srincp0
        if (k.eq.ij_albg) n2=ij_srincg
        n1=ij_srnfp0
        if (k.eq.ij_albg) n1=ij_srnfg
        do j=1,jm
        do i=1,im
          adenom(i,j)=aij(i,j,n2)
          anum(i,j)=100.*(adenom(i,j)-aij(i,j,n1))
        end do
        end do

c**** ratios: albedos from reflected radiation
      else if (k.eq.ij_albv) then
        do j=1,jm
        do i=1,im
          adenom(i,j)=aij(i,j,ij_srincp0)
          anum(i,j)=100.*aij(i,j,ij_srref)
        end do
        end do

c**** precomputed fields: northward tranports by eddies
      else if (k.eq.ij_ntdsese) then                   ! standing eddies
        byiacc=1./(idacc(ia_ij(ij_dsev))+teeny)   ; irange = ir_m95_265
        anum=SENTDSE*(byiacc*scale_ij(ij_dsev))  ;  jgrid = 2
        adenom(1,1) = undef

      else if (k.eq.ij_ntdsete) then                  ! transient eddies
        byiacc=1./(idacc(ia_ij(ij_dsev))+teeny)   ; irange = ir_m1_3
        anum=TENTDSE*(byiacc*scale_ij(ij_dsev))  ;  jgrid = 2
        adenom(1,1) = undef

c**** group of kgz_max-1 thickness temperatures (from heights)
      else if (k.ge.ij_dzt1 .and. k.le.ij_dzt1+kgz_max-2) then
        byiacc = 1./(idacc(ia_ij(ij_phi1k))+teeny) ; irange = ir_m80_28
        k1 = k-ij_dzt1+1  ; k2 = ij_phi1k + k1
        scalek = 1./(rgas*log(pmb(k1)/pmb(k1+1)))
        off = (ght(k1+1)-ght(k1)) * grav
        do j=1,jm
        do i=1,im
          anum(i,j)=((aij(i,j,k2)-aij(i,j,k2-1))*byiacc + off)*scalek-tf
        end do
        end do

c**** length of growing season   (not quite right ???)
      else if (k.eq.ij_grow) then
        byiacc = 1./(idacc(ia_inst)+teeny) ; irange = ir_0_180
        do j=1,jm
        do i=1,im
          anum(i,j)=(tsfrez(i,j,2)-tsfrez(i,j,1))*byiacc
        end do
        end do

      else  ! should not happen
        write (6,*) 'no field defined for ij_index',k
        stop 'ij_mapk: undefined extra ij_field'
      end if

c**** Find final field and zonal, global, and hemispheric means
  100 call ij_avg (anum,adenom,wt_ij(1,1,iwt),jgrid,       ! in
     *             smap,smapj,gm,nh,sh)                    ! out

c**** fill in some key numbers
      if (k .eq. IJ_RSIT) call keyij(gm,nh)

      return

      end subroutine ij_mapk


      subroutine ij_avg (anum,aden,wtij,jgrid, smap,smapj,gm,nh,sh)
!@sum ij_avg finds num/den and various averages from num and den
!@auth R.Ruedy
!@ver  1.0
      USE MODEL_COM, only :  im,jm,fim,jeq
      USE CONSTANT, only :  undef
      USE GEOM, only : wtj

      IMPLICIT NONE

      real*8, dimension(im,jm) :: anum,aden,wtij,smap
      real*8, dimension(jm) :: smapj
      real*8, dimension(2) :: znumh,zdenh
      real*8  gm,nh,sh, sumj,wt
      integer k,i,j,jgrid,jhemi,jmin,jmax,sum_avg

c**** find smap,smapj  from the numerators and denominators
      sum_avg = 2  ! find zonal/global means
      if (aden(1,1).eq.undef) sum_avg = 1  ! find zonal/global sums
      smap = undef ; smapj = undef
      znumh = 0. ; zdenh = 0.
      do jhemi=1,2
        jmin=jeq
        jmax=jm
        if (jhemi.eq.1) then
          jmin=jgrid          !  1 or 2
          jmax=jmin+jeq-2
        end if
        do j=jmin,jmax
          sumj = 0. ; wt = 0.
          do i=1,im
            sumj = sumj + anum(i,j)*wtij(i,j)
            wt  = wt  + aden(i,j)*wtij(i,j)
            if (aden(i,j)*wtij(i,j).gt.0.) smap(i,j)=anum(i,j)/aden(i,j)
            if (sum_avg .eq. 1) smap(i,j) = anum(i,j)
          end do
          znumh(jhemi) = znumh(jhemi) + sumj*wtj(j,sum_avg,jgrid)
          zdenh(jhemi) = zdenh(jhemi) + wt*wtj(j,sum_avg,jgrid)
          if (wt .gt. 0.) smapj(j) = sumj/wt
          if (aden(1,1).eq.undef) smapj(j) = sumj
        end do
      end do
c**** find hemispheric and global means
      nh = undef ; sh = undef ; gm = undef
      if (zdenh(1).gt.0.) sh = znumh(1)/zdenh(1)
      if (zdenh(2).gt.0.) nh = znumh(2)/zdenh(2)
      if (zdenh(1)+zdenh(2).gt.0.) gm = (znumh(1)+znumh(2))/
     /                                   (zdenh(1)+zdenh(2))
      if (aden(1,1).eq.undef) then
        sh = znumh(1) ; nh = znumh(2) ; gm = (znumh(1)+znumh(2))
      end if

      return
      end subroutine ij_avg


      SUBROUTINE DIAGIJ
!@sum  DIAGIJ produces lat-lon fields as maplets (6/page) or full-page
!@+    digital maps, and binary (netcdf etc) files (if qcheck=true)
!@auth Gary Russell,Maxwell Kelley,Reto Ruedy
!@ver   1.0
c     USE DAGCOM, only : QCHECK,acc_period,iu_ij
      USE CONSTANT, only : sha,teeny
      USE MODEL_COM, only :
     &     im,jm,lm,byim,
     &     FLAND,FLICE,FEARTH,FOCEAN,vdata,
     &     JHOUR,JHOUR0,JDATE,JDATE0,AMON,AMON0,JYEAR,JYEAR0,
     &     NDAY,Itime,Itime0,XLABEL,LRUNID
      USE LAKES_COM, only : flake
      USE GEOM, only : DXV
c     USE DAGCOM, only :
c    &     ajk,kdiag,aij,kaij,aijk,tsfrez,ia_ij,
c    *     ijk_v,ijk_dse,ijk_u,ijk_dp,
c    &     IJ_DSEV,IJ_TRNFP0,IJ_SRNFP0,IJ_SLP,IJ_TS !not a generic subr.
      USE DAGCOM
      USE BDIJ

      IMPLICIT NONE

!@var Qk: if Qk(k)=.true. field k still has to be processed
      logical, dimension (kaijx) :: Qk
!@var Iord: Index array, fields are processed in order Iord(k), k=1,2,..
!@+     only important for fields 1->nmaplets+nmaps (appear in printout)
!@+     Iord(k)=0 indicates that a blank space replaces a maplet
      INTEGER Iord(kaijx+10),nmaplets,nmaps ! 10 extra blank maplets
      INTEGER kmaplets
      REAL*8, DIMENSION(IM,JM) :: SMAP
      REAL*8, DIMENSION(JM) :: SMAPJ
      CHARACTER xlb*32,title*48,lname*80,name*30,units*30
!@var LINE virtual half page (with room for overstrikes)
      CHARACTER*133 LINE(53)

      INTEGER ::   I,J,K,L,M,N,kcolmn,nlines,jgrid,irange,iu_Iij

      DOUBLE PRECISION ::
     &     DAYS,ZNDE16,ZS,DPTI,PVTI,gm,
     &     DE4TI,BYDPK,SZNDEG

C**** OPEN PLOTTABLE OUTPUT FILE IF DESIRED
      IF(QCHECK) call open_ij(trim(acc_period)//'.ij'//XLABEL(1:LRUNID))

C**** INITIALIZE CERTAIN QUANTITIES
      call ij_titlex
C**** standard printout
      kmaplets = 51
      nmaplets = kmaplets+(kgz_max-1)*2 ; nmaps = 2
      iord(1:kmaplets) = (/
     *  ij_topo,    ij_fland,   ij_rsoi,     ! pg  1  row 1
     *  ij_rsnw,    ij_snow,    ij_rsit,     !        row 2
     *  ij_prec,    ij_evap,    ij_shdt,     ! pg  2  row 1
     *  ij_beta,    ij_rune,    ij_tg1,      !        row 2
     *  ij_ws,      ij_jet ,    ij_wsmn,     ! pg  3  row 1
     *  ij_sdrag,   ij_jetdir,  ij_wsdir,    !        row 2
     *  ij_cldcv,   ij_pmccld,  ij_cldtppr,  ! pg  4  row 1
     *  ij_pcldl,   ij_pcldm,   ij_pcldh,    !        row 2
     *  ij_netrdp,  ij_srtr,    ij_btmpw,    ! pg  5  row 1
     *  ij_albp,    ij_albg,    ij_albv,     !        row 2
     *  ij_trnfp0,  ij_neth,    ij_dtdp,     ! pg  6  row 1
     *  ij_dsev,    ij_ntdsese, ij_ntdsete,  !        row 2
     *  ij_gwtr,    ij_wmsum,   ij_dcnvfrq,  ! pg  7  row 1
     *  ij_scnvfrq, ij_pdcld,   ij_pscld,    !        row 2
     *  ij_gw1,     ij_gw2,     ij_gw3,      ! pg  8  row 1
     *  ij_gw4,     ij_gw5,     ij_gw6,      !        row 2
     *  ij_gw7,     ij_gw8,     ij_gw9/)     ! pg  9  row 1

C**** Fill in maplet indices for geoptential heights and thickness T's
      do k=1,kgz_max-1
        iord(k+kmaplets) = ij_phi1k+k  !i.e. first entry is ij_phi850
        iord(k+kmaplets+kgz_max-1) = ij_dzt1+k-1
      end do

C**** Add the full-page maps (nmaps)
      iord(nmaplets+1:nmaplets+nmaps) = (/ij_slp,ij_ts/)

c**** always skip unused fields
      Qk = .true.
      do k=1,kaijx
        if (index(lname_ij(k),'unused').gt.0) Qk(k) = .false.
      end do

      if (kdiag(3).gt.0) call set_ijout (kdiag(3),nmaplets,nmaps,
     *  Iord,Qk,iu_Iij)

      xlb=acc_period(1:3)//' '//acc_period(4:12)//' '//XLABEL(1:LRUNID)
C****
      DAYS=(Itime-Itime0)/DFLOAT(nday)
C**** Collect the appropriate weight-arrays in WT_IJ
      do j=1,jm
      do i=1,im
        wt_ij(i,j,1) = 1.
        wt_ij(i,j,2) = focean(i,j)
        wt_ij(i,j,3) = flake(i,j)
        wt_ij(i,j,4) = flice(i,j)
        wt_ij(i,j,5) = fearth(i,j)
        wt_ij(i,j,6) = fearth(i,j)*(vdata(i,j,1)+vdata(i,j,10))
        wt_ij(i,j,7) = fearth(i,j)*(1.-(vdata(i,j,1)+vdata(i,j,10)))
      end do
      end do
C**** CACULATE STANDING AND TRANSIENT EDDY NORTHWARD TRANSPORT OF DSE
      SENTDSE = 0
      TENTDSE = 0
      DO J=2,JM
      DO K=1,LM
        DPTI=0.
        PVTI=0.
        DE4TI=0.
        DO I=1,IM
          IF (AIJK(I,J,K,IJK_DP).GT.0.) THEN
            DPTI=DPTI+AIJK(I,J,K,IJK_DP)
            BYDPK=1./(AIJK(I,J,K,IJK_DP)+teeny)
            PVTI=PVTI+AIJK(I,J,K,IJK_V)
            DE4TI=DE4TI+AIJK(I,J,K,IJK_DSE)
            SENTDSE(I,J)=SENTDSE(I,J)
     *        +(AIJK(I,J,K,IJK_DSE)*AIJK(I,J,K,IJK_V)*BYDPK)
          END IF
        END DO
        SZNDEG=DE4TI*PVTI/(DPTI+teeny)
        DO I=1,IM
          SENTDSE(I,J)=SENTDSE(I,J)-SZNDEG*byim
        END DO
      END DO
      END DO
      DO J=2,JM
        ZNDE16=0.
        DO L=1,LM
          ZNDE16=ZNDE16+(SHA*AJK(J,L,JK_ZMFNTSH)+AJK(J,L,JK_ZMFNTGEO))
        END DO
        ZNDE16=4.*ZNDE16*DXV(J)*byim
        DO I=1,IM
          SENTDSE(I,J)=4.*SENTDSE(I,J)*DXV(J)
          TENTDSE(I,J)=AIJ(I,J,IJ_DSEV)-ZNDE16-SENTDSE(I,J)
        END DO
      END DO

C**** Fill in the undefined pole box duplicates
      DO N=1,KAIJ
      IF (JGRID_ij(N).EQ.2) CYCLE
      DO I=1,IM
        AIJ(I,1,N)=AIJ(1,1,N)
        AIJ(I,JM,N)=AIJ(1,JM,N)
      END DO
      END DO

C**** Print out 6-map pages
      do n=1,nmaplets
        if (mod(n-1,6) .eq. 0) then
c**** print header lines
          WRITE (6,901) XLABEL
          WRITE (6,902) JYEAR0,AMON0,JDATE0,JHOUR0,
     *      JYEAR,AMON,JDATE,JHOUR,ITIME,DAYS
        end if
        kcolmn = 1 + mod(n-1,3)
        if (kcolmn .eq. 1) line=' '
c**** Find, then display the appropriate array
        k = Iord(n)
        if (k .gt. 0 .and. Qk(k)) then
          call ij_mapk (k,smap,smapj,gm,jgrid,irange,name,lname,units)
          title=trim(lname)//' ('//trim(units)//')'
          call cmap(smap,smapj,gm,irange,title,line,kcolmn,nlines)
          if(qcheck) call pout_ij(title//xlb,name,lname,units,
     *                            smap,smapj,gm,jgrid)
          Qk(k) = .false.
        end if
c**** copy virtual half-page to paper if appropriate
        if (kcolmn.eq.3 .or. n.eq.nmaplets) then
          do k=1,nlines
            write (6,'(a133)') line(k)
          end do
        end if
      end do

C**** Print out full-page digital maps
      do n=nmaplets+1,nmaplets+nmaps
        k = Iord(n)
      if (k.le.0 .or. .not.Qk(k)) cycle
        call ij_mapk (k,smap,smapj,gm,jgrid,irange,name,lname,units)
        title=trim(lname)//' ('//trim(units)//')'
        call ijmap (title//xlb,smap,smapj,jgrid)
        if(qcheck) call pout_ij(title//xlb,name,lname,units,smap,smapj,
     *                          gm,jgrid)
        Qk(k) = .false.
      end do

      if (.not.qcheck) RETURN
C**** produce binary files of remaining fields if appropriate
      do k=1,kaijx
        if (Qk(k)) then
          call ij_mapk (k,smap,smapj,gm,jgrid,irange,name,lname,units)
          title=trim(lname)//' ('//trim(units)//')'
          call pout_ij(title//xlb,name,lname,units,smap,smapj,gm,jgrid)
        end if
      end do
      call close_ij
      if (kdiag(3).lt.7) CALL IJKMAP (iu_Iij)

      RETURN
C****
  901 FORMAT ('1',A)
  902 FORMAT ('0',15X,'From:',I6,A6,I2,',  Hr',I3,
     *  6X,'To:',I6,A6,I2,', Hr',I3,'  Model-Time:',I9,5X,
     *  'Dif:',F7.2,' Days')
      END SUBROUTINE DIAGIJ


      subroutine cmap (smap,smapj,gm,irange,title,line,kcol,nlines)
!@sum cmap prints a maplet onto 1/3 of a virtual half-page line(1-51)
!@auth R.Ruedy
!@ver  1.0
      use model_com, only : im,jm,fland
      use bdij
      use constant, only : undef
      use dagpcom, only : inci,incj

      IMPLICIT NONE

      real*8, dimension(im,jm) :: smap
      real*8, dimension(jm) :: smapj
      real*8  gm,val,zmax,zmin
      integer irange,kcol,k,k1,ifm,nz1,nz2,n1,n2,ibar,i,j
     *  ,nlines
      character(len=133), dimension(53) :: line
      character(len=36) :: title

c**** find first and last column for zonal means nz1,nz2 and maps n1,n2
      nz1 = 2 + (kcol-1)*(9+im/inci) ; nz2 = nz1 + 4
      n1  = nz2 + 2  ;  n2 = n1 + im/inci - 1
c**** pick color bar and format for zonal mean according to the range
      ibar = ib_of_legnd(irange)
      zmax = 999.9 ; zmin = -99.9
      if (kcol.gt.1) then
        zmax = 99999.9 ; zmin = -9999.9
        nz1 = nz1 - 2
      end if
      ifm = 2
      do j=1,jm
        if (smapj(j).eq.undef) cycle
        if (smapj(j).gt.zmax .or. smapj(j).lt.zmin) ifm = 1
      end do

c**** title on line 1
      line(1)(n1:n2) = title(1:36) ; line(1)(1:1) = '0'  ! line feed
c**** use 2 lines for each shown latitude because of overstrike
      k=0
      do j=jm,1,-incj
        k = k+2
c**** zonal mean
        val=0.
        if (smapj(j) .ne. undef) val = smapj(j)
        if (kcol.eq.1) then
          if (ifm.eq.1) write(line(k)(nz1:nz2),'(i5)') nint(val)
          if (ifm.eq.2) write(line(k)(nz1:nz2),'(f5.1)') val
        else
          if (ifm.eq.1) write(line(k)(nz1:nz2),'(i7)') nint(val)
          if (ifm.eq.2) write(line(k)(nz1:nz2),'(f7.1)') val
        end if
c**** mark selected longitudes for each selected latitude
        k1 = n1-1
        do i=1,im,inci
          k1 = k1+1
          val = undef
          if (smap(i,j).ne.undef) val = smap(i,j)*fac_legnd(irange)
c**** for angles, change range from -180->180 to 0-360 (fac=.1)
          if (irange.eq.ir_angl .and. val.lt.-.5) val = val+36.
          if (irange.eq.ir_angl .and. val.le..1) val = .1
          line(k)(k1:k1) = mark(val,ibar,undef)
          if (fland(i,j).gt..5) line(k+1)(k1:k1) = line(k)(k1:k1)
          line(k+1)(1:1) = '+'             !  overstrike
        end do
      end do
c**** below map, show global mean and mark each quarter with a '+'
      k = k+2
      if (gm .ge. zmin .and. gm .le. zmax) then
        if (kcol.eq.1) write(line(k)(nz1:nz2),'(f5.1)') gm
        if (kcol.gt.1) write(line(k)(nz1:nz2),'(f7.1)') gm
      else
        if (kcol.eq.1) write(line(k)(nz1:nz2),'(i5)') nint(gm)
        if (kcol.gt.1) write(line(k)(nz1:nz2),'(i7)') nint(gm)
      end if
      do k1 = n1,n2,im/(4*inci)
        line(k)(k1:k1) = '+'
      end do
c**** last line: legend (there is a little more space in col 1 and 2)
      if(kcol.lt.3) n2 = n1 + 39
      line(k+1)(n1:n2) = legend(irange)
      nlines = k+1

      return
      end subroutine cmap


      SUBROUTINE IJMAP (title,smap,smapj,jgrid)
      USE MODEL_COM, only :
     &     im,jm,FLAND,NDAY,JHOUR,JHOUR0,JDATE,JDATE0,AMON,AMON0,
     &     JYEAR,JYEAR0,Itime,Itime0,XLABEL,lrunid
      USE GEOM, only :
     &     LAT_DG,LON_DG
      use dagpcom, only : inc=>inci
      IMPLICIT NONE

      CHARACTER*48 TITLE

      CHARACTER(LEN=3), DIMENSION(IM) :: LINE
      CHARACTER(LEN=9) :: AVG

      REAL*8, DIMENSION(IM,JM) :: SMAP
      REAL*8, DIMENSION(JM) :: SMAPJ
      DOUBLE PRECISION :: DAYS
      INTEGER :: I,J,jgrid

C**** WRITE HEADER LINES
      DAYS=(Itime-Itime0)/DFLOAT(nday)
      WRITE(6,901)XLABEL
      WRITE (6,902) JYEAR0,AMON0,JDATE0,JHOUR0,
     *  JYEAR,AMON,JDATE,JHOUR,ITIME,DAYS
      WRITE(6,900) TITLE(1:48)
      DO I=1,IM
        WRITE(LINE(I),'(I3)') I
      end do
      AVG='     MEAN'
      WRITE (6,910) (LINE(I),I=1,IM,INC),AVG
      WRITE(6,940)
      WRITE(6,940)

C**** PRINT MAP
      DO J=JM,jgrid,-1
        do i=1,im
          IF (SMAP(I,J).LT.999.5.AND.SMAP(I,J).GE.-99.5) then
            line(i) = '   '
            write (line(i),'(i3)') nint(SMAP(I,J))
          else
            line(i)='***'
          end if
        end do
        WRITE(AVG,'(F9.2)') SMAPJ(J)
        WRITE (6,920) NINT(LAT_DG(J,jgrid)),J,(LINE(I),I=1,IM,INC),AVG
        DO I=1,IM
          IF (FLAND(I,J).lt..5) LINE(I)='   '
        end do
        WRITE (6,925) (LINE(I),I=1,IM,INC)
        WRITE (6,925) (LINE(I),I=1,IM,INC)
        IF (JM.LE.24) WRITE (6,940)
      END DO
      WRITE (6,930) (LON_DG(I,jgrid),I=1,IM,INC*2)
      RETURN
C****
  900 FORMAT('0',45X,A48)
  901 FORMAT ('1',A)
  902 FORMAT ('0',15X,'From:',I6,A6,I2,',  Hr',I3,
     *  6X,'To:',I6,A6,I2,', Hr',I3,'  Model-Time:',I9,5X,
     *  'Dif:',F7.2,' Days')
  910 FORMAT('0LAT  J/I  ',36A3,A9)
  920 FORMAT(2I4,3X,36A3,A9)
  925 FORMAT('+',10X,36A3)
  930 FORMAT('0  LONG',2X,19F6.1)
  940 FORMAT(' ')
      END SUBROUTINE IJMAP


      subroutine set_ijout (iopt,nmaplets,nmaps,Iord,Qk,iu_Iij)
!@sum set_ijout either lists or sets the fields to be processed
!@auth Reto A. Ruedy
!@ver  1.0
      USE DAGCOM
      USE BDIJ
      use filemanager

      IMPLICIT NONE
      character*80 line
      logical Qk(kaijx)
      INTEGER Iord(kaijx+10),nmaplets,nmaps,iopt,iu_Iij,k,
     *   n,kmap(3)

c**** Just list what's available - then do same for ijk-fields ???
      if (iopt .eq. 8) then
        call openunit('Iij',iu_Iij,.false.,.false.)
        write (iu_Iij,'(a)') 'List of fields shown as maplets'
        do n=1,nmaplets
          k = Iord(n)
          Qk(k) = .false.
          if (k.le.0) then
             write (iu_Iij,'(i3,1x,a)') k, '  blank maplet'
          else
             write (iu_Iij,'(i3,1x,a)') k,lname_ij(k)
          end if
        end do
        write (iu_Iij,'(a)') 'List of fields shown as 1-pg maps'
        do n=nmaplets+1,nmaplets+nmaps
          k = Iord(n)
          Qk(k) = .false.
          if (k.le.0) then
             cycle
          else
             write (iu_Iij,'(i3,1x,a)') k,lname_ij(k)
          end if
        end do
        write (iu_Iij,'(a)') 'List of other fields in binary output'
        do k=1,kaijx
          if (.not.Qk(k)) cycle
          write (iu_Iij,'(i3,1x,a)') k,lname_ij(k)
        end do
        Qk = .false.
        CALL IJKMAP (iu_Iij)
      end if

c**** Redefine nmaplets,nmaps,Iord,Qk if 0 < kdiag(3) < 8
      call openunit('Iij',iu_Iij,.false.,.true.)

      nmaplets = 0 ; nmaps = 0 ; Iord = 0 ; Qk = .false.

      kmap = 0 ; n=0 ; k=0
   10 read (iu_Iij,'(a)',end=20) line
      if (line(1:1) .eq. 'l') go to 20
      if (line(1:1) .eq. 'L') then
        n=n+1
        go to 10
      end if
      k = k+1
      read(line,'(i3)') Iord(k)
      if (Iord(k).gt.0) Qk(Iord(k)) = .true.
      kmap(n) = kmap(n) + 1
      go to 10

   20 nmaplets = kmap(1) ; nmaps = kmap(2)
      if (.not.qcheck .or. kdiag(3).eq.1) call closeunit(iu_Iij)
      return
      end subroutine set_ijout


      SUBROUTINE DIAGCP
!@sum  DIAGCP produces tables of the conservation diagnostics
!@auth Gary Russell/Gavin Schmidt
!@ver  1.0
c      USE PRTCOM, only :
      USE CONSTANT, only :
     &     twopi
      USE MODEL_COM, only :
     &     jm,fim,idacc,jhour,jhour0,jdate,jdate0,amon,amon0,
     &     jyear,jyear0,nday,jeq,itime,itime0,xlabel
      USE GEOM, only :
     &     areag,dlon,dxyp,dxyv,LAT_DG
      USE DAGCOM, only :
     &     consrv,kcon,scale_con,title_con,nsum_con,ia_con,kcmx
      use dagpcom, only : inc=>incj
      IMPLICIT NONE

      INTEGER, DIMENSION(JM) :: MAREA
      DOUBLE PRECISION, DIMENSION(KCON) :: FGLOB
      DOUBLE PRECISION, DIMENSION(2,KCON) :: FHEM
      INTEGER, DIMENSION(JM,KCON) :: MLAT
      DOUBLE PRECISION, DIMENSION(JM+3,KCON) :: CNSLAT
      CHARACTER*4, PARAMETER :: HEMIS(2) = (/' SH ',' NH '/),
     *     DASH = ('----')
      DOUBLE PRECISION, PARAMETER :: XWON = TWOPI/(DLON*FIM)

      INTEGER :: j,jhemi,jnh,jp1,jpm,jsh,jv1,jvm,jx,n
      DOUBLE PRECISION :: aglob,ahem,feq,fnh,fsh,days
C**** CALCULATE SCALING FACTORS
      IF (IDACC(12).LT.1) IDACC(12)=1
C**** CALCULATE SUMMED QUANTITIES
C**** LOOP BACKWARDS SO THAT INITIALISATION IS DONE BEFORE SUMMATION!
      DO J=1,JM
        DO N=KCMX,1,-1
          IF (NSUM_CON(N).eq.0) THEN
            CONSRV(J,N)=0.
          ELSEIF (NSUM_CON(N).gt.0) THEN
            CONSRV(J,NSUM_CON(N))=CONSRV(J,NSUM_CON(N))+CONSRV(J,N)
     *           *SCALE_CON(N)*IDACC(12)/(IDACC(IA_CON(N))+1d-20)
          END IF
        END DO
      END DO
C**** CALCULATE FINAL ANGULAR MOMENTUM + KINETIC ENERGY ON VELOCITY GRID
      DO N=1,23
        FEQ=CONSRV(JEQ,N)*SCALE_CON(N)/(IDACC(IA_CON(N))+1d-20)
        FGLOB(N)=FEQ
        FHEM(1,N)=.5*FEQ
        FHEM(2,N)=.5*FEQ
        CNSLAT(JEQ,N)=FEQ/(FIM*DXYV(JEQ))
        DO JSH=2,JEQ-1
          JNH=2+JM-JSH
          FSH=CONSRV(JSH,N)*SCALE_CON(N)/(IDACC(IA_CON(N))+1d-20)
          FNH=CONSRV(JNH,N)*SCALE_CON(N)/(IDACC(IA_CON(N))+1d-20)
          FGLOB(N)=FGLOB(N)+(FSH+FNH)
          FHEM(1,N)=FHEM(1,N)+FSH
          FHEM(2,N)=FHEM(2,N)+FNH
          CNSLAT(JSH,N)=FSH/(FIM*DXYV(JSH))
          CNSLAT(JNH,N)=FNH/(FIM*DXYV(JNH))
        END DO
        FGLOB(N)=FGLOB(N)/AREAG
        FHEM(1,N)=FHEM(1,N)/(.5*AREAG)
        FHEM(2,N)=FHEM(2,N)/(.5*AREAG)
      END DO
C**** CALCULATE ALL OTHER CONSERVED QUANTITIES ON TRACER GRID
      DO N=24,KCMX
        FGLOB(N)=0.
        FHEM(1,N)=0.
        FHEM(2,N)=0.
        DO JSH=1,JEQ-1
          JNH=1+JM-JSH
          FSH=CONSRV(JSH,N)*SCALE_CON(N)/(IDACC(IA_CON(N))+1d-20)
          FNH=CONSRV(JNH,N)*SCALE_CON(N)/(IDACC(IA_CON(N))+1d-20)
          FGLOB(N)=FGLOB(N)+(FSH+FNH)*DXYP(JSH)
          FHEM(1,N)=FHEM(1,N)+FSH*DXYP(JSH)
          FHEM(2,N)=FHEM(2,N)+FNH*DXYP(JNH)
          CNSLAT(JSH,N)=FSH/FIM
          CNSLAT(JNH,N)=FNH/FIM
        END DO
        FGLOB(N)=FGLOB(N)/AREAG
        FHEM(1,N)=FHEM(1,N)/(.5*AREAG)
        FHEM(2,N)=FHEM(2,N)/(.5*AREAG)
      END DO
      AGLOB=1.D-10*AREAG*XWON
      AHEM=1.D-10*(.5*AREAG)*XWON
C**** LOOP OVER HEMISPHERES
      DAYS=(Itime-Itime0)/DFLOAT(nday)
      DO N=1,KCMX
        DO J=1,JM
          MLAT(J,N)=NINT(CNSLAT(J,N))
        END DO
        CNSLAT(JM+1,N)=FHEM(1,N)
        CNSLAT(JM+2,N)=FHEM(2,N)
        CNSLAT(JM+3,N)=FGLOB(N)
      END DO
      DO JHEMI=2,1,-1
        WRITE (6,901) XLABEL
        WRITE (6,902) JYEAR0,AMON0,JDATE0,JHOUR0,
     *       JYEAR,AMON,JDATE,JHOUR,ITIME,DAYS
        JP1=1+(JHEMI-1)*(JEQ-1)
        JPM=JHEMI*(JEQ-1)
        JV1=2+(JHEMI-1)*(JEQ-2)
        JVM=JEQ+(JHEMI-1)*(JEQ-2)
C**** PRODUCE TABLES FOR ANGULAR MOMENTUM AND KINETIC ENERGY
        WRITE (6,903) (DASH,J=JV1,JVM,INC)
        WRITE (6,904) HEMIS(JHEMI),(NINT(LAT_DG(JX,2)),JX=JVM,JV1,-INC)
        WRITE (6,903) (DASH,J=JV1,JVM,INC)
        DO N=1,23
          WRITE (6,905) TITLE_CON(N),FGLOB(N),FHEM(JHEMI,N),
     *         (MLAT(JX,N),JX=JVM,JV1,-INC)
        END DO
        DO J=JV1,JVM
          MAREA(J)=1.D-10*XWON*FIM*DXYV(J)+.5
        END DO
        WRITE (6,906) AGLOB,AHEM,(MAREA(JX),JX=JVM,JV1,-INC)
C**** PRODUCE TABLES FOR OTHER CONSERVED QUANTITIES
        WRITE (6,907)
        WRITE (6,903) (DASH,J=JP1,JPM,INC)
        WRITE (6,904) HEMIS(JHEMI),(NINT(LAT_DG(JX,1)),JX=JPM,JP1,-INC)
        WRITE (6,903) (DASH,J=JP1,JPM,INC)
        DO N=24,KCMX
          WRITE (6,905) TITLE_CON(N),FGLOB(N),FHEM(JHEMI,N),
     *         (MLAT(JX,N),JX=JPM,JP1,-INC)
        END DO
        DO J=JP1,JPM
          MAREA(J)=1.D-10*XWON*FIM*DXYP(J)+.5
        END DO
        WRITE (6,906) AGLOB,AHEM,(MAREA(JX),JX=JPM,JP1,-INC)
      END DO
      RETURN
C****
  901 FORMAT ('1',A)
  902 FORMAT ('0Conservation Quantities       From:',
     *  I6,A6,I2,',  Hr',I3,  6X,  'To:',I6,A6,I2,', Hr',I3,
     *  '  Model-Time:',I9,5X,'Dif:',F7.2,' Days')
  903 FORMAT (1X,25('--'),13(A4,'--'))
  904 FORMAT (35X,'GLOBAL',A7,2X,13I6)
  905 FORMAT (A32,2F9.2,1X,13I6)
  906 FORMAT ('0AREA (10**10 M**2)',F22.1,F9.1,1X,13I6)
  907 FORMAT ('0')
      END SUBROUTINE DIAGCP


      SUBROUTINE DIAG5P
!@sum  DIAG5P PRINTS THE SPECTRAL ANALYSIS TABLES
!@auth Gary Russell
!@ver  1.0
c      USE PRTCOM, only :
      USE CONSTANT, only :
     &     grav,rgas,twopi
      USE MODEL_COM, only :
     &     im,jm,lm,fim,
     &     DT,IDACC,JHOUR,JHOUR0,JDATE,JDATE0,
     &     AMON,AMON0,JYEAR,JYEAR0,LS1,JEQ,XLABEL,istrat
      USE GEOM, only :
     &     DLON,DXYV
      USE DAGCOM, only :
     &     speca,atpe,ajk,aijk,kspeca,ktpe,nhemi,nspher,ijk_u,klayer
     &    ,JK_DPB
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(IM) :: X
      DOUBLE PRECISION, DIMENSION(KSPECA) :: SCALET,F0,FNSUM
      INTEGER, DIMENSION(KSPECA) :: MN

      DOUBLE PRECISION, DIMENSION(KTPE,NHEMI) :: FATPE

      INTEGER, PARAMETER :: IZERO=0

      INTEGER, DIMENSION(KTPE), PARAMETER ::
     &     MAPEOF=(/3,8,10,11,13,15,17,20/)

      CHARACTER*8 :: LATITD(4) = (/
     *     'SOUTHERN','NORTHERN',' EQUATOR','45 NORTH'/)
      CHARACTER*16 :: SPHERE(4)=
     *     (/'TROPOSPHERE     ','LOW STRATOSPHERE',
     *       'MID STRATOSPHERE','UPP STRATOSPHERE'/)

      INTEGER ::
     &     I,IUNITJ,IUNITW,J,J45N,
     &     K,KPAGE,KROW,KSPHER,L,
     &     M,MAPE,MTPE,N,NM,NM1

      DOUBLE PRECISION ::
     &     FACTOR,FNM

      NM=1+IM/2
      IF (IDACC(12).LT.1) IDACC(12)=1
      J45N=2.+.75*(JM-1.)
C****
C**** STANDING KINETIC ENERGY
C****
      DO 710 K=1,NSPHER
      DO 710 N=1,NM
  710 SPECA(N,1,K)=0.
      DO 770 L=1,LM
        KSPHER=KLAYER(L)
      DO 770 J=2,JM
      IF (AJK(J,L,JK_DPB).LE.1.D-20) GO TO 770
      FACTOR=FIM*DXYV(J)/AJK(J,L,JK_DPB)
      DO 769 K=IZERO,LM,LM
      DO 720 I=1,IM
  720 X(I)=AIJK(I,J,L+K,IJK_U)
      CALL FFTE (X,X)
      IF (J.EQ.JEQ) GO TO 750
      DO 730 N=1,NM
  730 SPECA(N,1,KSPHER)=SPECA(N,1,KSPHER)+X(N)*FACTOR
      IF (J.NE.J45N) GO TO 769
      DO 740 N=1,NM
  740 SPECA(N,1,KSPHER+2)=SPECA(N,1,KSPHER+2)+X(N)*FACTOR
      GO TO 769
  750 DO 760 N=1,NM
      SPECA(N,1,KSPHER+2)=SPECA(N,1,KSPHER+2)+X(N)*FACTOR
      SPECA(N,1,KSPHER)=SPECA(N,1,KSPHER)+.5*X(N)*FACTOR
  760 SPECA(N,1,KSPHER+1)=SPECA(N,1,KSPHER+1)+.5*X(N)*FACTOR
      IF (K.EQ.LM) KSPHER=KSPHER+1
  769 CONTINUE
  770 CONTINUE
C****
  600 SCALET(1)=100.D-17/(GRAV*IDACC(4)+1.D-20)
      SCALET(19)=100.D-17/(GRAV*IDACC(12))
      SCALET(20)=SCALET(19)*RGAS
      SCALET(2)=SCALET(19)*IDACC(12)/(IDACC(7)+1.D-20)
      SCALET(3)=SCALET(2)*RGAS
      SCALET(4)=100.D-12/(GRAV*DT*IDACC(6)+1.D-20)
      SCALET(5)=SCALET(4)
      SCALET(6)=SCALET(4)
      SCALET(7)=100.D-12/(GRAV*DT*(IDACC(7)+1.D-20))
      SCALET(8)=SCALET(7)*RGAS
      SCALET(9)=100.D-12/(GRAV*DT*(IDACC(8)+1.D-20))
      SCALET(10)=SCALET(9)*RGAS
      SCALET(11)=SCALET(10)
      SCALET(12)=SCALET(9)
      SCALET(13)=SCALET(10)
      SCALET(14)=100.D-12/(GRAV*DT*(IDACC(10)+1.D-20))
      SCALET(15)=SCALET(14)*RGAS
      SCALET(16)=100.D-12/(GRAV*DT*(.5*IDACC(9)+1.D-20))
      SCALET(17)=SCALET(16)*RGAS
      SCALET(18)=100.D-17/(GRAV*IDACC(4)+1.D-20)
      DO 605 K=1,KSPECA
  605 SCALET(K)=(TWOPI/(DLON*FIM))*SCALET(K)
      IUNITJ=17
      IUNITW=12
      DO 690 KPAGE=1,4   ! one for each loc. SH/NH/EQ/45N
C**** WRITE HEADINGS
      WRITE (6,901) XLABEL
      WRITE (6,902) JYEAR0,AMON0,JDATE0,JHOUR0,JYEAR,AMON,JDATE,JHOUR,
     *  IUNITJ,IUNITW
      DO 670 KROW=1,2+ISTRAT
      IF (JM.GE.25.AND.KROW.EQ.2) WRITE (6,901)
      WRITE (6,903) LATITD(KPAGE),SPHERE(KROW)
      KSPHER=4*(KROW-1)+KPAGE
C**** WRITE KINETIC AND AVAILABLE POTENTIAL ENERGY BY WAVE NUMBER
      DO 610 M=1,KSPECA
      F0(M)=SPECA(1,M,KSPHER)*SCALET(M)
      MN(M)=NINT(F0(M))
  610 FNSUM(M)=0.
      WRITE (6,904) MN
      DO 630 N=2,NM
      KSPHER=4*(KROW-1)+KPAGE
      DO 620 M=1,KSPECA
      FNM=SPECA(N,M,KSPHER)*SCALET(M)
      MN(M)=NINT(FNM)
  620 FNSUM(M)=FNSUM(M)+FNM
      NM1=N-1
  630 WRITE (6,905) NM1,MN
      DO 640 M=1,KSPECA
  640 MN(M)=NINT(FNSUM(M))
      WRITE (6,906) MN
      DO 650 M=1,KSPECA
  650 MN(M)=NINT(FNSUM(M)+F0(M))
      WRITE (6,907) MN
  670 CONTINUE
      IF (KPAGE.GE.3) GO TO 690
C**** WRITE TOTAL POTENTIAL ENERGY
      DO 680 MTPE=1,KTPE
      MAPE=MAPEOF(MTPE)
         FATPE(MTPE,KPAGE)=ATPE(MTPE,KPAGE)*SCALET(MAPE)/RGAS
  680 MN(MTPE)=NINT(FATPE(MTPE,KPAGE))
      WRITE (6,909) (MN(MTPE),MTPE=1,8)
      IF (KPAGE.NE.2) GO TO 690
      DO 685 M=1,KSPECA
  685 SCALET(M)=SCALET(M)*10.
      IUNITJ=16
      IUNITW=11
  690 CONTINUE
      RETURN
C****
  901 FORMAT ('1',A)
  902 FORMAT ('0**  Spectral Analysis **      From:',
     *  I6,A6,I2,',  Hr',I3,  6X,  'To:',I6,A6,I2,', Hr',I3,
     *  '       UNITS 10**',I2,' JOULES AND 10**',I2,' WATTS')
  903 FORMAT ('0',50X,A8,1X,A16/
     *  13X,'MEAN',19X,'DYNAMICS',25X,'SOURCES',16X,'FILTER',8X,
     *     'DAILY',4X,'PR SURF',5X,'LAST'/
     *'   N    SKE   KE   APE    KADV  KCOR   P-K  KDYN  PDYN   ',
     *     'KCNDS PCNDS   PRAD KSURF PSURF   KFIL  PFIL   KGMP  PGMP',
     *     '    KE',6X,'KE   APE')
  904 FORMAT ( '0  0',I7,I5,I6,I8,4I6,I8,I6,I7,2I6,I7,I6,I7,2I6,I8,I6/)
  905 FORMAT (     I4,I7,I5,I6,I8,4I6,I8,I6,I7,2I6,I7,I6,I7,2I6,I8,I6)
  906 FORMAT (' EDDY',I6,I5,I6,I8,4I6,I8,I6,I7,2I6,I7,I6,I7,2I6,I8,I6)
  907 FORMAT ('0TOTL',I6,I5,I6,I8,4I6,I8,I6,I7,2I6,I7,I6,I7,2I6,I8,I6)
  909 FORMAT (/'0TPE',I18,I32,I14,I7,I12,2I13,I20)
      END SUBROUTINE DIAG5P


      SUBROUTINE DIAGDD
!@sum  DIAGDD prints out diurnal cycle diagnostics
!@auth G. Russell
!@ver  1.0
c      USE PRTCOM, only :
      USE MODEL_COM, only :
     &     idacc,JDATE,JDATE0,AMON,AMON0,JYEAR,JYEAR0,XLABEL
      USE DAGCOM, only :   kdiag,
     &     adiurn,ijdd,namdd,ndiuvar,hr_in_day,scale_dd,lname_dd,name_dd
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(HR_IN_DAY+1) :: XHOUR
      INTEGER, DIMENSION(HR_IN_DAY+1) :: MHOUR
      DOUBLE PRECISION :: AVE,AVED,AVEN,BYIDAC
      INTEGER :: I,IH,IREGF,IREGL,IS,K,KQ,KR,NDAYS
C****
      NDAYS=IDACC(9)/2
      IF (NDAYS.LE.0) RETURN
      BYIDAC=1./NDAYS
C****
      IREGF=1
      IREGL=4
      IF (KDIAG(6).GT.0) IREGL=4-KDIAG(6)
      IF (KDIAG(6).LT.0.AND.KDIAG(6).GT.-5) IREGF=-KDIAG(6)
      IF (KDIAG(6).LT.0) IREGL=IREGF
      DO KR=IREGF,IREGL
        WRITE (6,901) XLABEL(1:105),JDATE0,AMON0,JYEAR0,JDATE,AMON,JYEAR
        WRITE (6,903) NAMDD(KR),IJDD(1,KR),IJDD(2,KR),(I,I=1,HR_IN_DAY)
        DO KQ=1,NDIUVAR
          IF (MOD(KQ-1,5).eq.0) WRITE(6,*)
          IF (LNAME_DD(KQ).eq."unused") CYCLE
          SELECT CASE (NAME_DD(KQ))
          CASE DEFAULT
C**** NORMAL QUANTITIES
            AVE=0.
            DO IH=1,HR_IN_DAY
              AVE=AVE+ADIURN(IH,KQ,KR)
              XHOUR(IH)=ADIURN(IH,KQ,KR)*SCALE_DD(KQ)*BYIDAC
            END DO
            XHOUR(HR_IN_DAY+1)=AVE/FLOAT(HR_IN_DAY)*SCALE_DD(KQ)*BYIDAC
C**** RATIO OF TWO QUANTITIES
          CASE ('LDC')
            AVEN=0.
            AVED=0.
            DO IH=1,HR_IN_DAY
              AVEN=AVEN+ADIURN(IH,KQ,KR)
              AVED=AVED+ADIURN(IH,KQ-1,KR)
              XHOUR(IH)=ADIURN(IH,KQ,KR)*SCALE_DD(KQ)/
     *             (ADIURN(IH,KQ-1,KR)+1D-20)
            END DO
            XHOUR(HR_IN_DAY+1)=AVEN*SCALE_DD(KQ)/(AVED+1D-20)
          END SELECT
          DO IS=1,HR_IN_DAY+1
CB        FHOUR(IS,KQ,KR)=XHOUR(IS)
            MHOUR(IS)=NINT(XHOUR(IS))
          END DO
          WRITE (6,904) LNAME_DD(KQ),MHOUR
        END DO
      END DO
      RETURN
C****
  901 FORMAT ('1',A,I3,1X,A3,I5,' - ',I3,1X,A3,I5)
  903 FORMAT ('0',A4,I2,',',I2,' ',I2,23I5,'  AVE')
  904 FORMAT (A8,25I5)
      END SUBROUTINE DIAGDD


      SUBROUTINE DIAG4
!@sum  DIAG4 prints out a time history of the energy diagnostics
!@auth G. Russell
!@ver  1.0
c      USE PRTCOM, only :
      USE CONSTANT, only :
     &     grav,rgas,twopi,bygrav
      USE MODEL_COM, only :
     &     im,jm,lm,fim,
     &     IDACC,JHOUR,JHOUR0,JDATE,JDATE0,AMON,AMON0,
     &     JYEAR,JYEAR0,NDA4,NDAY,Itime0,XLABEL,istrat
      USE GEOM, only :
     &     DLON
      USE DAGCOM, only :
     &     energy,ned,nehist,hist_days
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(2) :: FAC
      DOUBLE PRECISION, DIMENSION(NED) :: SCALET
      DOUBLE PRECISION, DIMENSION(2*NED) :: SUME
      INTEGER, DIMENSION(2*NED) :: IK
      DOUBLE PRECISION, DIMENSION(NEHIST,HIST_DAYS+1) :: EHIST

      INTEGER ::
     &     I,IDACC5,ItimeX,IDAYX,IDAYXM,K,K0,KS,KN,KSPHER
      DOUBLE PRECISION ::
     &     TOFDYX

      IDACC5=IDACC(5)
      IF (IDACC5.LE.0) RETURN
      IF (IDACC(12).LT.1) IDACC(12)=1
      SCALET(1)=100.D-18*BYGRAV
      SCALET(2)=SCALET(1)
      SCALET(3)=SCALET(1)
      SCALET(4)=SCALET(1)
      SCALET(5)=.5*SCALET(1)
      SCALET(6)=SCALET(5)
      SCALET(7)=SCALET(1)*RGAS
      SCALET(8)=SCALET(7)
      SCALET(9)=SCALET(7)
      SCALET(10)=SCALET(7)
      DO K=1,NED
        SCALET(K)=(TWOPI/(DLON*FIM))*SCALET(K)/IDACC(12)
      END DO
C****
      DO K0=1,MIN(1+ISTRAT,2)
        WRITE (6,901) XLABEL
        IF (K0.eq.1) THEN
          FAC(1) = 1.
          FAC(2) = 10.
          WRITE (6,902) JYEAR0,AMON0,JDATE0,JHOUR0,JYEAR,AMON,JDATE
     *         ,JHOUR
          WRITE (6,903)
        ELSE
          FAC(1) = 10.  !100.  do they correspond to titles?
          FAC(2) = 10.  !1000.
          WRITE (6,906) JYEAR0,AMON0,JDATE0,JHOUR0,JYEAR,AMON,JDATE
     *         ,JHOUR
          WRITE (6,907)
        END IF
        SUME(:)=0.
        DO I=1,IDACC5
          ItimeX=Itime0+I*NDA4-1
          IDAYX=1+ItimeX/NDAY
          IDAYXM=MOD(IDAYX,100000)
          TOFDYX=MOD(ItimeX,NDAY)*24./NDAY
          DO KSPHER=1,2
            DO K=1,NED
              KS=K+(KSPHER-1)*NED
              KN=KS+(K0-1)*NED
              IF (KN.le.NEHIST) THEN
                EHIST(KN,I)=ENERGY(KN,I)*SCALET(K)*FAC(KSPHER)
                IK(KS)=EHIST(KN,I)+.5
                SUME(KS)=SUME(KS)+ENERGY(KN,I)
              ELSE
                IK(KS)=-999
              END IF
            END DO
          END DO
          WRITE (6,904) IDAYXM,TOFDYX,IK
        END DO
        DO KSPHER=1,2
          DO K=1,NED
            KS=K+(KSPHER-1)*NED
            KN=KS+(K0-1)*NED
            IF (KN.le.NEHIST) THEN
              EHIST(KN,HIST_DAYS+1)=SUME(KS)*SCALET(K)*FAC(KSPHER)
     *             /IDACC5
              IK(KS)=EHIST(KN,HIST_DAYS+1)+.5
            ELSE
              IK(KS)=-999
            END IF
          END DO
        END DO
        WRITE (6,905) IK
        IF (K0.eq.1) CALL KEYD4 (IK)
      END DO
      RETURN
C****
  901 FORMAT ('1',A)
  902 FORMAT ('0** ENERGY HISTORY **      From:',I6,A6,I2,',  Hr',I3,
     *  6X,'To:',I6,A6,I2,', Hr',I3,
     *  '    UNITS OF 10**18 JOULES')
  903 FORMAT ('0',15X,21('-'),' TROPOSPHERE ',22('-'),5X,21('-'),
     *  '  LOW STRAT. * 10 ',17('-')/8X,2(11X,'ZKE',8X,'EKE',7X,
     *     'SEKE',9X,
     * 'ZPE',10X,'EPE')/3X,'DAY  HOUR     SH   NH    SH   NH     1    2
     *    SH    NH     SH    NH      SH   NH    SH   NH    SH   NH     S
     *H    NH     SH    NH'/1X,132('='))
  904 FORMAT (I6,F6.1,1X,3(I6,I5),2(I7,I6),2X,3(I6,I5),2(I7,I6))
  905 FORMAT (1X,132('=')/8X,'MEAN ',3(I6,I5),2(I7,I6),2X,3(I6,I5),
     *  2(I7,I6))
  906 FORMAT ('0** ENERGY HISTORY **      From:',I6,A6,I2,',  Hr',I3,
     *  6X,'To:',I6,A6,I2,', Hr',I3,
     *  '    UNITS OF 10**17 JOULES')
  907 FORMAT ('0',15X,19('-'),' MID STRATOSPHERE ',19('-'),5X,18('-'),
     *  ' HIGH STRAT. * 10  ',19('-')/8X,2(11X,'ZKE',8X,'EKE',7X,
     *   'NHKE',9X,
     *  'ZPE',10X,'EPE')/3X,'DAY  HOUR     SH   NH    SH   NH    1    2
     *    SH    NH     SH    NH      SH   NH    SH   NH    1    2      S
     *H    NH     SH    NH'/1X,132('='))
  920 FORMAT (1X)
      END SUBROUTINE DIAG4


      subroutine KEYVSUMS (QUANT,GSUM,HSUM,ASUM,SUMFAC)
      USE MODEL_COM, only : jm
      implicit none
!@var quant string designating the quantity for which to save keynrs
      CHARACTER(LEN=30) :: QUANT
      DOUBLE PRECISION, DIMENSION(JM) :: ASUM
      DOUBLE PRECISION, DIMENSION(2) :: HSUM
      DOUBLE PRECISION :: GSUM,SUMFAC
      if(quant.eq.'temp') CALL KEYJKT (GSUM,ASUM)
      if(quant.eq.'eddy_ke') CALL KEYJKE (18,HSUM,ASUM)
      if(quant.eq.'tot_ke') CALL KEYJKE (19,HSUM,ASUM)
      if(quant.eq.'nt_dse_stand_eddy') CALL KEYJKN (22,ASUM,SUMFAC)
      if(quant.eq.'nt_dse_eddy') CALL KEYJKN (23,ASUM,SUMFAC)
      if(quant.eq.'tot_nt_dse') CALL KEYJKN (24,ASUM,SUMFAC)
      if(quant.eq.'nt_lh_e') CALL KEYJKN (25,ASUM,SUMFAC)
      if(quant.eq.'tot_nt_lh') CALL KEYJKN (26,ASUM,SUMFAC)
      if(quant.eq.'nt_se_eddy') CALL KEYJKN (27,ASUM,SUMFAC)
      if(quant.eq.'tot_nt_se') CALL KEYJKN (28,ASUM,SUMFAC)
      if(quant.eq.'tot_nt_ke') CALL KEYJKN (30,ASUM,SUMFAC)
      if(quant.eq.'nt_am_stand_eddy') CALL KEYJKN (31,ASUM,SUMFAC)
      if(quant.eq.'nt_am_eddy') CALL KEYJKN (32,ASUM,SUMFAC)
      if(quant.eq.'tot_nt_am') CALL KEYJKN (33,ASUM,SUMFAC)
      RETURN
      end subroutine keyvsums


      subroutine keynrl(quant,l,flat)
      USE MODEL_COM, only : jm
      implicit none
      integer :: l
      DOUBLE PRECISION, DIMENSION(JM) :: FLAT
!@var quant string designating the quantity for which to save keynrs
      CHARACTER(LEN=30) :: QUANT
      if(quant.eq.'u') CALL KEYJKJ (L,FLAT)
      if(quant.eq.'psi_cp') CALL KEYJLS (L,FLAT)
      return
      end subroutine keynrl


      SUBROUTINE DIAGKS
C****
C**** THIS ROUTINE PRODUCES A SUMMARY OF KEY NUMBERS CALCULATED IN
C**** OTHER DIAGNOSTIC ROUTINES
C****
C**** CONTENTS OF KEYNR
C****
C   K   N
C****
C***1  1 MONTH
C***2  2 TOTAL CLOUD COVER (PERCENT)
C****  3 SNOW COVER--NORTHERN HEMSIPHERE (PERCENT)
C****  4 ICE COVER--NORTHERN HEMISPHERE (PERCENT)
C****  5 PLANETARY ALBEDO (PERCENT)
C****  6 SOLAR RADIATION ABSORBED BY ATMOSPHERE (WT/M**2)
C****  7 SOLAR RADIATION ABSORBED BY PLANET (WT/M**2)
C****  8 NET HEAT AT GROUND (WT/M**2)
C****  8 ANGULAR MOMENTUM PER UNIT AREA (10**10 J*SEC/M**2)
C****  9 EVAPORATION (.1 MM/DAY)
C****  9 PRECIPITATION (.1 MM/DAY)
C**** 10 SENSIBLE HEAT FLUX INTO GROUND (ABS.VALUE)
C**** 11 LATENT HEAT FLUX INTO GROUND (ABS.VALUE)
C**** 12 MEAN GROUND TEMPERATURE (DEGREES K)
C**** 13 MEAN GLOBAL ATMOSPHERIC TEMPERATURE (DEGREES K)
C**** 14 MERID. TEMPERATURE GRADIENT (N.HEMISPHERE)
C**** 15 MERID. TEMPERATURE GRADIENT (S.HEMISPHERE)
C**** 16 MEAN TROPOSPHERIC EKE-NORTHERN HEMISPHERE
C**** 17 MEAN TROPOSPHERIC EKE-SOUTHERNN HEMISPHERE
C**** 18 MEAN TROPOSPHERIC ZKE-NORTHERN HEMISPHERE
C**** 19 MEAN TROPOSPHERIC ZKE-SOUTHERN HEMISPHERE
C**** 20 MEAN TROPOSPHERIC EPE-NORTHERN HEMISPHERE
C**** 21 MEAN TROPOSPHERIC ZPE-NORTHERN HEMISPHERE
C**** 22 MEAN EDDY KINETIC ENERGY AT EQUATOR
C**** 23 MAX. MEAN EDDY KINETIC ENERGY IN MID NORTH LATITUDES
C**** 24 MAX. ZONAL WIND (U COMPONENT) IN TROPOSPHERE (NH), M/SEC
C**** 25 LATITUDE CORRESPONDING TO 24
C**** 26 MAX. ZONAL WIND (U COMPONENT) IN TROPOSPHERE (SH), M/SEC
C**** 27 LATITUDE CORRESPONDING TO 26
C**** 28-30: 29 IS LARGEST VALUE OF STREAM FUNCTION, POSITIVE OR
C****    NEGATIVE; 28 AND 30 ARE THE MAGNITUDES OF THE LARGEST VALUES OF
C****    OPPOSITE SIGN TO THE NORTH AND SOUTH RESPECTIVELY
C***3 31 SNOW AND ICE COVERAGE OF GLOBE (PERCENT)
C***4 32 SNOW AND ICE COVERAGE OF NORTHERN HEMISPHERE (PERCENT)
C****   33-39 REFER TO NORTHERN HEMISPHERE ONLY
C**** 33 MAX.NORTHWARD TRANS. OF DRY STATIC ENERGY BY STANDING EDDIES
C**** 34 MAX.NORTHWARD TRANS. OF DRY STATIC ENERGY BY EDDIES
C**** 35 MAX. TOTAL NORTH. TRANS. OF DRY STATIC ENERGY
C**** 36 MAX.NORTHWARD TRANS. OF STATIC ENERGY BY EDDIES
C**** 37 MAX.TOTAL NORTH. TRANS. OF STATIC ENERGY
C**** 38 LATITUDE CORRESPONDING TO 37
C**** 39 MAX. NORTH. TRANS. OF ANGULAR MOMENTUM BY STANDING EDDIES
C**** 40 MAX. NORTH. TRANS. OF ANGULAR MOMENTUM BY EDDIES
C**** 41 MAX. TOTAL NORTH. TRANS. OF ANGULAR MOMENTUM
C**** 42 LATITUDE CORRESPONDING TO 41
C****
c      USE PRTCOM, only :
      USE CONSTANT, only :
     &     twopi
      USE MODEL_COM, only :
     &     jm,lm,jeq, JHOUR,JHOUR0,
     &     JDATE,JDATE0,JMON,JMON0,AMON,AMON0,JYEAR,JYEAR0,
     &     Itime,ItimeI,Itime0,PTOP,SIG,XLABEL,
     &     PSFMPT,AMONTH,nday
      USE GEOM, only :
     &     DLAT,DXYP,LAT_DG
      USE DAGCOM, only :
     &     keyct,keynr,nehist,nkeynr
      USE PARAM
      IMPLICIT NONE
      SAVE
      INTEGER*4 :: NDEX(42) = (/
     *     1,2,31,32,3,   4,5,6,7,8,   9,10,11,12,13,
     *             14,15,16,17,18,  19,20,21,22,23,  24,25,26,27,28,
     *             29,30,33,34,35,  36,37,38,39,40,  41,42/)

      DOUBLE PRECISION, DIMENSION(JM) :: FLAT
      DOUBLE PRECISION, DIMENSION(JM,LM) :: FKEY
      DOUBLE PRECISION, DIMENSION(JM) :: ASUM
      DOUBLE PRECISION, DIMENSION(2) :: HSUM
      INTEGER, DIMENSION(NEHIST) :: IK

      INTEGER ::
     &     I,I35,I70,IEND,ISIGN,
     &     J,J60,JMAX,JEQP1,JNDEX,JSTART,
     &     K,KEYMAX,KNDEX,
     &     L,LL,LMAX,LNLM,LNM,LSLM,LSM,
     &     N,NT,NTDIF

      DOUBLE PRECISION ::
     &     A,BIG,CEPT,CHECK,
     &     FGLOB,FNH,GSUM,HN,HS,PISG,PISN,
     &     SAVE,SUMFAC,DAYS,TEQ,TNOR,TSOU,
     &     UNLM,UNM,USLM,USM,X60

C****
C**** ENTRIES CALLED FROM DIAGJ
C****
      ENTRY KEYDJ (N,FGLOB,FNH)
      GO TO (                100,100,100,110,100, 100,100,100,100,115,
     *  100,100,120,125,100, 100,100,130,100,135, 100,100,100,100,100,
     *  100,100,100,100,140, 145,100,100,100,100, 100,100,100,100,100,
     *  100,100,100,150,100, 100,100,100,100,100, 100,100,100,100,100,
     *  100,100,100,155),N
  100 RETURN
  110 KEYNR(6,KEYCT)=NINT(FGLOB)
      RETURN
  115 KEYNR(7,KEYCT)=NINT(FGLOB)
      RETURN
  120 KEYNR(10,KEYCT)=NINT(-FGLOB)
      RETURN
  125 KEYNR(11,KEYCT)=NINT(-FGLOB)
      RETURN
  130 KEYNR(12,KEYCT)=NINT(.1*FGLOB)
      RETURN
  135 KEYNR(9,KEYCT)=NINT(10.*FGLOB)
      RETURN
  140 KEYNR(4,KEYCT)=NINT(FNH)
      RETURN
  145 KEYNR(3,KEYCT)=NINT(FNH)
      RETURN
  150 KEYNR(8,KEYCT)=NINT(FGLOB)
      RETURN
  155 KEYNR(2,KEYCT)=NINT(FGLOB)
      RETURN
C****
      ENTRY KEYDJA (FGLOB)
      KEYNR(5,KEYCT)=NINT(10.*FGLOB)
      RETURN
C****
C**** ENTRIES CALLED FROM DIAGJL VIA JLMAP OR FROM DIAGJK VIA JKMAP
C****
      ENTRY KEYJKT (GSUM,ASUM)
C**** TEMPERATURES
C      JEQ=2.+.5*(JM-1.)
      TEQ=.5*(ASUM(JEQ-1)+ASUM(JEQ))
      X60=TWOPI/(12.*DLAT)
      J60=.5+X60
      A=DXYP(J60+1)*(X60+.5-J60)
      TSOU=ASUM(J60+1)*A
      TNOR=ASUM(JM-J60)*A
      DO 210 J=1,J60
      A=A+DXYP(J)
      TSOU=TSOU+ASUM(J)*DXYP(J)
  210 TNOR=TNOR+ASUM(JM+1-J)*DXYP(J)
      KEYNR(14,KEYCT)=NINT(TEQ-TNOR/A)
      KEYNR(15,KEYCT)=NINT(TEQ-TSOU/A)
      KEYNR(13,KEYCT)=NINT(.1*GSUM)
      RETURN
C****
      ENTRY KEYJKJ (L,FLAT)
C**** JET STREAMS
      IF (L.LT.LM) GO TO 220
      DO 216 LL=1,LM
      IF (PSFMPT*SIG(LL)+PTOP.LT.200.) GO TO 218
  216 CONTINUE
  218 LMAX=LL-1
  220 IF (L.GT.LMAX) RETURN
      USLM=-999999.
      DO 222 J=3,JEQ
      IF (FLAT(J).LT.USLM) GO TO 222
      USLM=FLAT(J)
      JMAX=J
  222 CONTINUE
      CEPT=.5*(FLAT(JMAX-1)-FLAT(JMAX+1))/
     *  (FLAT(JMAX-1)-2.*FLAT(JMAX)+FLAT(JMAX+1))
      LSLM=INT((JMAX-1.5+CEPT)*DLAT*360/TWOPI+.5)-90
      UNLM=-999999.
      DO 224 J=JEQ,JM-1
      IF (FLAT(J).LT.UNLM) GO TO 224
      UNLM=FLAT(J)
      JMAX=J
  224 CONTINUE
      CEPT=.5*(FLAT(JMAX-1)-FLAT(JMAX+1))/
     *  (FLAT(JMAX-1)-2.*FLAT(JMAX)+FLAT(JMAX+1))
      LNLM=INT((JMAX-1.5+CEPT)*DLAT*360/TWOPI+.5)-90
      IF (L.LT.LMAX) GO TO 226
      USM=USLM
      LSM=LSLM
      UNM=UNLM
      LNM=LNLM
      RETURN
  226 IF (USLM.LT.USM) GO TO 228
      USM=USLM
      LSM=LSLM
  228 IF (UNLM.LT.UNM) GO TO 230
      UNM=UNLM
      LNM=LNLM
  230 IF (L.NE.1) RETURN
      KEYNR(24,KEYCT)=.1*UNM+.5
      KEYNR(25,KEYCT)=LNM
      KEYNR(26,KEYCT)=.1*USM+.5
      KEYNR(27,KEYCT)=-LSM
      RETURN
C****
      ENTRY KEYJLS (L,FLAT)
C**** STREAM FUNCTION
      DO 290 J=2,JM
  290 FKEY(J,L)=FLAT(J)
      IF (L.NE.1) RETURN
  300 SAVE=0.
      HS=0.
      HN=0.
      DO 310 K=1,LM
      DO 310 I=2,JM
      CHECK=ABS(FKEY(I,K))
      IF (CHECK.LT.SAVE) GO TO 310
      SAVE=CHECK
      JNDEX=I
      KNDEX=K
  310 CONTINUE
      SAVE=FKEY(JNDEX,KNDEX)
      ISIGN=1
      IF (SAVE.GT.0.0) ISIGN=-1
      IF (JNDEX.LT.4) GO TO 325
      IEND=JNDEX-1
      DO 320 K=1,LM
      DO 320 I=2,IEND
      CHECK=FKEY(I,K)*ISIGN
  320 IF (CHECK.GT.HS)HS=CHECK
  325 CONTINUE
      IF (JNDEX.GT.(JM-2))GO TO 335
      JSTART=JNDEX+1
      DO 330 K=1,LM
      DO 330 I=JSTART,JM
      CHECK=FKEY(I,K)*ISIGN
  330 IF (CHECK.GT.HN)HN=CHECK
  335 CONTINUE
      KEYNR(28,KEYCT)=ABS(HN)+0.5
      KEYNR(29,KEYCT)=NINT(SAVE)
      KEYNR(30,KEYCT)=ABS(HS)+0.5
      RETURN
C****
      ENTRY KEYJKE (NT,HSUM,ASUM)
C**** EDDY AND ZONAL KINETIC ENERGY
      IF (NT.EQ.19) GO TO 450
      KEYNR(16,KEYCT)=NINT(HSUM(2))
      KEYNR(17,KEYCT)=NINT(HSUM(1))
      KEYNR(18,KEYCT)=KEYNR(18,KEYCT)-NINT(HSUM(2))
      KEYNR(19,KEYCT)=KEYNR(19,KEYCT)-NINT(HSUM(1))
      KEYNR(22,KEYCT)=NINT(ASUM(JEQ))
      BIG=-99999.
      I35=2.+(JM-1.)*125./180.
      I70=2.+(JM-1.)*160./180.
      DO 440 I=I35,I70
      IF (ASUM(I).LT.BIG) GO TO 440
      BIG=ASUM(I)
  440 CONTINUE
      KEYNR(23,KEYCT)=NINT(BIG)
      RETURN
  450 KEYNR(18,KEYCT)=KEYNR(18,KEYCT)+NINT(HSUM(2))
      KEYNR(19,KEYCT)=KEYNR(19,KEYCT)+NINT(HSUM(1))
      RETURN
C****
      ENTRY KEYJKN (NT,ASUM,SUMFAC)
C**** NORTHWARD TRANSPORTS
  500 BIG=-99999.
      DO 510 I=JEQ+1,JM
      IF (ASUM(I).LT.BIG) GO TO 510
      BIG=ASUM(I)
      JNDEX=I
  510 CONTINUE
      BIG=BIG*SUMFAC
      NTDIF=NT-21
      GO TO (392,392,392,390,390,396,394,390,390,400,400,398),NTDIF
  390 CONTINUE
  392 KEYNR(NT+11,KEYCT)=NINT(BIG)
      RETURN
  394 KEYNR(38,KEYCT)=NINT(LAT_DG(JNDEX,2))
  396 KEYNR(NT+9,KEYCT)=NINT(BIG)
      RETURN
  398 KEYNR(42,KEYCT)=NINT(LAT_DG(JNDEX,2))
  400 KEYNR(NT+8,KEYCT)=NINT(BIG)
      RETURN
C****
C**** ENTRY CALLED FROM DIAGIJ
C****
      ENTRY KEYIJ(PISG,PISN)
      KEYNR(31,KEYCT)=NINT(PISG)
      KEYNR(32,KEYCT)=NINT(PISN)
      RETURN
C****
C**** ENTRY CALLED FROM DIAG4
C****
      ENTRY KEYD4 (IK)
      KEYNR(20,KEYCT)=(IK(10)+IK(20)+5)/10
      KEYNR(21,KEYCT)=(IK(8)+IK(18)+5)/10
      RETURN
C****
      ENTRY DIAGKN
C**** PRINTS THE TABLE OF KEY NUMBERS
C****
      DAYS=(Itime-Itime0)/DFLOAT(nday)
      KEYNR(1,KEYCT)=JMON0
      IF (Itime.eq.ItimeI+1) KEYNR(1,KEYCT)=0
      IF (KEYCT.GE.2.AND.KEYNR(1,KEYCT-1).EQ.JMON0) KEYCT=KEYCT-1
      WRITE(6,901) XLABEL
      WRITE(6,910) JYEAR0,AMON0,JDATE0,JHOUR0,
     *  JYEAR,AMON,JDATE,JHOUR,ITIME,DAYS
      WRITE(6,902)
      DO 810 K=1,KEYCT
      IF (KEYNR(1,K).EQ.1) WRITE (6,905)
  810 WRITE(6,905) AMONTH(KEYNR(NDEX(1),K)),(KEYNR(NDEX(I),K),I=2,42)
      WRITE (6,915)
CB       DO 815 I=1,42
CB815    FKEYDS(I)=KEYNR(I,KEYCT)
      KEYCT=KEYCT+1
      KEYMAX=49
      IF (KEYNR(1,1).NE.0) KEYMAX=48
      IF (KEYCT.LE.KEYMAX) RETURN
C**** ROLL UP KEY NUMBERS 1 YEAR AT A TIME
      DO 820 K=1,36
      DO 820 I=1,NKEYNR
  820 KEYNR(I,K)=KEYNR(I,K+KEYMAX-36)
      DO 880 K=37,50
      DO 880 I=1,NKEYNR
  880 KEYNR(I,K)=0
      KEYCT=37
      RETURN
  901 FORMAT('1',A)
  902 FORMAT ('0',7X,'SN+IC NH NH AL AB NT NT PR        T   T-OF-ATM  EK
     *E   ZKE           EKE   JET-STREAMS STREAM-FN NOR-TRAN NOR-TRAN NO
     *RTH-TRANS'/
     *         5X,'CL GL    SN OI BE BY RD HT EC SN LAT OF  GL  GRAD ---
     *-- ----- EPE ZPE ------ NORTH SOUTH --------- DRY-STAT STAT-ENR AN
     *G MOMENTM'/
     *         5X,'CV OB NH CV CV DO AT P0 Z0 IP HT  HT GD  OB NH SH NH
     *SH NH SH  NH  NH EQ  ML VL LT VL LT NH MAX SH SE ED TL ED TL LT SE
     * ED TL LT'/)
  905 FORMAT (1X,A3,4I3,I2,I4,5I3,I4,I3,I4,6I3,2I4,I3,I4,5I3,I4,11I3)
  910 FORMAT ('0',15X,'From:',I6,A6,I2,',  Hr',I3,
     *  6X,'To:',I6,A6,I2,', Hr',I3,
     *  '  Model-Time:',I9,5X,'Dif:',F7.2,' Days')
  915 FORMAT('0')
      END SUBROUTINE DIAGKS


      SUBROUTINE IJK_TITLES
!@sum  IJK_TITLES extra titles for composite ijk output
!@auth G. Schmidt/M. Kelley
!@ver  1.0
      USE CONSTANT, only : bygrav
      USE DAGCOM, only : kaijk,kaijkx,UNITS_IJK,NAME_IJK,LNAME_IJK,
     *   units_ijk,name_ijk,lname_ijk,scale_ijk  !  dagcom
      IMPLICIT NONE
      INTEGER :: K
c
      k = kaijk
c
      k = k + 1
      name_ijk(k) = 'HEIGHT'
      lname_ijk(k) = 'HEIGHT'
      units_ijk(k) = 'm'
      scale_ijk(k) = 0.25*BYGRAV

      return
      END SUBROUTINE IJK_TITLES


      SUBROUTINE IJKMAP (iu_Iij)
!@sum  IJKMAP output 3-D constant pressure output fields
!@auth G. Schmidt
!@ver  1.0
C**** Note that since all IJK diags are weighted w.r.t pressure, all
C**** diagnostics must be divided by the accumulated pressure
C**** All titles/names etc. implicitly assume that this will be done.
      USE CONSTANT, only :
     &     grav,sha,undef
      USE MODEL_COM, only :
     &     im,jm,lm,
     &     PTOP,SIG,PSFMPT,XLABEL,LRUNID
      USE DAGCOM, only : kdiag,
     &     aijk,acc_period,ijk_u,ijk_v,ijk_t,ijk_q,ijk_dp,ijk_dse
     *     ,scale_ijk,off_ijk,name_ijk,lname_ijk,units_ijk,kaijk,kaijkx
      use filemanager
      IMPLICIT NONE

      CHARACTER XLB*24,TITLEX*56
      CHARACTER*80 TITLEL(LM)
      REAL*8 SMAP(IM,JM-1,LM),SMAPJK(JM-1,LM)
      REAL*8 flat,press,dp
      CHARACTER*4 CPRESS(LM)
      INTEGER i,j,l,kxlb,ni,kcomp,k,iu_Iij
      logical, dimension (kaijkx) :: Qk

C**** OPEN PLOTTABLE OUTPUT FILE
      call open_ijk(trim(acc_period)//'.ijk'//XLABEL(1:LRUNID))
C****
C**** INITIALIZE CERTAIN QUANTITIES
C****
      call ijk_titles

      Qk = .true.
      do k=1,kaijkx
        if (lname_ijk(k).eq.'unused') Qk(k) = .false.
      end do
      if (kdiag(3).eq.8) then
         write (iu_Iij,'(a)') 'list of 3-d fields'
         do k=1,kaijkx
           if (lname_ijk(k).ne.'unused')
     *        write (iu_Iij,'(i3,1x,a)') k,lname_ijk(k)
         end do
         stop 'created file Iij'
      else if (kdiag(3).gt.1) then
         Qk = .false.
   10    read (iu_Iij,'(i3)',end=20) k
         Qk(k) = .true.
         go to 10
   20    continue
         call closeunit(iu_Iij)
      end if

      SMAP(:,:,:) = UNDEF
      KXLB = INDEX(XLABEL(1:11),'(')-1
      IF(KXLB.le.0) KXLB = 10
      XLB = ' '
      XLB(1:13)=acc_period(1:3)//' '//acc_period(4:12)
      XLB(15:14+KXLB) = XLABEL(1:KXLB)
C****
C**** Complete 3D-field titles
C****
      DO L=1,LM
        PRESS = PTOP+PSFMPT*SIG(L)
        IF (PRESS.GE.1.) THEN
          WRITE(CPRESS(L),'(I4)') NINT(PRESS)
        ELSE
          WRITE(CPRESS(L),'(F4.3)') PRESS
        END IF
      END DO

C**** Select fields
      DO K=1,KAIJKx
        if (.not.Qk(k).or.k.eq.ijk_dp) cycle
        TITLEX = lname_ijk(k)(1:17)//"   at        mb ("//
     *       trim(units_ijk(k))//", UV grid)"
        IF (K.le.kaijk) THEN     !  simple cases
          DO L=1,LM
            DO J=2,JM
              NI = 0
              FLAT = 0.
              SMAPJK(J-1,L) = UNDEF
              DO I=1,IM
                DP=AIJK(I,J,L,IJK_DP)
                SMAP(I,J-1,L) = UNDEF
                IF(DP.GT.0.) THEN
                  SMAP(I,J-1,L)=SCALE_IJK(K)*AIJK(I,J,L,K)/DP+OFF_IJK(K)
                  FLAT = FLAT+SMAP(I,J-1,L)
                  NI = NI+1
                END IF
              END DO
              IF (NI.GT.0) SMAPJK(J-1,L) = FLAT/NI
            END DO
            WRITE(TITLEX(27:30),'(A)') CPRESS(L)
            TITLEL(L) = TITLEX//XLB
          END DO
        ELSEIF (name_ijk(K).eq.'HEIGHT') THEN ! special compound case
          DO L=1,LM
            DO J=2,JM
              NI = 0
              FLAT = 0.
              SMAPJK(J-1,L) = UNDEF
              DO I=1,IM
                DP=AIJK(I,J,L,IJK_DP)
                SMAP(I,J-1,L) = UNDEF
                IF(DP.GT.0.) THEN
                  SMAP(I,J-1,L) = SCALE_IJK(k)*(AIJK(I,J,L,IJK_DSE)-
     *                 SHA*AIJK(I,J,L,IJK_T))/DP
                  FLAT = FLAT+SMAP(I,J-1,L)
                  NI = NI+1
                END IF
              END DO
              IF (NI.GT.0) SMAPJK(J-1,L) = FLAT/NI
            END DO
            WRITE(TITLEX(27:30),'(A)') CPRESS(L)
            TITLEL(L) = TITLEX//XLB
          END DO
        END IF
        CALL POUT_IJK(TITLEL,name_ijk(k),lname_ijk(k),units_ijk(k)
     *       ,SMAP,SMAPJK,UNDEF)
      END DO
C****
      call close_ijk
C****
      RETURN
      END SUBROUTINE IJKMAP
