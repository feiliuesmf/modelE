      MODULE CLD01
!@sum  CLD01 column physics of moist conv. and large-scale condensation
!@auth M.S.Yao/A. Del Genio (modifications by Gavin Schmidt)
!@ver  1.0 (taken from CB265)
!@cont MSTCNV,LSCOND
      USE CONSTANT, only : rgas,grav,lhe,lhs,lhm,kapa,sha,bysha
     *     ,by3,tf,bytf,rvap,bygrav,deltx,bymrat
      USE MODEL_COM, only : IM,LM,DTsrc
      USE QUSDEF, only : nmom,xymoms,zmoms,zdir
      USE RANDOM
      IMPLICIT NONE
      SAVE
C**** parameters and constants
      REAL*8, PARAMETER :: TI=233.16d0   !@param TI pure ice limit
!@param WMU critical cloud water content for rapid conversion (g m**-3)
      REAL*8, PARAMETER :: WMU=.25       
      REAL*8, PARAMETER :: WMUL=.5       !@param WMUL WMU over land
      REAL*8, PARAMETER :: WMUI=.1d0     !@param WMUI WMU over ice
      REAL*8, PARAMETER :: BRCLD=.2d0    !@param BRCLD for cal. BYBR
      REAL*8, PARAMETER :: SLHE=LHE*BYSHA
      REAL*8, PARAMETER :: SLHS=LHS*BYSHA

      REAL*8 :: BYBR,BYDTsrc,XMASS
!@var BYBR factor for converting cloud particle radius to effect. radius
!@var XMASS dummy variable

C**** Set-able variables from NAMELIST
!@var LMCM max level for originating MC plumes (set in init_CLD)
      INTEGER :: LMCM
!@var U00wtr critical humidity for water cloud condensation (default)
      REAL*8 :: U00wtr = .7d0
!@var U00ice critical humidity for ice cloud condensation (default)
      REAL*8 :: U00ice = .7d0

C**** input variables
!@var RA ratio of primary grid box to secondary gridbox
      REAL*8, DIMENSION(IM) :: RA 
!@var UM,VM,U_0,V_0 velocity related variables (UM,VM)=(U,V)*AIRM
      REAL*8, DIMENSION(IM,LM) :: UM,VM 
      REAL*8, DIMENSION(IM,LM) :: U_0,V_0

!@var Miscellaneous vertical arrays set in driver
      REAL*8, DIMENSION(LM+1) :: PLE    !@var PLE pressure at layer edge
      REAL*8, DIMENSION(LM) :: PL,PLK,AIRM,BYAM,ETAL,TL,QL,TH,RH,WMX
     *     ,VSUBL,AJ8,AJ13,AJ50,AJ51,AJ52,AJ57,AQ,DPDT
!@var PL layer pressure (mb)
!@var PLK PL**KAPA
!@var AIRM the layer's pressure depth (mb)
!@var BYAM 1./AIRM
!@var ETAL fractional entrainment rate
!@var TL, QL temperature, specific humidity of the layer
!@var TH potential temperature (K)
!@var RH relative humidity
!@var WMX cloud water mixing ratio (kg/Kg)
!@var VSUBL downward vertical velocity due to cumulus subsidence (cm/s)
!@var AJ8, AJ13, AJ50, AJ52, AJ57 dummy variables
!@var AQ time change rate of specific humidity (s**-1)
!@var DPDT time change rate of pressure (mb/s)
      REAL*8, DIMENSION(LM+1) :: PRECNVL
!@var PRECNVL convective precip entering the layer top
!new arrays must be set to model arrays in driver (before MC)
      REAL*8, DIMENSION(LM) :: SDL,WML
!@var SDL vertical velocity in sigma coordinate
!@var WML cloud water mixing ratio (kg/Kg)
C**** new arrays must be set to model arrays in driver (after MC)
      REAL*8, DIMENSION(LM) :: TAUMCL,SVLATL,CLDMCL,SVLHXL,SVWMXL
!@var TAUMCL convective cloud optical thickness
!@var SVLATL saved LHX for convective cloud
!@var CLDMCL convective cloud cover
!@var SVLHXL saved LHX for large-scale cloud
!@var SVWMXL saved detrained convective cloud water
      REAL*8, DIMENSION(LM) :: CSIZEL
!@var CSIZEL cloud particle radius (micron)
C**** new arrays must be set to model arrays in driver (before COND)
      REAL*8, DIMENSION(LM) :: TTOLDL,CLDSAVL
!@var TTOLDL previous potential temperature
!@var CLDSAVL saved large-scale cloud cover
C**** new arrays must be set to model arrays in driver (after COND)
      REAL*8, DIMENSION(LM) :: AJ11,AJ55,TAUSSL,CLDSSL
!@var AJ11, AJ55 dummy variables
!@var TAUSSL large-scale cloud optical thickness
!@var CLDSSL large-scale cloud cover

!@var SM,QM Vertical profiles of T/Q
      DOUBLE PRECISION, DIMENSION(LM) :: SM,QM
      DOUBLE PRECISION, DIMENSION(NMOM,LM) :: SMOM,QMOM

!@var KMAX index for surrounding velocity
      INTEGER ::  KMAX
!@var PEARTH fraction of land in grid box
!@var TS average surface temperture (C)
!@var RI1, RI2 Richardson numbers
      REAL*8 :: PEARTH,TS,QS,US,VS,DCL,RI1,RI2
!@var LPBL max level of planetary boundary layer
      INTEGER :: LPBL

C**** output variables
      REAL*8 :: PRCPMC,PRCPSS,HCNDSS,WMSUM
!@var PRCPMC precip due to moist convection
!@var PRCPSS precip due to large-scale condensation
!@var HCNDSS heating due to large-scale condensation
!@var WMSUM cloud liquid water path
      REAL*8 :: CLDSLWIJ,CLDDEPIJ
!@var CLDSLWIJ shallow convective cloud cover
!@var CLDDEPIJ deep convective cloud cover
      INTEGER :: LMCMAX,LMCMIN
!@var LMCMAX upper-most convective layer
!@var LMCMIN lowerest convective layer
!@var AIRXL is convective mass flux (kg/m*m)
      REAL*8 AIRXL

C**** functions
      REAL*8 :: QSAT, DQSATDT
!@var QSAT saturation humidity
!@var DQSATDT dQSAT/dT

      CONTAINS

      SUBROUTINE MSTCNV
!@sum  MSTCNV moist convective processes (precip, convective clouds,...)
!@auth M.S.Yao/A. Del Genio (modularisation by Gavin Schmidt)
!@ver  1.0 (taken from CB265)
!@calls adv1d,QSAT,DQSATDT,THBAR
      IMPLICIT NONE
      REAL*8 LHX,MPLUME,MCLOUD,MPMAX,MPO,SENV,QENV
!@var LHX latent heat of evaporation (J/Kg)
!@var MPLUME mass of convective plume (mb)
!@var MCLOUD air mass available for re-evaporation of precip
!@var MPMAX convective plume at the detrainment level
!@var MPO old MPLUME
!@var SENV,QENV dummy variables

      REAL*8, DIMENSION(0:LM) :: CM     !@var CM air mass of subsidence
      REAL*8, DIMENSION(IM) :: UMP,VMP,UMDN,VMDN,UMPO,VMPO  !@var
!@var UMP, VMP momentum carried by convective plumes
!@var UMDN,VMDN dummy variables
!@var UMPO, VMPO old UMP, VMP
!@var DQM,DSM,DQMR,DSMR Vertical profiles of T/Q and changes
      REAL*8, DIMENSION(LM) ::
     * SMOLD,QMOLD, DQM,DSM,DQMR,DSMR
!@var SMOLD,QMOLD old SM, QM
      DOUBLE PRECISION, DIMENSION(LM) :: F,CMNEG
      DOUBLE PRECISION, DIMENSION(NMOM,LM) :: FMOM
      DOUBLE PRECISION, DIMENSION(NMOM,LM) :: SMOMOLD,QMOMOLD
      DOUBLE PRECISION, DIMENSION(NMOM,LM) :: DSMOM,DQMOM,DSMOMR,DQMOMR
      DOUBLE PRECISION, DIMENSION(NMOM) ::
     &     SMOMP,QMOMP, SMOMPMAX,QMOMPMAX, SMOMDN,QMOMDN

      REAL*8, DIMENSION(LM) ::
     *     DM,COND,CDHEAT,CCM,SM1,QM1,DMR,ML,SMT,QMT,TPSAV,SVTP,DDM
!@var DM change in air mass
!@var COND condensate
!@var CDHEAT heating due to condensation
!@var CCM convective plume mass
!@var SM1, QM1 dummy variables
!@var DMR change of air mass
!@var ML layer air mass
!@var SMT, QMT dummy variables
!@var TPSAV, SVTP  arrays to save plume temperature
!@var DDM downdraft mass

      INTEGER LDRAFT,LEVAP,LMAX,LMIN,MCCONT,MAXLVL
     *     ,MINLVL,ITER,IC,LFRZ,NSUB,LDMIN
!@var LDRAFT the layer the downdraft orginates
!@var LEVAP the layer evaporation of precip starts
!@var LMAX, LMIN the lowest, the highest layer of a convective event
!@var MCCONT interger to count convective events
!@var MAXLVL, MINLVL the lowest, the highest layer of convective events
!@var ITER number for iteration
!@var IC interger for cloud types
!@var LFRZ freezing level
!@var nsub = LMAX - LMIN + 1
!@var LDMIN the lowest layer the downdraft drops
      REAL*8 TERM1,FMP0,SMO1
     *     ,QMO1,SMO2,QMO2,SDN,QDN,SUP,QUP,SEDGE,QEDGE,WMDN,WMUP,SVDN
     *     ,SVUP,WMEDG,SVEDG,DMSE,FPLUME,DFP,FMP2,FRAT1,FRAT2,SMN1
     *     ,QMN1,SMN2,QMN2,SMP,QMP,TP,GAMA,DQSUM,TNX
     *     ,DQ,DMSE1,FCTYPE,BETAU,ALPHAU
     *     ,CDHDRT,DDRAFT,DELTA
     *     ,ALPHA,BETA,CDHM,CDHSUM,CLDM,CLDREF,CONSUM,DQEVP
     *     ,DQRAT,EPLUME,ETADN,ETAL1,EVPSUM,FCDH
     *     ,FCDH1,FCLD,FCLOUD,FDDL,FDDP,FENTR,FENTRA,FEVAP,FLEFT
     *     ,FQCOND,FSEVP,FSSUM,HEAT1
     *     ,PRHEAT,PRCP
     *     ,QMDN,QMIX,QMPMAX,QMPT,QNX,QSATC,QSATMP
     *     ,RCLD,RCLDE,SLH,SMDN,SMIX,SMPMAX,SMPT,SUMAJ
     *     ,SUMDP,DDRUP,EDRAFT
     *     ,TOLD,TOLD1,TEMWM,TEM,WTEM,WCONST,WORK,SMPO,QMPO
!@var TERM1 contribution to non-entraining convective cloud
!@var FMP0 non-entraining convective mass
!@var SMO1,QMO1,SMO2,QMO2,SDN,QDN,SUP,QUP,SEDGE,QEDGE dummy variables
!@var WMDN,WMUP,SVDN,SVUP,WMEDG,SVEDG,DDRUP dummy variables
!@var DMSE difference in moist static energy
!@var FPLUME fraction of convective plume
!@var DFP an iterative increment
!@var FMP2,FRAT1,FRAT2,SMN1,QMN1,SMN2,QMN2 dummy variables
!@var SMP, QMP plume's SM, QM
!@var TP plume's temperature (K)
!@var GAMA,DQSUM,TNX,DQ,CONSUM,BETAU,ALPHAU dummy variables
!@var DMSE1 difference in moist static energy
!@var FCTYPE fraction for convective cloud types
!@var CDHDRT,DELTA.ALPHA,BETA,CDHM,CDHSUM,CLDREF dummy variables
!@var DDRAFT downdraft mass
!@var DELTA fraction of plume stays in the layer
!@var CLDM subsidence due to convection
!@var DQEVP amount of condensate that evaporates in downdrafts
!@var DQRAT fraction for the condensate to evaporate
!@var EPLUME air mass of entrainment
!@var ETADN fraction of the downdraft
!@var ETAL1 fractional entrainment rate
!@var EVPSUM,FCDH,FCDH1,FCLD,FCLOUD,FDDL,FDDP dummy variables
!@var FENTR fraction of entrainment
!@var FENTRA fraction of entrainment
!@var FEVAP fraction of air available for precip evaporation
!@var FLEFT fraction of plume after removing downdraft mass
!@var FQCOND
!@var FSEVP
!@var FSSUM
!@var HEAT1 heating due to phase change
!@var PRHEAT energy of condensate
!@car PRCP precipipation
!@var QMDN,QMIX,SMDN,SMIX dummy variables
!@var QMPMAX,SMPMAX detrainment of QMP, SMP
!@var QMPT,SMPT dummy variables
!@var QNX,SUMAJ,SUMDP dummy variables
!@var QSATC saturation vapor mixing ratio
!@var QSATMP plume's saturation vapor mixing ratio
!@var RCLD,RCLDE cloud particle's radius, effective radius
!@var SLH LHX/SHA
!@var EDRAFT entrainment into downdrafts
!@var TOLD,TOLD1 old temperatures
!@var TEMWM,TEM,WTEM,WCONST dummy variables
!@var WORK work done on convective plume
!@var SMPO,QMPO old SMP, QMP

      LOGICAL MC1  !@var MC1 true for the first convective event

      REAL*8,  PARAMETER :: CK1 = 1.       !@param CK1 a tunning const.

      INTEGER K,L,N  !@var K,L,N loop variables
      INTEGER ITYPE  !@var convective cloud types
!@var DUM, DVM changes of UM,VM
      REAL*8, DIMENSION(IM,LM) :: DUM,DVM 

      REAL*8 THBAR  !@var THBAR virtual temperature at layer edge
C****
C**** MOIST CONVECTION
C****
C**** CONVECTION OCCURS AT THE LOWEST MOIST CONVECTIVELY UNSTABLE
C**** LEVEL AND CONTINUES UNTIL A STABLE LAYER PAIR IS REACHED.  RE-
C**** EVAPORATION AND PRECIPITATION ARE COMPUTED AT THE END OF THIS
C**** CYCLE.  THE WHOLE PROCESS MAY BE REPEATED FROM A NEW LOWEST
C**** UNSTABLE LEVEL.
C****
      LMCMIN=0
      LMCMAX=0
      MCCONT=0
C**** initiallise arrays of computed ouput
      TAUMCL=0
      SVWMXL=0
      SVLATL=0
      VSUBL=0
      PRECNVL=0
      CLDMCL=0
      CLDSLWIJ=0
      CLDDEPIJ=0
      PRCPMC=0.
      SVTP=0
      CSIZEL=10.    !  effective droplet radius in stem of convection
C**** zero out diagnostics
         AJ8 =0.
         AJ13=0.
         AJ50=0.
         AJ51=0.
         AJ52=0.
         AJ57=0.
C**** save initial values
      SM1=SM
      QM1=QM
C**** OUTER MC LOOP OVER BASE LAYER
      DO 600 LMIN=1,LMCM-1
      MAXLVL=0
      MINLVL=LM
C****
C**** COMPUTE THE CONVECTIVE MASS OF THE NON-ENTRAINING PART
C****
      TERM1=-10.*CK1*SDL(LMIN+1)*BYGRAV
      FMP0=TERM1*XMASS
      IF(FMP0.LE.0.) FMP0=0.
C**** CREATE A PLUME IN THE BOTTOM LAYER
C****
C**** ITERATION TO FIND FPLUME WHICH RESTORES THE ATM TO NEUTRAL STATE
C****
      SMO1=SM(LMIN)
      QMO1=QM(LMIN)
      SMO2=SM(LMIN+1)
      QMO2=QM(LMIN+1)
      SDN=SMO1*BYAM(LMIN)
      SUP=SMO2*BYAM(LMIN+1)
      SEDGE=THBAR(SUP,SDN)
      QDN=QMO1*BYAM(LMIN)
      QUP=QMO2*BYAM(LMIN+1)
      WMDN=WML(LMIN)
      WMUP=WML(LMIN+1)
      SVDN=SDN*(1.+DELTX*QDN-WMDN)
      SVUP=SUP*(1.+DELTX*QUP-WMUP)
      QEDGE=.5*(QUP+QDN)
      WMEDG=.5*(WMUP+WMDN)
      SVEDG=SEDGE*(1.+DELTX*QEDGE-WMEDG)
      LHX=LHE
      SLH=LHX*BYSHA
      DMSE=(SVUP-SVEDG)*PLK(LMIN+1)+(SVEDG-SVDN)*PLK(LMIN)+
     *  SLHE*(QSAT(SUP*PLK(LMIN+1),LHX,PL(LMIN+1))-QDN)
      IF(DMSE.GT.-1.E-10) GO TO 600
C****
      FPLUME=.25
      DFP = .25
      DO ITER=1,8
      DFP=DFP*0.5
      FMP2=FPLUME*AIRM(LMIN)
      FRAT1=FMP2*BYAM(LMIN+1)
      FRAT2=FMP2*BYAM(LMIN+2)
      SMN1=SMO1*(1.-FPLUME)+FRAT1*SMO2
      QMN1=QMO1*(1.-FPLUME)+FRAT1*QMO2
      SMN2=SMO2*(1.-FRAT1)+FRAT2*SM(LMIN+2)
      QMN2=QMO2*(1.-FRAT1)+FRAT2*QM(LMIN+2)
      SMP=SMO1*FPLUME
      QMP=QMO1*FPLUME
      TP=SMO1*PLK(LMIN+1)*BYAM(LMIN)
      QSATMP=FMP2*QSAT(TP,LHX,PL(LMIN+1))
      GAMA=SLH*QSATMP*DQSATDT(TP,LHX)/FMP2
      DQSUM=(QMP-QSATMP)/(1.+GAMA)
      IF(DQSUM.LE.0.) GO TO 205
      FEVAP=.5*FPLUME
      MCLOUD=FEVAP*AIRM(LMIN+1)
      TNX=SMO2*PLK(LMIN+1)*BYAM(LMIN+1)
      QNX=QMO2*BYAM(LMIN+1)
      QSATC=QSAT(TNX,LHX,PL(LMIN+1))
      DQ=MCLOUD*(QSATC-QNX)/(1.+SLH*QSATC*DQSATDT(TNX,LHX))
      IF(DQ.GT.DQSUM) DQ=DQSUM
      SMN2=SMN2-SLH*DQ/PLK(LMIN+1)
      QMN2=QMN2+DQ
      DQSUM=DQSUM-DQ
      IF(DQSUM.LE.0.) GO TO 205
      MCLOUD=FEVAP*AIRM(LMIN)
      TNX=SMO1*PLK(LMIN)*BYAM(LMIN)
      QNX=QMO1*BYAM(LMIN)
      QSATC=QSAT(TNX,LHX,PL(LMIN))
      DQ=MCLOUD*(QSATC-QNX)/(1.+SLH*QSATC*DQSATDT(TNX,LHX)) ! new func'n
      IF(DQ.GT.DQSUM) DQ=DQSUM
      SMN1=SMN1-SLH*DQ/PLK(LMIN)
      QMN1=QMN1+DQ
  205 SDN=SMN1*BYAM(LMIN)
      SUP=SMN2*BYAM(LMIN+1)
      SEDGE=THBAR(SUP,SDN)
      QDN=QMN1*BYAM(LMIN)
      QUP=QMN2*BYAM(LMIN+1)
      SVDN=SDN*(1.+DELTX*QDN-WMDN)
      SVUP=SUP*(1.+DELTX*QUP-WMUP)
      QEDGE=.5*(QUP+QDN)
      SVEDG=SEDGE*(1.+DELTX*QEDGE-WMEDG)
      DMSE1=(SVUP-SVEDG)*PLK(LMIN+1)+(SVEDG-SVDN)*PLK(LMIN)+
     *  SLHE*(QSAT(SUP*PLK(LMIN+1),LHX,PL(LMIN+1))-QDN)
      IF (ABS(DMSE1).LE.1.d-3) GO TO 411
      IF(DMSE1.GT.1.d-3) FPLUME=FPLUME-DFP
      IF(DMSE1.LT.-1.d-3) FPLUME=FPLUME+DFP
      END DO
  411 IF(FPLUME.LE..001) GO TO 600
C****
C**** ITERATION THROUGH CLOUD TYPES
C****
      ITYPE=2                        ! always 2 types of clouds:
C     IF(LMIN.LE.2) ITYPE=2          ! entraining & non-entraining
      FCTYPE=1.
C**** SET PROFILE TO BE CONSTANT FOR BOTH TYPES OF CLOUDS
      SMOLD(:) = SM(:)
      SMOMOLD(:,:) = SMOM(:,:)
      QMOLD(:) = QM(:)
      QMOMOLD(:,:) = QMOM(:,:)
      DO 570 IC=1,ITYPE
C**** INITIALLISE VARIABLES USED FOR EACH TYPE
      DO L=1,LM
        COND(L)=0.
        CDHEAT(L)=0.
        DM(L)=0.
        DMR(L)=0.
        CCM(L)=0.
        DDM(L)=0.
      END DO
      DUM(1:KMAX,:)=0.
      DVM(1:KMAX,:)=0.
      DSM(:) = 0.
      DSMOM(:,:) = 0.
      DSMR(:) = 0.
      DSMOMR(:,:) = 0.
      DQM(:) = 0.
      DQMOM(:,:) = 0.
      DQMR(:) = 0.
      DQMOMR(:,:) = 0.
      MC1=.FALSE.
      LHX=LHE
      MPLUME=MIN(1.*AIRM(LMIN),1.*AIRM(LMIN+1))
      IF(MPLUME.GT.FMP2) MPLUME=FMP2
      IF(ITYPE.EQ.2) THEN
      FCTYPE=1.
      IF(MPLUME.GT.FMP0) FCTYPE=FMP0/MPLUME
      IF(IC.EQ.2) FCTYPE=1.-FCTYPE
      IF(FCTYPE.LT.0.001) GO TO 570
      END IF
      MPLUME=MPLUME*FCTYPE
C     FPLUM0=FMP1*BYAM(LMIN)
      FPLUME=MPLUME*BYAM(LMIN)
      SMP  =  SMOLD(LMIN)*FPLUME
      SMOMP(xymoms)=SMOMOLD(xymoms,LMIN)*FPLUME
      QMP  =  QMOLD(LMIN)*FPLUME
      QMOMP(xymoms)=QMOMOLD(xymoms,LMIN)*FPLUME
      TPSAV(LMIN)=SMP*PLK(LMIN)/MPLUME
      DMR(LMIN)=-MPLUME
        DSMR(LMIN)=-SMP
      DSMOMR(xymoms,LMIN)=-SMOMP(xymoms)
      DSMOMR(zmoms,LMIN)=-SMOMOLD(zmoms,LMIN)*FPLUME
        DQMR(LMIN)=-QMP
      DQMOMR(xymoms,LMIN)=-QMOMP(xymoms)
      DQMOMR(zmoms,LMIN)=-QMOMOLD(zmoms,LMIN)*FPLUME
      DO K=1,KMAX !vref
         UMP(K)=UM(K,LMIN)*FPLUME !vref
         DUM(K,LMIN)=-UMP(K) !vref
         VMP(K)=VM(K,LMIN)*FPLUME !vref
         DVM(K,LMIN)=-VMP(K) !vref
      ENDDO !vref
C****
C**** RAISE THE PLUME TO THE TOP OF CONVECTION AND CALCULATE
C**** ENTRAINMENT, CONDENSATION, AND SECONDARY MIXING
C****
      CDHSUM=0.
      CDHDRT=0.
      ETADN=0.
      LDRAFT=LM
      EVPSUM=0.
      DDRAFT=0.
      LFRZ=0
      LMAX=LMIN
 220  L=LMAX+1
      SVTP(L)=SMP*PLK(L)/(MPLUME+1.E-20)
C**** TEST FOR SUFFICIENT AIR, MOIST STATIC STABILITY AND ENERGY
C     IF(L.GT.LMIN+1.AND.SDL(L).GT.0.) GO TO 340
      IF(MPLUME.LE..001*AIRM(L)) GO TO 340
      SDN=SMP/MPLUME
      SUP=SM1(L)*BYAM(L)
      QDN=QMP/MPLUME
      QUP=QM1(L)*BYAM(L)
      WMDN=0.
      WMUP=WML(L)
      SVDN=SDN*(1.+DELTX*QDN-WMDN)
      SVUP=SUP*(1.+DELTX*QUP-WMUP)
      IF(LMAX.GT.LMIN) THEN
      SEDGE=THBAR(SUP,SDN)
      QEDGE=.5*(QUP+QDN)
      WMEDG=.5*(WMUP+WMDN)
      SVEDG=SEDGE*(1.+DELTX*QEDGE-WMEDG)
      LHX=LHE
      DMSE=(SVUP-SVEDG)*PLK(L)+(SVEDG-SVDN)*PLK(L-1)+
     *  SLHE*(QSAT(SUP*PLK(L),LHX,PL(L))-QDN)
      IF(DMSE.GT.-1.E-10) GO TO 340
      END IF
      IF(PLK(L-1)*(SVUP-SVDN)+SLHE*(QUP-QDN).GE.0.) GO TO 340
C****
C**** DEPOSIT PART OF THE PLUME IN LOWER LAYER
C****
      DELTA=0.
      SMPO =SMP
      QMPO =QMP
      MPO = MPLUME
      DO K=1,KMAX !vref
         UMPO(K)=UMP(K) !vref
         VMPO(K)=VMP(K) !vref
      ENDDO !vref
C     IF(MPLUME.GT.AIRM(L)) THEN
      IF(MPLUME.GT..95*AIRM(L)) THEN
C     DELTA=(MPLUME-AIRM(L))/MPLUME
      DELTA=(MPLUME-.95*AIRM(L))/MPLUME
      SMP = SMP  *(1.-DELTA)
      QMP = QMP  *(1.-DELTA)
C     MPLUME=AIRM(L)
      MPLUME=.95*AIRM(L)
      DO K=1,KMAX !vref
         UMP(K)=UMP(K)-UMP(K)*DELTA !vref
         VMP(K)=VMP(K)-VMP(K)*DELTA !vref
      ENDDO !vref
      END IF
C****
C**** CONVECTION IN UPPER LAYER   (WORK DONE COOLS THE PLUME)
C****
      WORK=MPLUME*(SUP-SDN)*(PLK(L-1)-PLK(L))/PLK(L-1)
C     SMP=SMP-WORK
      DSM(L-1)=DSM(L-1)-WORK
      CCM(L-1)=MPLUME
      DM(L-1)=DM(L-1)+DELTA*MPO
C**** TEST FOR CONDENSATION ALSO DETERMINES IF PLUME REACHES UPPER LAYER
      TP=SMP*PLK(L)/MPLUME
      TPSAV(L)=TP
      IF(TPSAV(L-1).GE.TF.AND.TPSAV(L).LT.TF) LFRZ=L-1
      IF(TP.LT.TI) LHX=LHS
      QSATMP=MPLUME*QSAT(TP,LHX,PL(L))
      IF(QMP.LT.QSATMP) GO TO 340
      IF(TP.GE.TF.OR.LHX.EQ.LHS) GO TO 290
      LHX=LHS
      QSATMP=MPLUME*QSAT(TP,LHX,PL(L))
  290 SLH=LHX*BYSHA
        DSM(L-1)=  DSM(L-1)+DELTA*SMPO
      DSMOM(xymoms,L-1)=DSMOM(xymoms,L-1)+DELTA*SMOMP(xymoms)
c      SMP = SMP *(1.-DELTA)     ! already set above
      SMOMP(xymoms) = SMOMP(xymoms)*(1.-DELTA)
        DQM(L-1)=  DQM(L-1)+DELTA*QMPO
      DQMOM(xymoms,L-1)=DQMOM(xymoms,L-1)+DELTA*QMOMP(xymoms)
c      QMP = QMP *(1.-DELTA)     ! already set above
      QMOMP(xymoms) = QMOMP(xymoms)*(1.-DELTA)
      DO K=1,KMAX !vref
         DUM(K,L-1)=DUM(K,L-1)+UMPO(K)*DELTA !vref
         DVM(K,L-1)=DVM(K,L-1)+VMPO(K)*DELTA !vref
      ENDDO !vref
C****
C**** ENTRAINMENT
C****
      IF(IC.EQ.2.OR.(IC.EQ.1.AND.PL(L).GE.800.)) THEN
      FENTR=ETAL(L)*FPLUME
      IF(FENTR+FPLUME.GT.1.) FENTR=1.-FPLUME
      ETAL1=FENTR/(FPLUME+1.E-20)
      FPLUME=FPLUME+FENTR
      EPLUME=MPLUME*ETAL1
C**** Reduce EPLUME so that mass flux is less than mass in box
      IF (EPLUME.GT.AIRM(L)*0.975d0-MPLUME) THEN
        EPLUME=AIRM(L)*0.975d0-MPLUME
      END IF
      MPLUME=MPLUME+EPLUME
      FENTRA = EPLUME*BYAM(L)
      DSMR(L)=DSMR(L)-EPLUME*SUP        ! = DSM(L)-SM(L)*FENTRA
      DSMOMR(:,L)=DSMOMR(:,L)-SMOM(:,L)*FENTRA
      DQMR(L)=DQMR(L)-EPLUME*QUP        ! = DQM(L)-QM(L)*FENTRA
      DQMOMR(:,L)=DQMOMR(:,L)-QMOM(:,L)*FENTRA
      DMR(L)=DMR(L)-EPLUME
      SMP=SMP+EPLUME*SUP
       SMOMP(xymoms)= SMOMP(xymoms)+ SMOM(xymoms,L)*FENTRA
      QMP=QMP+EPLUME*QUP
       QMOMP(xymoms)= QMOMP(xymoms)+ QMOM(xymoms,L)*FENTRA
      DO K=1,KMAX !vref
         UMP(K)=UMP(K)+U_0(K,L)*EPLUME !vref
         DUM(K,L)=DUM(K,L)-U_0(K,L)*EPLUME !vref
         VMP(K)=VMP(K)+V_0(K,L)*EPLUME !vref
         DVM(K,L)=DVM(K,L)-V_0(K,L)*EPLUME !vref
      ENDDO !vref
      END IF
C****
C**** CHECK THE DOWNDRAFT POSSIBILITY
C****
      IF(L-LMIN.LE.1) GO TO 291
      IF(ETADN.GT.1.E-10) GO TO 291
      SMIX=.5*(SUP+SMP/MPLUME)
      QMIX=.5*(QUP+QMP/MPLUME)
C     WMIX=.5*(WMUP+WMDN)
C     SVMIX=SMIX*(1.+DELTX*QMIX-WMIX)
C     DMMIX=(SVUP-SVMIX)*PLK(L)+
C    *  SLHE*(QSAT(SUP*PLK(L),LHX,PL(L))-QMIX)
C     IF(DMMIX.LT.1.E-10) GO TO 291
      IF(SMIX.GE.SUP) GO TO 291
      IF(PL(L).GT.700.) GO TO 291
      LDRAFT=L
      ETADN=BY3
      FLEFT=1.-.5*ETADN
      DDRAFT=ETADN*MPLUME
      FDDP = .5*DDRAFT/MPLUME
      FDDL = .5*DDRAFT*BYAM(L)
      MPLUME=FLEFT*MPLUME
      SMDN=DDRAFT*SMIX         ! = SM(L)*FDDL +  SMP*FDDP
      SMOMDN(xymoms)=SMOM(xymoms,L)*FDDL +  SMOMP(xymoms)*FDDP
      SMP=FLEFT*SMP
       SMOMP(xymoms)= SMOMP(xymoms)*FLEFT
      QMDN=DDRAFT*QMIX        ! = QM(L)*FDDL +  QMP*FDDP
      QMOMDN(xymoms)=QMOM(xymoms,L)*FDDL +  QMOMP(xymoms)*FDDP
      QMP=FLEFT*QMP
       QMOMP(xymoms)= QMOMP(xymoms)*FLEFT
      DMR(L) = DMR(L)-.5*DDRAFT
      DSMR(L)=DSMR(L)-.5*DDRAFT*SUP        ! = DSM(L)-SM(L)*FDDL
      DSMOMR(:,L)=DSMOMR(:,L) - SMOM(:,L)*FDDL
      DQMR(L)=DQMR(L)-.5*DDRAFT*QUP       ! = DQM(L)-QM(L)*FDDL
      DQMOMR(:,L)=DQMOMR(:,L) - QMOM(:,L)*FDDL
      DO K=1,KMAX !vref
         UMDN(K)=.5*(ETADN*UMP(K)+DDRAFT*U_0(K,L)) !vref
         UMP(K)=UMP(K)*FLEFT !vref
         DUM(K,L)=DUM(K,L)-.5*DDRAFT*U_0(K,L) !vref
         VMDN(K)=.5*(ETADN*VMP(K)+DDRAFT*V_0(K,L)) !vref
         VMP(K)=VMP(K)*FLEFT !vref
         DVM(K,L)=DVM(K,L)-.5*DDRAFT*V_0(K,L) !vref
      ENDDO !vref
c         TMIX=SMIX*PLK(L)
c         HMIX=SHA*TMIX+LHE*QMIX
c         TENV=SUP*PLK(L)
c         QENV=QUP
c         HENV=SHA*TENV+LHE*QENV
c         HSENV=SHA*TENV+LHE*QSAT(TENV,LHX,PL(L))
C        IF (I.EQ.35.AND.J.EQ.13) WRITE(6,299)
C    *     LDRAFT,TMIX,QMIX,HMIX,TENV,QENV,HENV,HSENV
c  299 FORMAT(1X,'LDD TMIX QMIX HMIX TENV QENV HENV HSENV =',
c     *  I3,7E12.4/)
C****
C**** CONDENSE VAPOR IN THE PLUME AND ADD LATENT HEAT
C****
  291 DQSUM=0.
      SMPT=SMP
      QMPT=QMP
      DO 292 N=1,3
      TP=SMP*PLK(L)/MPLUME
      QSATMP=MPLUME*QSAT(TP,LHX,PL(L))
      GAMA=SLH*QSATMP*DQSATDT(TP,LHX)/MPLUME
      DQ=(QMP-QSATMP)/(1.+GAMA)
      SMP=SMP+SLH*DQ/PLK(L)
      QMP=QMP-DQ
  292 DQSUM=DQSUM+DQ
      IF(DQSUM.GE.0.) THEN
      FQCOND = 0
      IF (QMPT.gt.1d-20) FQCOND = DQSUM/QMPT
       QMOMP(xymoms) =  QMOMP(xymoms)*(1.-FQCOND)
      ELSE  ! no change
        DQSUM=0.
        SMP=SMPT
        QMP=QMPT
      END IF
      COND(L)=DQSUM
      TAUMCL(L)=TAUMCL(L)+DQSUM
      CDHEAT(L)=SLH*COND(L)
      CDHSUM=CDHSUM+CDHEAT(L)
      IF(ETADN.GT.1.E-10) CDHDRT=CDHDRT+SLH*COND(L)
C****
C**** UPDATE ALL QUANTITIES CARRIED BY THE PLUME
C****
      MCCONT=MCCONT+1
      IF(MCCONT.EQ.1) MC1=.TRUE.
      IF(MC1.AND.PLE(LMIN)-PLE(L+2).GE.450.) SVLATL(L)=LHX
      SMPMAX=SMP
      SMOMPMAX(xymoms) =  SMOMP(xymoms)
      QMPMAX=QMP
      QMOMPMAX(xymoms) =  QMOMP(xymoms)
      MPMAX=MPLUME
      LMAX = LMAX + 1
      IF (LMAX.LT.LM) GO TO 220   ! CHECK FOR NEXT POSSIBLE LMAX
C**** UPDATE CHANGES CARRIED BY THE PLUME IN THE TOP CLOUD LAYER
  340 IF(LMIN.EQ.LMAX) GO TO 600
      IF(TPSAV(LMAX).GE.TF) LFRZ=LMAX
      DM(LMAX)=DM(LMAX)+MPMAX
      DSM(LMAX)=DSM(LMAX)+SMPMAX
      DSMOM(xymoms,LMAX)=DSMOM(xymoms,LMAX) + SMOMPMAX(xymoms)
      DQM(LMAX)=DQM(LMAX)+QMPMAX
      DQMOM(xymoms,LMAX)=DQMOM(xymoms,LMAX) + QMOMPMAX(xymoms)
      CCM(LMAX)=0.
      DO K=1,KMAX !vref
         DUM(K,LMAX)=DUM(K,LMAX)+UMP(K) !vref
         DVM(K,LMAX)=DVM(K,LMAX)+VMP(K) !vref
      ENDDO !vref
      CDHM=0.
      IF(MINLVL.GT.LMIN) MINLVL=LMIN
      IF(MAXLVL.LT.LMAX) MAXLVL=LMAX
      IF(LMCMIN.EQ.0) LMCMIN=LMIN
      IF(LMCMAX.LT.MAXLVL) LMCMAX=MAXLVL
C     IF(.NOT.MC1) GO TO 344
C     IF(IC.EQ.ITYPE) MC1=.FALSE.
C 344 CONTINUE
C****
C**** PROCESS OF DOWNDRAFTING
C****
      LDMIN=LMIN
      EDRAFT=0.
      IF(ETADN.GT.1.E-10) THEN
      CONSUM=0.
      DO 347 L=LMIN,LDRAFT
  347 CONSUM=CONSUM+COND(L)
      TNX=SMDN*PLK(LMIN)/DDRAFT
      QNX=QMDN/DDRAFT
      LHX=LHE
      IF(TPSAV(LMIN).LT.TF) LHX=LHS
      SLH=LHX*BYSHA
      QSATC=QSAT(TNX,LHX,PL(LMIN))
      DQ=(QSATC-QNX)/(1.+SLH*QSATC*DQSATDT(TNX,LHX))
      DQRAT=DQ*DDRAFT/(CONSUM+1.E-20)
      DO 346 L=LDRAFT,1,-1
C     TNX=SMDN*PLK(L)/DDRAFT
C     QNX=QMDN/DDRAFT
      LHX=LHE
      IF(TPSAV(L).LT.TF) LHX=LHS
      SLH=LHX*BYSHA
C     QSATC=QSAT(TNX,LHX,PL(L))
c     DQ=(QSATC-QNX)/(1.+SLH*QSATC*DQSATDT(TNX,LHX))
C     IF(DQ.LT.0.) DQ=0.
C     DQEVP=DQ*DDRAFT
      DQEVP=DQRAT*COND(L)
      IF(DQEVP.GT.COND(L)) DQEVP=COND(L)
      IF (L.LT.LMIN) DQEVP=0.
      FSEVP = 0
      IF (ABS(PLK(L)*SMDN).gt.1d-20) FSEVP = SLH*DQEVP/(PLK(L)*SMDN)
      SMDN=SMDN-SLH*DQEVP/PLK(L)
      SMOMDN(xymoms)=SMOMDN(xymoms)*(1.-FSEVP)
      QMDN=QMDN+DQEVP
      COND(L)=COND(L)-DQEVP
      TAUMCL(L)=TAUMCL(L)-DQEVP
      CDHEAT(L)=CDHEAT(L)-DQEVP*SLH
      EVPSUM=EVPSUM+DQEVP*SLH
C**** ENTRAINMENT OF DOWNDRAFTS
      IF(L.LT.LDRAFT.AND.L.GT.1) THEN
        DDRUP=DDRAFT
        DDRAFT=DDRAFT*(1.+ETAL(L))
C       IF(L.GT.LMIN.AND.DDRAFT.GT.CCM(L-1)) DDRAFT=CCM(L-1)
        IF(DDRUP.GT.DDRAFT) DDRAFT=DDRUP
        IF(DDRAFT.GT..95*(AIRM(L-1)+DMR(L-1)))
     *    DDRAFT=.95*(AIRM(L-1)+DMR(L-1))
        EDRAFT=DDRAFT-DDRUP
        FENTRA=EDRAFT*BYAM(L)
        SENV=SM1(L)/AIRM(L)
        QENV=QM1(L)/AIRM(L)
        SMDN=SMDN+EDRAFT*SENV
        QMDN=QMDN+EDRAFT*QENV
        SMOMDN(xymoms)= SMOMDN(xymoms)+ SMOM(xymoms,L)*FENTRA
        QMOMDN(xymoms)= QMOMDN(xymoms)+ QMOM(xymoms,L)*FENTRA
        DSMR(L)=DSMR(L)-EDRAFT*SENV
        DSMOMR(:,L)=DSMOMR(:,L)-SMOM(:,L)*FENTRA
        DQMR(L)=DQMR(L)-EDRAFT*QENV
        DQMOMR(:,L)=DQMOMR(:,L)-QMOM(:,L)*FENTRA
        DMR(L)=DMR(L)-EDRAFT
      ENDIF
      IF(L.GT.1) DDM(L-1)=DDRAFT
c         TMIX=SMDN*PLK(L)/DDRAFT
c         QMIX=QMDN/DDRAFT
c         HMIX=SHA*TMIX+LHX*QMIX
c         TENV=SM1(L)*PLK(L)*BYAM(L)
c         QENV=QM1(L)*BYAM(L)
c         HENV=SHA*TENV+LHX*QENV
c         HSENV=SHA*TENV+LHX*QSAT(TENV,LHX,PL(L))
C     IF(I.EQ.35.AND.J.EQ.13) WRITE (6,399)
C    *  L,TMIX,QMIX,HMIX,TENV,QENV,HENV,HSENV
c  399 FORMAT(1X,'L TMIX QMIX HMIX TENV QENV HENV HSENV=',
c     *  I5,7E12.4)
      LDMIN=L
C     ALLOW FOR DOWNDRAFT TO DROP BELOW LMIN, IF IT'S NEGATIVE BUOYANT
      IF (L.LE.LMIN.AND.L.GT.1) THEN
        SMIX=SMDN/DDRAFT
        IF (SMIX.GE.(SM1(L-1)/AIRM(L-1))) GO TO 345
      ENDIF
  346 CONTINUE
  345 DSM(LDMIN)=DSM(LDMIN)+SMDN
      DSMOM(xymoms,LDMIN)=DSMOM(xymoms,LDMIN) + SMOMDN(xymoms)
      DQM(LDMIN)=DQM(LDMIN)+QMDN
      DQMOM(xymoms,LDMIN)=DQMOM(xymoms,LDMIN) + QMOMDN(xymoms)
      DO K=1,KMAX !vref
      DUM(K,LDMIN)=DUM(K,LDMIN)+UMDN(K) !vref
      DVM(K,LDMIN)=DVM(K,LDMIN)+VMDN(K) !vref
      ENDDO !vref
      DM(LDMIN)=DM(LDMIN)+DDRAFT
      END IF
C****
C**** SUBSIDENCE AND MIXING
C****
C**** Calculate vertical mass fluxes (Note CM for subsidence is defined
C**** in opposite sense than normal (positive is down))
      DO L=0,LDMIN-1
        CM(L) = 0.
      END DO
      DO L=LDMIN,LMAX
        CM(L) = CM(L-1) - DM(L) - DMR(L)
        SMT(L)=SM(L)    ! Save profiles for diagnostics
        QMT(L)=QM(L)
      END DO
      DO L=LMAX+1,LM
        CM(L) = 0.
      END DO
C**** simple upwind scheme for momentum
      ALPHA=0.
      DO 380 L=LDMIN,LMAX
      CLDM=CCM(L)
      IF(L.LT.LDRAFT.AND.ETADN.GT.1.E-10) CLDM=CCM(L)-DDM(L)
      IF(MC1) VSUBL(L)=100.*CLDM*RGAS*TL(L)/(PL(L)*GRAV*DTsrc)
      BETA=CLDM*BYAM(L+1)
      IF(CLDM.LT.0.) BETA=CLDM*BYAM(L)
c         FCDH=0.
c         IF(L.EQ.LMAX) FCDH=CDHSUM-(CDHSUM-CDHDRT)*.5*ETADN+CDHM
C***     IF(L.EQ.LDEP) FCDH=FCDH+CDHM
c         FCDH1=0.
c         IF(L.EQ.LMIN) FCDH1=(CDHSUM-CDHDRT)*.5*ETADN-EVPSUM
c         AJ8(L)=AJ8(L)+CCM(L)
c       AJ13(L)=AJ13(L)+PLK(L)*(-ALPHA*SM(L)+BETA*SM(L+1)+DSM(L))
c   *                  -FCDH-FCDH1
c       AJ51(L)=AJ51(L)+SLHE*(-ALPHA*QM(L)+BETA*QM(L+1)+DQM(L)+
c   *           COND(L))
c       AJ57(L)=AJ57(L)+SLHE*(-ALPHA*QM(L)+BETA*QM(L+1)+DQM(L))
c     SM(L)=SM(L)*(1.-ALPHA)+BETA*SM(L+1)+DSM(L)
c     QM(L)=QM(L)*(1.-ALPHA)+BETA*QM(L+1)+DQM(L)
      BETAU=BETA
      ALPHAU=ALPHA
      IF(BETA.LT.0.) BETAU=0.
      IF(ALPHA.LT.0.) ALPHAU=0.
      DO K=1,KMAX !vref
       UM(K,L)=
     *    UM(K,L)+RA(K)*(-ALPHAU*UM(K,L)+BETAU*UM(K,L+1)+DUM(K,L)) !vre
       VM(K,L)=
     *    VM(K,L)+RA(K)*(-ALPHAU*VM(K,L)+BETAU*VM(K,L+1)+DVM(K,L)) !vref
      ENDDO !vref
  380 ALPHA=BETA
C**** Subsidence uses Quadratic Upstream Scheme for QM and SM
      ML(LDMIN:LMAX) = AIRM(LDMIN:LMAX) +   DMR(LDMIN:LMAX)
      SM(LDMIN:LMAX) =   SM(LDMIN:LMAX) +  DSMR(LDMIN:LMAX)
      SMOM(:,LDMIN:LMAX) =  SMOM(:,LDMIN:LMAX) + DSMOMR(:,LDMIN:LMAX)
C****
      nsub = lmax-ldmin+1
      cmneg(ldmin:lmax)=-cm(ldmin:lmax)
      cmneg(lmax)=0. ! avoid roundoff error (esp. for qlimit)
      call adv1d(sm(ldmin),smom(1,ldmin), f(ldmin),fmom(1,ldmin),
     &     ml(ldmin),cmneg(ldmin), nsub,.false.,1, zdir)
      SM(LDMIN:LMAX) =   SM(LDMIN:LMAX) +   DSM(LDMIN:LMAX)
      SMOM(:,LDMIN:LMAX) =  SMOM(:,LDMIN:LMAX) +  DSMOM(:,LDMIN:LMAX)
C****
      ML(LDMIN:LMAX) = AIRM(LDMIN:LMAX) +   DMR(LDMIN:LMAX)
      QM(LDMIN:LMAX) =   QM(LDMIN:LMAX) +  DQMR(LDMIN:LMAX)
      QMOM(:,LDMIN:LMAX) =  QMOM(:,LDMIN:LMAX) + DQMOMR(:,LDMIN:LMAX)
      call adv1d(qm(ldmin),qmom(1,ldmin), f(ldmin),fmom(1,ldmin),
     &     ml(ldmin),cmneg(ldmin), nsub,.true.,1, zdir)
      QM(LDMIN:LMAX) =   QM(LDMIN:LMAX) +   DQM(LDMIN:LMAX)
      QMOM(:,LDMIN:LMAX) =  QMOM(:,LDMIN:LMAX) +  DQMOM(:,LDMIN:LMAX)
C**** diagnostics
      DO L=LDMIN,LMAX
        FCDH=0.
        IF(L.EQ.LMAX) FCDH=CDHSUM-(CDHSUM-CDHDRT)*.5*ETADN+CDHM
        FCDH1=0.
        IF(L.EQ.LDMIN) FCDH1=(CDHSUM-CDHDRT)*.5*ETADN-EVPSUM
        AJ8(L)=AJ8(L)+CCM(L)
        AJ13(L)=AJ13(L)+PLK(L)*(SM(L)-SMT(L))-FCDH-FCDH1
        AJ51(L)=AJ51(L)+SLHE*(QM(L)-QMT(L)+COND(L))
        AJ57(L)=AJ57(L)+SLHE*(QM(L)-QMT(L))
      END DO
c      SUMOLD=0.
c      SUMNEW=0.
c      SUMDP=0.
c      DO 382 L=LMIN,LMAX
c      SUMDP=SUMDP+AIRM(L)
c      SUMOLD=SUMOLD+SM1(L)*PLK(L)
c  382 SUMNEW=SUMNEW+SM(L)*PLK(L)
C     DIFSUM=(SUMNEW-SUMOLD)-(CDHSUM+CDHM-EVPSUM)
C     DO 383 L=LMIN,LMAX
C        AJ13(L)=AJ13(L)-DIFSUM/SUMDP
C 383 SM(L)=SM(L)-DIFSUM/SUMDP/PLK(L)
      DO 381 L=1,LM
      SM1(L)=SM(L)
  381 QM1(L)=QM(L)
C****
C**** REEVAPORATION AND PRECIPITATION
C****
      LEVAP=LMAX-1
      IF(MC1.AND.PLE(LMIN)-PLE(LMAX+1).GE.450.) LEVAP=LMIN+2
      IF(LMAX-LEVAP.GT.1) THEN
      DO 488 L=LMAX,LEVAP+2,-1
C        AJ52(L)=AJ52(L)+CDHEAT(L)
      SVWMXL(L)=0.5000000*COND(L)*BYAM(L)
  488 COND(L)=COND(L)*0.5000000
      END IF
      LEVAP=LMAX-1
C     IF(LMAX-LMIN.LE.2) GO TO 700
      PRCP=COND(LEVAP+1)
      PRHEAT=CDHEAT(LEVAP+1)
C     FEVAP=.5*FPLUME
C     IF(FEVAP.GT.1.) FEVAP=1.
         AJ50(LMAX)=AJ50(LMAX)+CDHSUM-(CDHSUM-CDHDRT)*.5*ETADN+CDHM
C        IF(LMAX.EQ.LDEP) AJ50(LMAX)=AJ50(LMAX)+CDHM
C        AJ52(LEVAP+1)=AJ52(LEVAP+1)+CDHEAT(LEVAP+1)
      DO 540 L=LMAX-1,1,-1
      IF(L.LE.LEVAP.AND.PRCP.LE.0.) GO TO 530
      FCLOUD=CCM(L)*BYAM(L+1)
      IF(PLE(LMIN)-PLE(L+2).GE.450.) FCLOUD=5.*CCM(L)*BYAM(L+1)
      IF(L.LT.LMIN) FCLOUD=CCM(LMIN)*BYAM(LMIN+1)
      IF(PLE(LMIN)-PLE(LMAX+1).LT.450.) THEN
        IF(L.EQ.LMAX-1) FCLOUD=3.*CCM(L)*BYAM(L+1)
        IF(L.LT.LMIN) FCLOUD=0.
      ENDIF
      IF(FCLOUD.GT.1.) FCLOUD=1.
      FEVAP=.5*CCM(L)*BYAM(L+1)
      IF(L.LT.LMIN) FEVAP=.5*CCM(LMIN)*BYAM(LMIN+1)
      IF(FEVAP.GT..5) FEVAP=.5
C     IF(LMAX-LMIN.LE.2) FEVAP=1.
      CLDMCL(L+1)=CLDMCL(L+1)+FCLOUD
      CLDREF=CLDMCL(L+1)
      IF(PLE(LMAX+1).GT.700..AND.CLDREF.GT.CLDSLWIJ)
     *  CLDSLWIJ=CLDREF
      IF(PLE(LMIN)-PLE(LMAX+1).GE.450..AND.CLDREF.GT.CLDDEPIJ)
     *  CLDDEPIJ=CLDREF
C**** REEVAPORATE ALL PRECIPITATION FROM ABOVE
C     EVAP=PRCP
C     PRCP=0.
C**** FORWARD STEP COMPUTES HUMIDITY CHANGE BY RECONDENSATION
C**** Q = Q + F(TOLD,PRHEAT,QOLD+EVAP)
      IF(L.GT.LEVAP) GO TO 540
      PRECNVL(L+1)=PRECNVL(L+1)+PRCP*BYGRAV
      MCLOUD=0.
      IF(L.LE.LMIN) MCLOUD=2.*FEVAP*AIRM(L)
      TOLD=SMOLD(L)*PLK(L)*BYAM(L)
      TOLD1=SMOLD(L+1)*PLK(L+1)*BYAM(L+1)
      HEAT1=0.
      IF(L.EQ.LFRZ.OR.(L.LE.LMIN.AND.TOLD.GE.TF.AND.TOLD1.LT.TF))
     *  HEAT1=LHM*PRCP*BYSHA
C     TN=TOLD-PRHEAT/MCLOUD
C     QN=QMOLD(L)*BYAM(L)+EVAP/MCLOUD
C     LHX=LHE
C     IF(TOLD.LT.TF) LHX=LHS
C     SLH=LHX*BYSHA
C     QSATC=QSAT(TN,LHX,PL(L))
C     IF(QN-QSATC.LE.0.) GO TO 520
      TNX=TOLD
      QNX=QMOLD(L)*BYAM(L)
      LHX=LHE
      IF(TNX.LT.TF) LHX=LHS
      IF(L.GT.LMIN.AND.TPSAV(L).GE.TF) LHX=LHE
      SLH=LHX*BYSHA
      DQSUM=0.
      DO 510 N=1,3
      QSATC=QSAT(TNX,LHX,PL(L))
      DQ=(QSATC-QNX)/(1.+SLH*QSATC*DQSATDT(TNX,LHX))
      TNX=TNX-SLH*DQ
      QNX=QNX+DQ
  510 DQSUM=DQSUM+DQ*MCLOUD
      IF(DQSUM.LT.0.) DQSUM=0.
      IF(DQSUM.GT.PRCP) DQSUM=PRCP
      PRCP=PRCP-DQSUM
C     IF(PRCP.GT.EVAP) PRCP=EVAP
C**** UPDATE TEMPERATURE AND HUMIDITY DUE TO NET REVAPORATION IN CLOUDS
      FSSUM = 0
      IF (ABS(PLK(L)*SM(L)).gt.1d-20) FSSUM = (SLH*DQSUM+HEAT1)/
     *     (PLK(L)*SM(L))
      SM(L)=SM(L)-(SLH*DQSUM+HEAT1)/PLK(L)
       SMOM(:,L) =  SMOM(:,L)*(1.-FSSUM)
      QM(L)=QM(L)+DQSUM
         FCDH1=0.
C        IF(L.EQ.LDEP) FCDH1=FCDH1+CDHM
         IF(L.EQ.LDMIN) FCDH1=(CDHSUM-CDHDRT)*.5*ETADN-EVPSUM
         AJ50(L)=AJ50(L)-SLH*DQSUM+FCDH1-HEAT1
         AJ52(L)=AJ52(L)-SLH*DQSUM
C**** ADD PRECIPITATION AND LATENT HEAT BELOW
  530 PRHEAT=CDHEAT(L)+SLH*PRCP
      PRCP=PRCP+COND(L)
  540 CONTINUE
C****
      IF(PRCP.GT.0.) CLDMCL(1)=CLDMCL(1)+CCM(LMIN)*BYAM(LMIN+1)
      PRCPMC=PRCPMC+PRCP
      IF(LMCMIN.GT.LDMIN) LMCMIN=LDMIN
C****
C**** END OF LOOP OVER CLOUD TYPES
C****
  570 CONTINUE
C****
C**** END OF OUTER LOOP OVER CLOUD BASE
C****
  600 CONTINUE
      IF(LMCMIN.GT.0) THEN
C**** ADJUSTMENT TO CONSERVE CP*T
        SUMAJ=0.
        SUMDP=0.
        DO L=LMCMIN,LMCMAX 
          SUMDP=SUMDP+AIRM(L)
          SUMAJ=SUMAJ+AJ13(L)
        END DO
        DO L=LMCMIN,LMCMAX
          AJ13(L)=AJ13(L)-SUMAJ*AIRM(L)/SUMDP
          SM(L)=SM(L)-SUMAJ*AIRM(L)/(SUMDP*PLK(L))
        END DO
C**** LOAD MASS EXCHANGE ARRAY FOR GWDRAG                               
        AIRXL = 0.
        DO L=LMCMIN,LMCMAX
          AIRXL = AIRXL+AJ8(L)
        END DO
      END IF
C**** CALCULATE OPTICAL THICKNESS
      WCONST=1.d-3*(WMU*(1.-PEARTH)+WMUL*PEARTH)
      WMSUM=0.
      DO L=1,LMCMAX
         TL(L)=(SM(L)*BYAM(L))*PLK(L)
         TEMWM=(TAUMCL(L)-SVWMXL(L)*AIRM(L))*1.d2*BYGRAV
         IF(TL(L).GE.TF) WMSUM=WMSUM+TEMWM
         IF (CLDMCL(L).GT.0.) TAUMCL(L)=AIRM(L)*.08d0
         IF(PLE(LMCMIN)-PLE(LMCMAX+1).LT.450..AND.(CLDMCL(L).GT.0.))THEN
            IF(L.EQ.LMCMAX) TAUMCL(L)=AIRM(L)*.02d0
         ENDIF
         IF(PLE(LMCMIN)-PLE(LMCMAX+1).GE.450..AND.(CLDMCL(L).GT.0.))THEN
            IF(L.EQ.LMCMIN) TAUMCL(L)=AIRM(L)*.02d0
         ENDIF
         IF(SVLATL(L).EQ.0.) THEN
            SVLATL(L)=LHE
            IF(SVTP(L).LT.TF) SVLATL(L)=LHS
         ENDIF
         IF(SVWMXL(L).GT.0.) THEN
            FCLD=CLDMCL(L)+1.E-20
            TEM=1.d5*SVWMXL(L)*AIRM(L)*BYGRAV
            WTEM=1.d5*SVWMXL(L)*PL(L)/(FCLD*TL(L)*RGAS)
            IF(SVLATL(L).EQ.LHE.AND.SVWMXL(L)/FCLD.GE.WCONST)
     *           WTEM=1.d5*WCONST*PL(L)/(TL(L)*RGAS)
            IF(SVLATL(L).EQ.LHS.AND.SVWMXL(L)/FCLD.GE.WMUI*1.d-3)
     *           WTEM=1.E5*WMUI*1.d-3*PL(L)/(TL(L)*RGAS)
            IF(WTEM.LT.1.d-10) WTEM=1.d-10
            IF(SVLATL(L).EQ.LHE) RCLD=(10.*(1.-PEARTH)+7.0*PEARTH)*
     *           (WTEM*4.)**BY3
            IF(SVLATL(L).EQ.LHS) RCLD=25.0*(WTEM/4.2d-3)**BY3
            RCLDE=RCLD/BYBR   !  effective droplet radius in anvil
            CSIZEL(L)=RCLDE   !  effective droplet radius in anvil
            TAUMCL(L)=1.5*TEM/(FCLD*RCLDE+1.E-20)
            IF(TAUMCL(L).GT.100.) TAUMCL(L)=100.
         END IF
      END DO

      RETURN
      END SUBROUTINE MSTCNV

      SUBROUTINE LSCOND(IERR,WMERR,LERR)
!@sum  LSCOND column physics of large scale condensation
!@auth M.S.Yao/A. Del Genio (modularisation by Gavin Schmidt)
!@ver  1.0 (taken from CB265)
!@calls CTMIX,QSAT,DQSATDT,THBAR
      IMPLICIT NONE

!@var IERR,WMERR,LERR error reporting
      INTEGER, INTENT(OUT) :: IERR,LERR
      REAL*8, INTENT(OUT) :: WMERR
      REAL*8 LHX,LHXUP

      REAL*8, PARAMETER :: CM00=1.d-4
!@param CM00 upper limit for autoconversion rate
      REAL*8, PARAMETER :: WM0=.5d-3   !@param WM0
      REAL*8, PARAMETER :: EPS=.622d0  !@param EPS

      REAL*8, DIMENSION(IM) :: UMO1,UMO2,UMN1,UMN2 !@var dummy variables
      REAL*8, DIMENSION(IM) :: VMO1,VMO2,VMN1,VMN2 !@var dummy variables
!@var Miscellaneous vertical arrays
      REAL*8, DIMENSION(LM) ::
     *     QSATL,RHF,RH1,ATH,SQ,ER,QHEAT,
     *     CAREA,PREP,RH00,EC,WMXM
!@var QSATL saturation water vapor mixing ratio
!@var RHF environmental relative humidity
!@var RH1 relative humidity to compare with the threshold humidity
!@var ATH change in potential temperature
!@var SQ dummy variables
!@var ER evaporation of precip
!@var QHEAT change of latent heat
!@var CAREA fraction of clear region
!@var PREP precip conversion rate
!@var RH00 threshold relative humidity
!@var EC cloud evaporation rate
!@var WMXM cloud water mass (mb)
      REAL*8, DIMENSION(LM+1) :: PREBAR,PREICE
!@var PREBAR,PREICE precip entering layer top for total, snow

      REAL*8 Q1,AIRMR,BETA,BMAX
     *     ,CBF,CBFC0,CK,CKIJ,CK1,CK2,CKM,CKR,CM,CM0,CM1,DFX,DQ,DQSDT
     *     ,DQSUM,DQUP,DRHDT,DSE,DSEC,DSEDIF,DWDT,DWDT1,DWM,ECRATE,EXPST
     *     ,FCLD,FMASS,FMIX,FPLUME,FPMAX,FQTOW,FRAT,FUNI,FUNIL,FUNIO
     *     ,HCHANG,HDEP,HPHASE,OLDLAT,OLDLHX,PFR,PMI,PML,HDEP1
     *     ,PRATIO,QCONV,QHEATC,QLT1,QLT2,QMN1,QMN2,QMO1,QMO2,QNEW,QNEWU
     *     ,QOLD,QOLDU,QSATC,QSATE,RANDNO,RCLDE,RHI,RHN,RHO,RHT1,RHW
     *     ,SEDGE,SIGK,SLH,SMN1,SMN2,SMO1,SMO2,TEM,TEMP,TEVAP,THT1,THT2
     *     ,TLT1,TNEW,TNEWU,TOLD,TOLDU,TOLDUP,VDEF,WCONST,WMN1,WMN2
     *     ,WMNEW,WMO1,WMO2,WMT1,WMT2,WMX1,WTEM,VVEL,XY,RCLD,FCOND
!@var Q1,BETA,BMAX,CBFC0,CKIJ,CK1,CK2 dummy variables
!@var AIRMR
!@var CBF enhancing factor for precip conversion
!@var CK ratio of cloud top jumps in moist static energy and total water
!@var CKM CTEI threshold in MacVean and Mason theory
!@var CKR CTEI threshold in Randall theory
!@var CM conversion rate for large cloud water content
!@var CM0,CM1 limiting autoconversion rate for large cloud water content
!@var DFX iteration increment
!@var DQ condensed water vapor
!@var DQSDT derivative of saturation vapor pressure w.r.t. temperature
!@var DQSUM,DQUP dummy variables
!@var DRHDT time change of relative humidity
!@var DSE moist static energy jump at cloud top
!@var DSEC critical DSE for CTEI to operate
!@var DSFDIF DSE-DSEC
!@var DWDT,DWDT1 time change rates of cloud water
!@var FQTOWump of water and vapor at cloud top
!@var ECRATE cloud droplet evaporation rate
!@var EXPST exponential term in determining the fraction in CTEI
!@var FCLD cloud fraction
!@var FMIX,FRAT fraction of mixing in CTEI
!@var FMASS mass of mixing in CTEI
!@var FPLUME fraction of mixing in CTEI
!@var FPMAX max fraction of mixing in CTEI
!@var FQTOW
!@var FUNI the probablity for ice cloud to form
!@var FUNIL,FUNIO FUNI over land, ocean
!@var HCHANG,HPHASE latent heats for changing phase
!@var HDEP,HDEP1 layer depth (Km)
!@var OLDLAT,OLDLHX previous LHX
!@var PFR PROBABLITY OF GLACIATION OF SUPER-COOLED WATER
!@var PMI icy precip entering the layer top
!@var PML layer's cloud water devided by GRAV
!@var PRATIO PMI/PML
!@var QCONV convergence of latent heat
!@var QHEATC,QLT1,QLT2,QMN1,QMN2,QMO1,QMO2 dummy variables
!@var QNEW,QNEWU updated specific humidity
!@var QOLD,QOLDU previous specific humidity
!@var QSATC saturation vapor mixing ratio
!@var QSATE saturation vapor mixing ratio w.r.t. water
!@var RANDNO random number
!@var RCLD,RCLDE cloud particle's radius, effective radius
!@var RHI relative humidity w.r.t. ice
!@var RHN dummy variable
!@var RHO air density
!@var RHT1,RHW,SIGK dummy variables
!@var SEDGE potential temperature at layer edge
!@var SMN1,SMN2,SMO1,SMO2,TEM,TEMP,TEVAP,THT1,THT2,TLT1 dummy variables
!@var TNEW,TNEWU updated tempertures
!@var TOLD,TOLDU,TOLDUP previous temperatures
!@var VDEF = VVEL - VSUB
!@var WCONST,WMN1,WMN2,WMO1,WMO2,WMT1,WMT2,WMX1 dummy variables
!@var WMNEW updated cloud water mixing ratio
!@var WTEM cloud water density (g m**-3)
!@var VVEL vertical velocity (cm/s)
!@var XY,FCOND dummy variables
      INTEGER LN,ITER !@var LN,ITER loop variables
      LOGICAL BANDF  !@var BANDF true if Bergero-Findeisen proc. occurs

      INTEGER K,L,N  !@var K,L,N loop variables

      REAL*8 THBAR !@var THBAR potential temperature at layer edge
C****
C**** LARGE-SCALE CLOUDS AND PRECIPITATION
C**** THE LIQUID WATER CONTENT IS PREDICTED
C****
      IERR=0
C****
      PRCPSS=0.
      HCNDSS=0.
      CKIJ=1.
C**** initialise vertical arrays
      ER=0.
      EC=0.
      PREP=0.
      PREBAR=0.
      QHEAT=0.
      CLDSSL=0
      TAUSSL=0
      DO L=1,LM
         CAREA(L)=1.-CLDSAVL(L)
         IF(WMX(L).LE.0.) CAREA(L)=1.
C     CAREA(L)=(1.-RH(L))/(1.-RHF(L)*0.999+1.E-20)
C     IF(CAREA(L).GT.1.) CAREA(L)=1.
C     IF(RH(L).GT.1.) CAREA(L)=0.
      END DO
      DQUP=0.
      LHXUP=LHE
      TOLDUP=TL(LM)
      PREICE(LM+1)=0.
      WCONST=WMU*(1.-PEARTH)+WMUL*PEARTH
         AJ11=0.
         AJ55=0.
C****
C**** MAIN L LOOP FOR LARGE-SCALE CONDENSATION, PRECIPITATION AND CLOUDS
C****
      DO 304 L=LM,1,-1
      TOLD=TL(L)
      QOLD=QL(L)
      OLDLHX=SVLHXL(L)
      OLDLAT=SVLATL(L)
C**** COMPUTE VERTICAL VELOCITY IN CM/S
      TEMP=100.*RGAS*TL(L)/(PL(L)*GRAV)
      IF(L.EQ.1) VVEL=-SDL(L+1)*TEMP
      IF(L.EQ.LM) VVEL=-SDL(L)*TEMP
      IF(L.GT.1.AND.L.LT.LM)
     *     VVEL=-.5*(SDL(L)+SDL(L+1))*TEMP
C**** COMPUTE THE LIMITING AUTOCONVERSION RATE FOR CLOUD WATER CONTENT
      CM0=CM00
      VDEF=VVEL-VSUBL(L)
      IF(VDEF.GT.0.) CM0=CM00*10.**(-VDEF)
      FCLD=1.-CAREA(L)+1.E-20
C**** COMPUTE THE PROBABILITY OF ICE FORMATION, FUNI, AND
C**** THE PROBABLITY OF GLACIATION OF SUPER-COOLED WATER, PFR
C**** DETERMINE THE PHASE MOISTURE CONDENSES TO
C**** DETERMINE THE POSSIBILITY OF B-F PROCESS
      BANDF=.FALSE.
      LHX= LHE
      PMI=PREICE(L+1)*DTsrc
      PML=WMX(L)*AIRM(L)*BYGRAV
      PRATIO=PMI/(PML+1.E-20)
      IF(PRATIO.GT.10.) PRATIO=10.
      CBF=1.+1.*EXP(-((TL(L)-258.16d0)/10.)**2)
      CBFC0=.5*CM0*CBF*DTsrc
      PFR=(1.-EXP(-(PRATIO*PRATIO)))*(1.-EXP(-(CBFC0*CBFC0)))
      FUNIO=1.-EXP(-((TL(L)-269.16d0)/15.)**2)
      FUNIL=1.-EXP(-((TL(L)-263.16d0)/15.)**2)
      IF(TL(L).GT.269.16) FUNIO=0.
      IF(TL(L).GT.263.16) FUNIL=0.
      FUNI=FUNIO*(1.-PEARTH)+FUNIL*PEARTH
      IF(TL(L).LE.TI) FUNI=1.
      IF(TL(L).LT.TI) LHX=LHS
      IF(TL(L).GE.TI.AND.RANDU(XY).LT.FUNI) LHX=LHS
      IF((OLDLHX.EQ.LHS.OR.OLDLAT.EQ.LHS).AND.TL(L).LT.TF) THEN
        IF(LHX.EQ.LHE) BANDF=.TRUE.
        LHX=LHS
      ENDIF
      IF(L.LT.LM) THEN
      RANDNO=RANDU(XY)
      IF(PFR.GT.RANDNO.AND.TL(L).LT.TF) THEN
        IF(LHX.EQ.LHE) BANDF=.TRUE.
        LHX=LHS
      ENDIF
      IF(LHXUP.EQ.LHE) LHX=LHE
      IF(LHXUP.EQ.LHE) BANDF=.FALSE.
      END IF
C**** COMPUTE RELATIVE HUMIDITY
      QSATL(L)=QSAT(TL(L),LHX,PL(L))
      RH1(L)=QL(L)/QSATL(L)
      IF(LHX.EQ.LHS) THEN
      QSATE=QSAT(TL(L),LHE,PL(L))
      RHW=.00536d0*TL(L)-.276d0
      RH1(L)=QL(L)/QSATE
      IF(TL(L).LT.238.16) RH1(L)=QL(L)/(QSATE*RHW)
      END IF
C**** PHASE CHANGE OF CLOUD WATER CONTENT
      HCHANG=0.
      IF(LHX.EQ.LHS) THEN
        IF(OLDLHX.EQ.LHE) HCHANG=WML(L)*LHM
        IF(OLDLHX.EQ.LHE.OR.OLDLAT.EQ.LHE) BANDF=.TRUE.
      END IF
      IF(OLDLHX.EQ.LHS.AND.LHX.EQ.LHE) HCHANG=-WML(L)*LHM
      IF(OLDLAT.EQ.LHE.AND.LHX.EQ.LHS) HCHANG=HCHANG+SVWMXL(L)*LHM
      IF(OLDLAT.EQ.LHS.AND.LHX.EQ.LHE) HCHANG=HCHANG-SVWMXL(L)*LHM
      SVLHXL(L)=LHX
      TL(L)=TL(L)+HCHANG/SHA
      TH(L)=TL(L)/PLK(L)
      ATH(L)=(TH(L)-TTOLDL(L))*BYDTsrc
C**** COMPUTE RH IN THE CLOUD-FREE AREA, RHF
      RHI=QL(L)/QSAT(TL(L),LHS,PL(L))
      RH00(L)=U00wtr
      IF(LHX.EQ.LHS) RH00(L)=U00ice
      IF(L.EQ.1) THEN
        HDEP=AIRM(L)*TL(L)*RGAS/(1000.*GRAV*PL(L))
        RH00(L)=1.-9.8d0*LHE*HDEP/(RVAP*TS*TS)
        IF(DCL.LE.1) THEN
         IF(RI1.LT..25) HDEP1=.5*HDEP
         IF(RI1.GE..25.AND.RI1.LT.1.) HDEP1=.01-(.5*HDEP-.01)*
     *     (RI1-1.)/.75
         IF(RI1.GE.1.) HDEP1=.01
         IF(RI2.LT..25) HDEP1=HDEP
         IF(RI2.GE..25.AND.RI2.LT.1.) HDEP1=.5*HDEP-.5*HDEP*(RI2-1.)/.75
         IF(RI2.GE.1.) HDEP1=.5*HDEP
         RH00(L)=1.-9.8*LHE*HDEP1/(RVAP*TS*TS)
        ENDIF
        IF(RH00(L).LT.0.) RH00(L)=0.
      ENDIF
      IF(L.GT.1.AND.PLE(L+1).GT.930.) THEN
        HDEP=0.
        DO 216 LN=L,1,-1
  216   HDEP=HDEP+AIRM(LN)*TL(LN)*RGAS/(1000.*GRAV*PL(LN))
        RH00(L)=1.-9.8d0*LHE*HDEP/(RVAP*TS*TS)
        IF(RH00(L).LT.0.) RH00(L)=0.
      ENDIF
C     IF(L.LE.LPBL) RH00(L)=.75
      IF(RH00(L).GT.1.) RH00(L)=1.
      RHF(L)=RH00(L)+(1.-CAREA(L))*(1.-RH00(L))
      IF(WMX(L).GT.0.) THEN
C**** COMPUTE THE AUTOCONVERSION RATE OF CLOUD WATER TO PRECIPITATION
      RHO=1.E5*PL(L)/(RGAS*TL(L))
      TEM=RHO*WMX(L)/(WCONST*FCLD+1.E-20)
      IF(LHX.EQ.LHS ) TEM=RHO*WMX(L)/(WMUI*FCLD+1.E-20)
      TEM=TEM*TEM
      IF(TEM.GT.10.) TEM=10.
      CM1=CM0
      IF(BANDF) CM1=CM0*CBF
      IF(LHX.EQ.LHS) CM1=CM0
      CM=CM1*(1.-1./EXP(TEM*TEM))+1.*100.*(PREBAR(L+1)+
     *   PRECNVL(L+1)*BYDTsrc)
      IF(CM.GT.BYDTsrc) CM=BYDTsrc
      PREP(L)=WMX(L)*CM
      END IF
C**** FORM CLOUDS ONLY IF RH GT RH00
  219 IF(RH1(L).LT.RH00(L)) GO TO 220
C**** COMPUTE THE CONVERGENCE OF AVAILABLE LATENT HEAT
      SQ(L)=LHX*QSATL(L)*DQSATDT(TL(L),LHX)*BYSHA
      TEM=-LHX*DPDT(L)/PL(L)
      QCONV=LHX*AQ(L)-RH(L)*SQ(L)*SHA*PLK(L)*ATH(L)
     *  -TEM*QSATL(L)*RH(L)
      IF(QCONV.LE.0.0.AND.WMX(L).LE.0.) GO TO 220
C**** COMPUTE EVAPORATION OF RAIN WATER, ER
      RHN=RHF(L)
      IF(RHF(L).GT.RH(L)) RHN=RH(L)
C     QCONV0=LHX*AQ(L)-RHN*SQ(L)*SHA*PLK(L)*ATH(L)
C    *  -TEM*QSATL(L)*RHN
      IF(WMX(L).GT.0.) ER(L)=(1.-RHN)*LHX*PREBAR(L+1)*GRAV*BYAM(L)
      IF(WMX(L).LE.0.) ER(L)=(1.-RH(L))*LHX*PREBAR(L+1)*GRAV*BYAM(L)
      IF(WMX(L).LE.0..AND.PREICE(L+1).GT.0..AND.TL(L).LT.TF)
     *  ER(L)=(1.-RHI)*LHX*PREBAR(L+1)*GRAV*BYAM(L)
      IF(ER(L).LT.0.) ER(L)=0.
C**** COMPUTATION OF CLOUD WATER EVAPORATION
      IF (CAREA(L).GT.0.) THEN
      WTEM=1d5*WMX(L)*PL(L)/(FCLD*TL(L)*RGAS+1.E-20)
      IF(LHX.EQ.LHE.AND.WMX(L)/FCLD.GE.WCONST*1.E-3)
     *  WTEM=1d5*WCONST*1d-3*PL(L)/(TL(L)*RGAS)
      IF(LHX.EQ.LHS.AND.WMX(L)/FCLD.GE.WMUI*1.E-3)
     *  WTEM=1d5*WMUI*1d-3*PL(L)/(TL(L)*RGAS)
      IF(WTEM.LT.1.E-10) WTEM=1d-10
      IF(LHX.EQ.LHE) RCLD=1d-6*(10.*(1.-PEARTH)+7.0*PEARTH)*
     *  (WTEM*4.)**BY3
      IF(LHX.EQ.LHS) RCLD=25.d-6*(WTEM/4.2d-3)**BY3
      CK1=1000.*LHX*LHX/(2.4d-2*RVAP*TL(L)*TL(L))
      CK2=1000.*RVAP*TL(L)/(2.4d-3*QSATL(L)*PL(L)/.622d0)
c      CK2=1000.*RGAS*TL(L)/(2.4d-3*QSATL(L)*PL(L))    ! new
      TEVAP=1000.*(CK1+CK2)*RCLD*RCLD
      WMX1=WMX(L)-PREP(L)*DTsrc
      ECRATE=(1.-RHF(L))/(TEVAP*FCLD+1.E-20)
      IF(ECRATE.GT.BYDTsrc) ECRATE=BYDTsrc
      EC(L)=WMX1*ECRATE*LHX
      END IF
C**** COMPUTE NET LATENT HEATING DUE TO STRATIFORM CLOUD PHASE CHANGE,
C**** QHEAT, AND NEW CLOUD WATER CONTENT, WMNEW
      DRHDT=2.*CAREA(L)*CAREA(L)*(1.-RH00(L))*(QCONV+ER(L)+EC(L))/LHX/
     *  (WMX(L)/(FCLD+1.E-20)+2.*CAREA(L)*QSATL(L)*(1.-RH00(L))+1.E-20)
      IF(ER(L).EQ.0..AND.EC(L).EQ.0..AND.WMX(L).LE.0.) DRHDT=0.
C     IF(RH1(L).GT.1.) DRHDT=0.
      QHEAT(L)=(QCONV-LHX*DRHDT*QSATL(L))/(1.+RH(L)*SQ(L))
      DWDT=QHEAT(L)/LHX-PREP(L)+CAREA(L)*ER(L)/LHX
      WMNEW =WMX(L)+DWDT*DTsrc
      IF(WMNEW.LT.0.) THEN
      WMNEW=0.
      QHEAT(L)=(-WMX(L)*BYDTsrc+PREP(L))*LHX-CAREA(L)*ER(L)
      END IF
      GO TO 230
C**** UNFAVORABLE CONDITIONS FOR CLOUDS TO EXIT, PRECIP OUT CLOUD WATER
  220 Q1=0.
      IF(WMX(L).GT.0.) PREP(L)=WMX(L)*BYDTsrc
      ER(L)=(1.-RH(L))*LHX*PREBAR(L+1)*GRAV*BYAM(L)
      IF(PREICE(L+1).GT.0..AND.TL(L).LT.TF)
     *  ER(L)=(1.-RHI)*LHX*PREBAR(L+1)*GRAV*BYAM(L)
      IF(ER(L).LT.0.) ER(L)=0.
      QHEAT(L)=-CAREA(L)*ER(L)+Q1
      WMNEW=0.
  230 CONTINUE
C**** PHASE CHANGE OF PRECIPITATION, FROM ICE TO WATER
      HPHASE=0.
      IF(L.LT.LM.AND.TL(L).GT.TF)
     *  HPHASE=LHM*PREICE(L+1)*GRAV*BYAM(L)
C**** COMPUTE THE PRECIP AMOUNT ENTERING THE LAYER TOP
      IF(TL(L).GT.TF) PREICE(L+1)=0.
      PREICE(L)=PREICE(L+1)-AIRM(L)*ER(L)*CAREA(L)*PREICE(L+1)/
     *          (GRAV*LHX*PREBAR(L+1)+1.E-20)
      IF(LHX.EQ.LHS) PREICE(L)=PREICE(L)+AIRM(L)*PREP(L)*BYGRAV
      PREBAR(L)=PREBAR(L+1)+
     *          AIRM(L)*(PREP(L)-ER(L)*CAREA(L)/LHX)*BYGRAV
C**** UPDATE NEW TEMPERATURE AND SPECIFIC HUMIDITY
      QNEW =QL(L)-DTsrc*QHEAT(L)/LHX
      IF(QNEW.LT.0.) THEN
      QNEW=0.
      QHEAT(L)=QL(L)*LHX*BYDTsrc
      DWDT1=QHEAT(L)/LHX-PREP(L)+CAREA(L)*ER(L)/LHX
      WMNEW=WMX(L)+DWDT1*DTsrc
C**** IF WMNEW .LT. 0., THE COMPUTATION IS UNSTABLE
      IF(WMNEW.LT.0.) THEN
        IERR=1
        LERR=L
        WMERR=WMNEW
        WMNEW=0.
      END IF
      END IF
C**** Only Calculate fractional changes of Q to W
c      FPR=0.
c      IF (WMX(L).gt.1d-20) FPR=PREP(L)*DTsrc/WMX(L)          ! CLW->P
c      FER=0.
c      IF (PREBAR(L+1).gt.1d-20) FER=CAREA(L)*ER(L)*AIRM(L)/(
c     *     GRAV*LHX*PREBAR(L+1))                              ! P->Q
      FQTOW=0.                                                ! Q->CLW
c      FWTOQ=0.                                                ! CLW->Q
      IF (QHEAT(L)+CAREA(L)*ER(L).gt.0) THEN
        IF (LHX*QL(L)+DTsrc*CAREA(L)*ER(L).gt.1d-20) FQTOW=(QHEAT(L)
     *       +CAREA(L)*ER(L))*DTsrc/(LHX*QL(L)+DTsrc*CAREA(L)*ER(L))
c      ELSE
c        IF (WMX(L)-PREP(L)*DTsrc.gt.1d-20) FWTOQ=-(QHEAT(L)+CAREA(L)
c     *       *ER(L))*DTsrc/(LHX*(WMX(L)-PREP(L)*DTsrc))
      END IF
      QL(L)=QNEW
C**** adjust gradients down if Q decreases
       QMOM(:,L)= QMOM(:,L)*(1.-FQTOW)
      WMX(L)=WMNEW
      TL(L)=TL(L)+DTsrc*(QHEAT(L)-HPHASE)/SHA
      TH(L)=TL(L)/PLK(L)
      TNEW=TL(L)
      QSATC=QSAT(TL(L),LHX,PL(L))
      RH(L)=QL(L)/QSATC
C**** CONDENSE MORE MOISTURE IF RELATIVE HUMIDITY .GT. 1
      IF(RH(L).GT.1.) THEN
      SLH=LHX*BYSHA
      DQSUM=0.
      DO 231 N=1,3
      IF(N.NE.1) QSATC=QSAT(TL(L),LHX,PL(L))
      DQ=(QL(L)-QSATC)/(1.+SLH*QSATC*DQSATDT(TL(L),LHX))
      TL(L)=TL(L)+SLH*DQ
      QL(L)=QL(L)-DQ
  231 DQSUM=DQSUM+DQ
      IF(DQSUM.GT.0.) THEN
      WMX(L)=WMX(L)+DQSUM
      FCOND=DQSUM/QNEW
C**** adjust gradients down if Q decreases
       QMOM(:,L)= QMOM(:,L)*(1.-FCOND)
      ELSE
      TL(L)=TNEW
      QL(L)=QNEW
      END IF
      RH(L)=QL(L)/QSAT(TL(L),LHX,PL(L))
      TH(L)=TL(L)/PLK(L)
      TNEW=TL(L)
      END IF
      IF(RH(L).LE.RHF(L)) THEN
C**** PRECIP OUT CLOUD WATER IF RH LESS THAN THE RH OF THE ENVIRONMENT
      PREBAR(L)=PREBAR(L)+WMX(L)*AIRM(L)*BYGRAV*BYDTsrc
      WMX(L)=0.
      END IF
C**** COMPUTE THE LARGE-SCALE CLOUD COVER
      IF(RH(L).LE.1.) CAREA(L)=DSQRT((1.-RH(L))/(1.-RH00(L)+1.E-20))
      IF(CAREA(L).GT.1.) CAREA(L)=1.
      IF(RH(L).GT.1.) CAREA(L)=0.
      IF(WMX(L).LE.0.) CAREA(L)=1.
      IF(CAREA(L).LT.0.) CAREA(L)=0.
      CLDSSL(L)=1.-CAREA(L)
      TOLDUP=TOLD
      LHXUP=LHX
C**** ACCUMULATE SOME DIAGNOSTICS
         HCNDSS=HCNDSS+(TNEW-TOLD)*AIRM(L)
  304    AJ11(L)=AJ11(L)+(TNEW-TOLD)*AIRM(L)
C****
C**** CLOUD-TOP ENTRAINMENT INSTABILITY
C****
C     DO 310 L=1,LM
      DO 382 L=LM-1,1,-1
      SM(L)=TH(L)*AIRM(L)
      QM(L)=QL(L)*AIRM(L)
      WMXM(L)=WMX(L)*AIRM(L)
      SM(L+1)=TH(L+1)*AIRM(L+1)
      QM(L+1)=QL(L+1)*AIRM(L+1)
      WMXM(L+1)=WMX(L+1)*AIRM(L+1)
      TOLD=TL(L)
      TOLDU=TL(L+1)
         QOLD=QL(L)
         QOLDU=QL(L+1)
      FCLD=1.-CAREA(L)+1.E-30
      IF(CAREA(L).EQ.1.) GO TO 382
      IF(CAREA(L).LT.1..AND.CAREA(L+1).LT.1.) GO TO 382
      SEDGE=THBAR(TH(L+1),TH(L))
      DSE=(TH(L+1)-SEDGE)*PLK(L+1)+(SEDGE-TH(L))*PLK(L)+
     *  SLHE*(QL(L+1)-QL(L))
      DWM=QL(L+1)-QL(L)+(WMX(L+1)-WMX(L))/FCLD
      DQSDT=DQSATDT(TL(L),LHE)*QL(L)/(RH(L)+1d-30)
      BETA=(1.+BYMRAT*TL(L)*DQSDT)/(1.+SLHE*DQSDT)
      CKM=(1.+SLHE*DQSDT)*(1.+.392d0*TL(L)/SLHE)/
     *  (2.+(1.+BYMRAT*TL(L)/SLHE)*SLHE*DQSDT)
      CKR=TL(L)/(BETA*SLHE)
      CK=DSE/(SLHE*DWM)
      SIGK=0.
      IF(CKR.GT.CKM) GO TO 382
      IF(CK.GT.CKR) SIGK=1d-3*((CK-CKR)/(CKM-CKR+1.E-20))**5.
      EXPST=EXP(-SIGK*DTsrc)
      IF(L.LE.1) CKIJ=EXPST
      DSEC=DWM*TL(L)/BETA
C     DSEC=.53*SLHE*DWM
C     TIME0=ABS(8./(DSE-DSEC+1.E-30))
C     FPMAX=1.-EXP(-1./TIME0)
C     FPMAX=1.
      IF(CK.LT.CKR) GO TO 382
      FPMAX=1.-EXPST
C     IF(CK.GT.CKM) FPMAX=1.
      IF(FPMAX.LE.0.) GO TO 382
      IF(FPMAX.GT.1.) FPMAX=1.
      IF(DSE.GE.DSEC) GO TO 382
C**** MIXING TO REMOVE CLOUD-TOP ENTRAINMENT INSTABILITY
      AIRMR=(AIRM(L+1)+AIRM(L))*BYAM(L+1)*BYAM(L)
      SMO1=SM(L)
      QMO1=QM(L)
      WMO1=WMXM(L)
      SMO2=SM(L+1)
      QMO2=QM(L+1)
      WMO2=WMXM(L+1)
      DO K=1,KMAX !vref
         UMO1(K)=UM(K,L) !vref
         VMO1(K)=VM(K,L) !vref
         UMO2(K)=UM(K,L+1) !vref
         VMO2(K)=VM(K,L+1) !vref
      ENDDO !vref
      FPLUME=FPMAX
      DFX=FPMAX
      DO 320 ITER=1,9
      DFX=DFX*0.5
      FMIX=FPLUME*FCLD
      FMASS=FMIX*AIRM(L)
      FMASS=MIN(FMASS,(AIRM(L+1)*AIRM(L))/(AIRM(L+1)+AIRM(L)))
      FMIX=FMASS*BYAM(L)
      FRAT=FMASS*BYAM(L+1)
      SMN1=SMO1*(1.-FMIX)+FRAT*SMO2
      QMN1=QMO1*(1.-FMIX)+FRAT*QMO2
      WMN1=WMO1*(1.-FMIX)+FRAT*WMO2
      SMN2=SMO2*(1.-FRAT)+FMIX*SMO1
      QMN2=QMO2*(1.-FRAT)+FMIX*QMO1
      WMN2=WMO2*(1.-FRAT)+FMIX*WMO1
      THT1=SMN1*BYAM(L)
      QLT1=QMN1*BYAM(L)
      TLT1=THT1*PLK(L)
      LHX=SVLHXL(L)
      RHT1=QLT1/(QSAT(TLT1,LHX,PL(L)))
      WMT1=WMN1*BYAM(L)
      THT2=SMN2*BYAM(L+1)
      QLT2=QMN2*BYAM(L+1)
      WMT2=WMN2*BYAM(L+1)
      SEDGE=THBAR(THT2,THT1)
      DSE=(THT2-SEDGE)*PLK(L+1)+(SEDGE-THT1)*PLK(L)+SLHE*(QLT2-QLT1)
      DWM=QLT2-QLT1+(WMT2-WMT1)/FCLD
      DQSDT=DQSATDT(TLT1,LHE)*QLT1/(RHT1+1d-30)
      BETA=(1.+BYMRAT*TLT1*DQSDT)/(1.+SLHE*DQSDT)
      CKM=(1.+SLHE*DQSDT)*(1.+.392d0*TLT1/SLHE)/
     *  (2.+(1.+BYMRAT*TLT1/SLHE)*SLHE*DQSDT)
      DSEC=DWM*TLT1/BETA
C     DSEC=.53*SLHE*DWM
      DSEDIF=DSE-DSEC
      IF(DSEDIF.GT.1.E-3) FPLUME=FPLUME-DFX
      IF(DSEDIF.LT.-1.E-3) FPLUME=FPLUME+DFX
      IF(ABS(DSEDIF).LE.1.E-3.OR.FPLUME.GT.FPMAX) GO TO 380
  320 CONTINUE
C**** UPDATE TEMPERATURE, SPECIFIC HUMIDITY AND MOMENTUM DUE TO CTEI
  380 TH(L)=SMN1*BYAM(L)
      TL(L)=TH(L)*PLK(L)
      QL(L)=QMN1*BYAM(L)
      LHX=SVLHXL(L)
      RH(L)=QL(L)/QSAT(TL(L),LHX,PL(L))
      WMX(L)=WMN1*BYAM(L)
      TH(L+1)=SMN2*BYAM(L+1)
      QL(L+1)=QMN2*BYAM(L+1)
      WMX(L+1)=WMN2*BYAM(L+1)
      CALL CTMIX (SM(L),SMOM(1,L),FMASS*AIRMR,FMIX,FRAT)
      CALL CTMIX (QM(L),QMOM(1,L),FMASS*AIRMR,FMIX,FRAT)
      DO K=1,KMAX !vref
         UMN1(K)=(UMO1(K)*(1.-FMIX)+FRAT*UMO2(K)) !vref
         VMN1(K)=(VMO1(K)*(1.-FMIX)+FRAT*VMO2(K)) !vref
         UMN2(K)=(UMO2(K)*(1.-FRAT)+FMIX*UMO1(K)) !vref
         VMN2(K)=(VMO2(K)*(1.-FRAT)+FMIX*VMO1(K)) !vref
         UM(K,L)=UM(K,L)+(UMN1(K)-UMO1(K))*RA(K) !vref
         VM(K,L)=VM(K,L)+(VMN1(K)-VMO1(K))*RA(K) !vref
         UM(K,L+1)=UM(K,L+1)+(UMN2(K)-UMO2(K))*RA(K) !vref
         VM(K,L+1)=VM(K,L+1)+(VMN2(K)-VMO2(K))*RA(K) !vref
      ENDDO !vref
         QNEW=QL(L)
         QNEWU=QL(L+1)
C**** RE-EVAPORATION OF LWC IN THE UPPER LAYER
      QL(L+1)=QL(L+1)+WMX(L+1)
      TH(L+1)=TH(L+1)-(LHX*BYSHA)*WMX(L+1)/PLK(L+1)
      TL(L+1)=TH(L+1)*PLK(L+1)
      RH(L+1)=QL(L+1)/QSAT(TL(L+1),LHX,PL(L+1))
      WMX(L+1)=0.
      IF(RH(L).LE.1.) CAREA(L)=DSQRT((1.-RH(L))/(1.-RH00(L)+1.E-20))
      IF(CAREA(L).GT.1.) CAREA(L)=1.
      IF(RH(L).GT.1.) CAREA(L)=0.
      CLDSSL(L)=1.-CAREA(L)
      TNEW=TL(L)
      TNEWU=TL(L+1)
         HCNDSS=HCNDSS+(TNEW-TOLD)*AIRM(L)+(TNEWU-TOLDU)*AIRM(L+1)
         AJ11(L)=AJ11(L)+(TNEW-TOLD)*AIRM(L)
         AJ11(L+1)=AJ11(L+1)+(TNEWU-TOLDU)*AIRM(L+1)
         AJ55(L)=AJ55(L)+(QNEW-QOLD)*AIRM(L)*LHX*BYSHA
         AJ55(L+1)=AJ55(L+1)+(QNEWU-QOLDU)*AIRM(L+1)*LHX*BYSHA
  382 CONTINUE
         WMSUM=0.
C**** COMPUTE CLOUD PARTICLE SIZE AND OPTICAL THICKNESS
      DO 388 L=1,LM
      FCLD=CLDSSL(L)+1.E-20
      WTEM=1.d5*WMX(L)*PL(L)/(FCLD*TL(L)*RGAS+1.E-20)
      LHX=SVLHXL(L)
C     IF(LHX.EQ.LHE.AND.WMX(L)/FCLD.GE.WCONST*1.E-3)
C    *  WTEM=1.E5*WCONST*1.E-3*PL(L)/(TL(L)*RGAS)
      IF(LHX.EQ.LHS.AND.WMX(L)/FCLD.GE.WMUI*1.E-3)
     *  WTEM=1d5*WMUI*1.d-3*PL(L)/(TL(L)*RGAS)
      IF(WTEM.LT.1.E-10) WTEM=1.d-10
      IF(LHX.EQ.LHE) RCLD=(10.*(1.-PEARTH)+7.0*PEARTH)*
     *  (WTEM*4.)**BY3
      IF(LHX.EQ.LHE) THEN
         QHEATC=(QHEAT(L)+CAREA(L)*(EC(L)+ER(L)))/LHX
         IF(RCLD.GT.20..AND.PREP(L).GT.QHEATC) RCLD=20.
      ENDIF
      IF(LHX.EQ.LHS) RCLD=25.0*(WTEM/4.2d-3)**BY3
      RCLDE=RCLD/BYBR
      CSIZEL(L)=RCLDE
      TEM=AIRM(L)*WMX(L)*1.d2*BYGRAV
      TAUSSL(L)=1.5d3*TEM/(FCLD*RCLDE+1.E-20)
      IF(TAUSSL(L).GT.100.) TAUSSL(L)=100.
  388    IF(LHX.EQ.LHE) WMSUM=WMSUM+TEM
      PRCPSS=PREBAR(1)*GRAV*DTsrc

C**** CALCULATE OPTICAL THICKNESS
      DO L=1,LM
      CLDSAVL(L)=CLDSSL(L)
      IF(TAUMCL(L).GT.0..AND.CKIJ.EQ.1.) GO TO 526
      BMAX=1.-EXP(-(CLDSAVL(L)/.3d0))
      IF(CLDSAVL(L).GE..95) BMAX=CLDSAVL(L)
      IF(L.EQ.1.OR.PLE(L+1).GT.930.) THEN
        CLDSSL(L)=CLDSSL(L)+
     *               (BMAX-CLDSSL(L))*CKIJ
        TAUSSL(L)=TAUSSL(L)*CLDSAVL(L)/
     *               (CLDSSL(L)+1.E-20)
      ENDIF
      IF(TAUSSL(L).LE.0.) CLDSSL(L)=0.
      IF(L.EQ.1.OR.PLE(L+1).GT.930..OR.TAUMCL(L).GT.0.) GO TO 526
      CLDSSL(L)=CLDSSL(L)**(2.*BY3)
      TAUSSL(L)=TAUSSL(L)*CLDSAVL(L)**BY3
  526 CONTINUE
      IF(WMX(L).LE.0.) SVLHXL(L)=0.
      END DO

      RETURN
      END SUBROUTINE LSCOND

      END MODULE CLD01
