#include "rundeck_opts.h"

      MODULE tracers_DRYDEP
c
!@sum  tracers_DRYDEP tracer dry deposition from Harvard CTM.
!@+    Current version only calculates the "bulk surface resistance
!@+    to deposition" component of the deposition velocity.
!@auth D.J. Jacob and Y.H. Wang, modularized by G.M. Gardner, 
!@+    adapted for GISS GCM by D. Koch, modelEified by G. Faluvegi
!@ver  1.0 (based on DRYDEP subroutines in DB396Tds3M23.f, based on 
!@+    Harvard version 3.1: 12/17/97)  
C*********************************************************************
C  Literature cited in drydep.f routines: 
C     Baldocchi, D.D., B.B. Hicks, and P. Camara, A canopy stomatal
C       resistance model for gaseous deposition to vegetated surfaces,
C       Atmos. Environ. 21, 91-101, 1987.
C     Guenther, A., and 15 others, A global model of natural volatile
C       organic compound emissions, J. Geophys. Res., 100, 8873-8892,
C       1995.  
C     Brutsaert, W., Evaporation into the Atmosphere, Reidel, 1982.
C     Dwight, H.B., Tables of integrals and other mathematical data,
C       MacMillan, 1957.
C     Hicks, B.B., and P.S. Liss, Transfer of SO2 and other reactive
C       gases across the air-sea interface, Tellus, 28, 348-354, 1976.
C     Jacob, D.J., and S.C. Wofsy, Budgets of reactive nitrogen,
C       hydrocarbons, and ozone over the Amazon forest during the wet
C       season, J.  Geophys. Res., 95, 16737-16754, 1990.
C     Jacob, D.J., and 9 others, Deposition of ozone to tundra,
C       J. Geophys. Res., 97, 16473-16479, 1992.
C     Levine, I.N., Physical Chemistry, 3rd ed., McGraw-Hill, New York,
C        1988.  
C     Munger, J.W., and 8 others, Atmospheric deposition of reactive
C       nitrogen oxides and ozone in a temperate deciduous forest and a
C       sub-arctic woodland, J. Geophys. Res., in press, 1996.
C     Walcek, C.J., R.A. Brost, J.S. Chang, and M.L. Wesely, SO2,
C        sulfate, and HNO3 deposition velocities computed using
C        regional landuse and meteorological data, Atmos. Environ.,
C        20, 949-964, 1986.   
C     Wang, Y.H., paper in preparation, 1996.   
C     Wesely, M.L, Improved parameterizations for surface resistance to
C       gaseous dry deposition in regional-scale numerical models,
C       Environmental Protection Agency Report EPA/600/3-88/025,   
C       Research Triangle Park (NC), 1988.   
C     Wesely, M.L., same title, Atmos. Environ., 23, 1293-1304, 1989.
C*********************************************************************
      USE TRACER_COM, only   : ntm
      USE MODEL_COM, only   : im,jm
c
      IMPLICIT NONE
      SAVE
C
!@param NPOLY ?
!@param NVEGTYPE number of Olson vegetation types
!@param NTYPE number of surface types
!@var IJREG # of landtypes in grid square
!@var IJLAND Land type ID for element LDT =1, IJREG(I,J) (could be
!@+          from any source - mapped to deposition surface ID in
!@+          input files LAI## where ## is a 2-digit month number.)
!@var IJUSE Fraction ((per mil) of gridbox area occupied by land type
!@+          element LDT
!@var IREG ?
!@var IDEP ?
!@var IRI read in ___ resistance
!@var IRLU Cuticular resistances per unit area of leaf
!@var IRAC read in ___ resistance
!@var IRGSS read in ___ resistance
!@var IRGSO read in ___ resistance
!@var IRCLS read in ___ resistance
!@var IRCLO read in ___ resistance
!@var IVSMAX maximum deposition velocity for aerosol from file
!@var XYLAI Leaf Area Index of land type element LDT
!@var XLAI leaf area index variable 
!@var XLAI2 leaf area index variable
!@var DRYCOEFF polynomial fittings coeffcients  
!@var CZ Altitude (m) at which deposition velocity would be computed
!@var dtr_dd to save drydep change for conservation quantities
c
      INTEGER, PARAMETER :: NPOLY   = 20,
     &                      NTYPE   = 16,
     &                      NVEGTYPE= 74
      REAL*8,  DIMENSION(IM,JM,NTYPE)     :: XYLAI,XLAI,XLAI2
      REAL*8,  DIMENSION(NPOLY)           :: DRYCOEFF
      REAL*8,  DIMENSION(JM,ntm)          :: dtr_dd
      INTEGER, DIMENSION(IM,JM)           :: IJREG,IREG
      INTEGER, DIMENSION(IM,JM,NTYPE)     :: IJLAND,IJUSE
      INTEGER, DIMENSION(NVEGTYPE)        :: IDEP
      INTEGER, DIMENSION(NTYPE)           :: IRI,IRLU,IRAC,IRGSS,
     &                                       IRGSO,IRCLS,IRCLO,IVSMAX
       
      END MODULE tracers_DRYDEP
C
C      
#ifdef TRACERS_DRYDEP
      SUBROUTINE get_dep_vel(I,J,ITYPE,OBK,ZHH,USTARR,TEMPK,DEP_VEL)
!@sum  get_dep_vel computes the Bulk surface reistance to
!@+    tracer dry deposition using a resistance-in-series model
!@+    from a portion of the Harvard CTM dry deposition routine.
!@+    The deposition velocity is the reciprocal of the
!@+    Bulk surface reistance...
!@auth D.J. Jacob and Y.H. Wang, modularized by G.M. Gardner, 
!@+    adapted for GISS GCM by D. Koch modelEified by G. Faluvegi
!@ver  1.0 (based on DRYDEP subroutines in DB396Tds3M23.f, based on
!@+    Harvard version 3.1: 12/17/97)  
C      uses functions: BIOFIT,DIFFG
c
C**** GLOBAL parameters and variables:  
C
      USE MODEL_COM,  only : im,jm
      USE GEOM,       only : imaxj
      USE CONSTANT,   only : tf     
      USE RADNCB,     only : COSZ1,cfrac,srdn
 ! tr_wd_TYPE,nPART are defined if drydep on
      USE TRACER_COM, only : ntm, tr_wd_TYPE, nPART, trname, tr_mm,
     &     dodrydep, F0_glob=>F0, HSTAR_glob=>HSTAR
#ifdef TRACERS_SPECIAL_Shindell
     &     , n_NOx
      USE TRCHEM_Shindell_COM, only : pNOx
#endif
      USE tracers_DRYDEP, only: NPOLY,IJREG,IJLAND,XYLAI,
     & DRYCOEFF,IJUSE,NTYPE,IDEP,IRI,IRLU,IRAC,IRGSS,IRGSO,
     & IRCLS,IRCLO,IVSMAX
     
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments
C
!@param XMWH2O molecular weight of water in KG/mole
      REAL*8, PARAMETER ::  XMWH2O = 18.d-3
!@var N the current tracer number
!@var LDT dummy loop variable for land type
!@var k,IW dummy loop variables
!@var RI ______ resistance
!@var RLU cuticular resistance for the bulk canopy
!@var RAC ______ resistance
!@var RGSS ______ resistance
!@var RGSO ______ resistance
!@var RCLS ______ resistance
!@var RCLO ______ resistance
!@var RSURFACE Bulk surface resistance for species K landtype LDT
!@var RI internal resistance (minimum stomatal resistance for
!@+   water vapor, per unit area of leaf)
!@var RAD0 downward solar radiation flux at surface (w/m2)
!@var RT correction term for bulk surface resistance for gases?
!@var RIX ______ resistance
!@var GFACT,GFACI,RGSX,RCLX,CZH ?
!@var RIXX corrected RIX
!@var RLUXX corrected RLU
!@var RDC aerodynamic resistance to elements in lower part of
!@+   the canopy or structure
!@var DTMP1,DTMP2,DTMP3,DTMP4 temp arrays for recipricol resistances
!@var VDS temp deposition velocity
!@var DUMMY1,DUMMY2,DUMMY3,DUMMY4 dummy temp variables
!@var TEMPK Surface air temperature (Kelvin)
!@var TEMPC Surface air temperature (oC)  
!@var byTEMPC 1/TEMPC
!@var I,J GCM grid box horizontal position
!@var ITYPE GCM surface type 1=ocean; 2=ocean ice; 3=land ice; 4=land
!@var TOTA total landuse area, excluding water and ice dep. types.
!@var problem_point logical that is true if ITYPE=4, but there are no
!@+   non-water, non-ice portions, according to Harvard dataset.
!@var OBK Monin-Obukhov length (m) (now is the PBL lmonin variable)
!@var ZHH boundary layer depth (m) ("mixing depth") (PBL dbl variable)
!@var USTARR Friction velocity (m s-1) (PBL ustar variable)
!@var tr_mm_temp temporary variable to hold tr_mm(k), etc.
!@var SUNCOS Cosine of solar zenith angle
!@var IOLSON integer index for olson surface types?
!@var VD deposition velocity temp array   (s m-1)
      REAL*8,  DIMENSION(ntm,NTYPE) :: RSURFACE
      REAL*8,  DIMENSION(ntm)   :: TOTA,VD,HSTAR,F0
      REAL*8,  DIMENSION(NTYPE) :: RI,RLU,RAC,RGSS,RGSO,RCLS,RCLO
      REAL*8 RT,RAD0,RIX,GFACT,GFACI,RDC,RIXX,RLUXX,RGSX,RCLX,VDS,
     &       DTMP1,DTMP2,DTMP3,DTMP4,CZH,DUMMY1,DUMMY2,DUMMY3,DUMMY4,
     &       TEMPC,byTEMPC,BIOFIT,DIFFG,tr_mm_temp,SUNCOS
      REAL*8, INTENT(IN) :: OBK,ZHH,USTARR,TEMPK
!@var dep_vel the deposition velocity = 1/bulk sfc. res. (m/s)
      REAL*8, INTENT(OUT), DIMENSION(NTM) :: dep_vel
      INTEGER k,n,LDT,II,IW,IOLSON
      INTEGER, INTENT(IN) :: I,J,ITYPE  
      LOGICAL problem_point
C
C Use cosine of the solar zenith angle from the radiation code,
C ...which seems to have a minumum of 0, like suncos used to have
C when defined in SCALERAD subroutine from Harvard CTM.
C
      SUNCOS = COSZ1(I,J)
c
C* Initialize VD and RSURFACE and reciprocal: 
      DO K = 1,ntm
       if(dodrydep(K))then
        DO LDT = 1,NTYPE 
          RSURFACE(K,LDT) = 0.
        END DO             
        VD(K)             = 0.
        dep_vel(K)        = 0.
       end if
      END DO    
C
C** TEMPK and TEMPC are surface air temperatures in K and in C  
ccc      TEMPK was BLDATA(I,J,2) now input as argument
      TEMPC = TEMPK-tf
      byTEMPC = 1.D0/TEMPC    
      RAD0 = srdn(I,J)
C               
C* Compute bulk surface resistance for gases.
C*   
C* Adjust external surface resistances for temperature;
C* from Wesely [1989], expression given in text on p. 1296.
C* Note: the sign of 4.0 was fixed after consulting w/ Harvard.
C* 
      RT = 1000.0*EXP(-TEMPC-4.0)
C
C  Check for species-dependant corrections to drydep parameters:      
      DO K = 1,ntm
       if(dodrydep(K)) then
         HSTAR(K)=HSTAR_glob(K)
         F0(K)=F0_glob(K)
#ifdef TRACERS_SPECIAL_Shindell
c       For NOx, sum deposition for NO and NO2
c       HSTAR: NO2 = 0.01, NO = 0.002, F0: NO2 = 0.1, NO = 0.0
         if(trname(K).eq.'NOx')then
           HSTAR(K)=pNOx(i,j,1)*0.01d0+(1.-pNOx(i,j,1))*.002d0
           F0(K)=pNOx(i,j,1)*0.1d0
         endif 
#endif
       end if
      END DO
C
      SELECT CASE(ITYPE)
      CASE(4)     ! LAND*************************************
C   
C    Get surface resistances - loop over land types LDT   
C**********************************************************************
C* The land types within each grid square are defined using the Olson
C* land-type database.  Each of the Olson land types is assigned a
C* corresponding "deposition land type" with characteristic values of
C* surface resistance components.  There are 74 Olson land-types but
C* only 11 deposition land-types (i.e., many of the Olson land types
C* share the same deposition characteristics).  Surface resistance 
C* components for the "deposition land types" are from Wesely [1989]
C* except for tropical forests [Jacob and Wofsy, 1990] 
C* and for tundra [Jacob et al., 1992].  All surface resistance
C* components are normalized to a leaf area index of unity.
C*   
C* Olson land types, deposition land types, and surface resistance
C* components are read from file DRYTBL='drydep.table'; check that file
C* for further details.   
C**********************************************************************
C* first check for any points that are ITYPE=4 but all ice/water:
        problem_point=.true.
        DO LDT=1, IJREG(I,J)
         II=IDEP(IJLAND(I,J,LDT)+1)
         if(II.ne.1.and.II.ne.11) problem_point=.false.
        END DO

        DO LDT = 1,IJREG(I,J)
          IF (IJUSE(I,J,LDT) .NE. 0) THEN
            IOLSON = IJLAND(I,J,LDT)+1  
            II = IDEP(IOLSON)
            ! exclude ice and water (except for problem points) :
            IF((II.ne.1.and.II.ne.11).or.problem_point) THEN
C**   
C** Here, we should probably put some provisions that if the GCM land
C** grid box is mostly snow-covered, set II=1 <<<<<<<<<<<<<<<<<<<<<<<
C**       
C* Read the internal resistance RI (minimum stomatal resistance for
C* water vapor, per unit area of leaf) from the IRI array; a '9999'
C* value means no deposition to stomata so we impose a very large   
C* value for RI.  
C*
            RI(LDT) = REAL(IRI(II))   
            IF (RI(LDT).GE. 9999.) RI(LDT)= 1.E12
C**
C** Cuticular resistances IRLU read in from 'drydep.table' 
C** are per unit area of leaf; divide them by the LAI to get
C** a cuticular resistance for the bulk canopy.  If IRLU is '9999' it
C** means there are no cuticular surfaces on which to deposit so we  
C** impose a very large value for RLU.  
C**  
            IF (IRLU(II).GE. 9999 .OR. XYLAI(I,J,LDT).LE.0.) THEN
              RLU(LDT)  = 1.E6
            ELSE
              RLU(LDT)= REAL(IRLU(II))/XYLAI(I,J,LDT) + RT
            ENDIF
C**  
C** The following are the remaining resistances for the Wesely   
C** resistance-in-series model for a surface canopy   
C** (see Atmos. Environ. paper, Fig.1).
C**
            RAC(LDT)  = MAX(REAL(IRAC(II)), 1.) 
            IF (RAC(LDT)  .GE. 9999.) RAC(LDT)  = 1.E12 
            RGSS(LDT) = MAX(REAL(IRGSS(II)) + RT ,1.D0)
            IF (RGSS(LDT) .GE. 9999.) RGSS(LDT) = 1.E12
            RGSO(LDT) = MAX(REAL(IRGSO(II)) + RT ,1.D0)
            IF (RGSO(LDT) .GE. 9999.) RGSO(LDT) = 1.E12
            RCLS(LDT) = REAL(IRCLS(II)) + RT
            IF (RCLS(LDT) .GE. 9999.) RCLS(LDT) = 1.E12
            RCLO(LDT) = REAL(IRCLO(II)) + RT
            IF (RCLO(LDT) .GE. 9999.) RCLO(LDT) = 1.E12
C                  
C**********************************************************************
C*  Adjust stomatal resistances for insolation and temperature:
C*  Temperature adjustment is from Wesely [1989], equation (3).
C*
C*  Light adjustment by the function BIOFIT described by Wang [1996].
C*  It combines
C*  - Local dependence of stomal resistance on the intensity I of light
C*    impinging the leaf; this is expressed as a mutliplicative
C*    factor I/(I+b) to the stomatal resistance where b = 50 W m-2
C*    (equation (7) of Baldocchi et al. [1987])
C*  - radiative transfer of direct and diffuse radiation in the
C*    canopy using equations (12)-(16) from Guenther et al. [1995]
C*  - separate accounting of sunlit and shaded leaves using
C*    equation (12) of Guenther et al. [1995]
C*  - partitioning of the radiation @ the top of the canopy into direct
C*    and diffuse components using a parameterization to results from
C*    an atmospheric radiative transfer model [Wang, 1996]
C*  The dependent variables of the function BIOFIT are the leaf area
C* index (XYLAI), the cosine of zenith angl (SUNCOS) and the fractional
C*  cloud cover (CFRAC).  The factor GFACI integrates the light
C*  dependence over the canopy depth; sp even though RI is input per 
C*  unit area of leaf it need not be scaled by LAI to yield a bulk
C*  canopy value because that's already done in the GFACI formulation.
C**********************************************************************

            RIX = RI(LDT)  
            IF (RIX .LT. 9999.) THEN
              IF (TEMPC.GT.0. .AND. TEMPC.LT.40.) THEN
                GFACT = 400.*byTEMPC/(40.0-TEMPC)  
              ELSE
                GFACT = 100.
              END IF
              IF (RAD0.GT.0. .AND. XYLAI(I,J,LDT).GT.0.) THEN  
                GFACI=1./BIOFIT
     *          (DRYCOEFF,XYLAI(I,J,LDT),SUNCOS,CFRAC(I,J))
              ELSE
                GFACI = 100.
              ENDIF   
              RIX = RIX*GFACT*GFACI   
            END IF
C*
C*    Compute aerodynamic resistance to lower elements in lower part
C*    of the canopy or structure, assuming level terrain -
C*    equation (5) of Wesely [1989].
C*  
            RDC = 100.*(1.+1000./(RAD0 + 10.))
C*
C*    Loop over species; species-dependent corrections to resistances
C*    are from equations (6)-(9) of Wesely [1989].
C*
            DO K = 1,ntm
             if(dodrydep(K))then
              tr_mm_temp = tr_mm(k)*1.D-3
#ifdef TRACERS_SPECIAL_Shindell
c             For NOx, use NO2 molecular weight to get collision diameter:
              if(trname(K).eq.'NOx') tr_mm_temp = 4.4d-2
#endif
C**           For non-aerosols:
              IF(tr_wd_TYPE(K).ne.nPART) THEN
                RIXX=RIX*DIFFG(TEMPK,XMWH2O)/DIFFG(TEMPK,tr_mm_temp)
     &          + 1./(HSTAR(K)/3.D3+100.*F0(K))
                IF(RLU(LDT).LT.9999.) THEN
                  RLUXX = RLU(LDT)/(1.D-5*HSTAR(K)+F0(K))
                ELSE
                  RLUXX = 1.E12
                END IF
C*
C* To prevent virtually zero resistance to species with huge HSTAR,
C* such as HNO3, a minimum value of RLUXX needs to be set. The
C* rationality of the existence of such a minimum is demonstrated by
C* the observed relationship between Vd(NOy-NOx) and Ustar in Munger
C* et al.[1996]; Vd(HNO3) never exceeds 2 cm s-1 in observations. The
C* corresponding minimum resistance is 50 s m-1.  This correction
C* was introduced by J.Y. Liang on 7/9/95.  
C* 
                IF(RLUXX .LT. 50.) RLUXX = 50.
                RGSX = 1./(1.D-5*HSTAR(K)/RGSS(LDT)+F0(K)/RGSO(LDT))   
                RCLX = 1./(1.D-5*HSTAR(K)/RCLS(LDT)+F0(K)/RCLO(LDT))      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Note from Greg & Gavin: it might be necessary to 
C limit some of these other resistances too:
                IF(RGSX.LT. 50.) RGSX= 50.
C               IF(RCLX.LT. 50.) RCLX= 50.
C               IF(RIXX.LT. 50.) RIXX= 50.
C               IF(RDC.LT. 50.) RDC= 50.
C               IF(RAC(LDT).LT. 50.) RAC(LDT)= 50.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C** Get bulk surface resistance of the canopy, from the network
C** of resistances in parallel and in series (Fig. 1 of Wesely [1989])
                DTMP1=1./RIXX
                DTMP2=1./RLUXX   
                DTMP3=1./(RAC(LDT)+RGSX)
                DTMP4=1./(RDC+RCLX)
                RSURFACE(K,LDT) = 1./(DTMP1+DTMP2+DTMP3+DTMP4)
              END IF                                 ! gases (above)
            
              IF (tr_wd_TYPE(K).eq.nPART) THEN ! aerosols (below)
C**             get surface deposition velocity for aerosols if needed;
C**             equations (15)-(17) of Walcek et al. [1986]            
                VDS = 0.002*USTARR
                IF(OBK.LT.0.)VDS = VDS*(1.+(-300./OBK)**0.6667)
                IF(OBK.EQ.0.) call stop_model('OBK=0 in TRDRYDEP',255)
                CZH  = ZHH/OBK
                IF(CZH.LT.-30.)VDS=0.0009*USTARR*(-CZH)**0.6667
C*
C* Set VDS to be less than VDSMAX (entry in input file divided by 1.D4)
C* VDSMAX is taken from Table 2 of Walcek et al. [1986].
C* Invert to get corresponding R
C*
                RSURFACE(K,LDT)=1.D0/MIN(VDS,1.D-4*REAL(IVSMAX(II)))
              END IF ! aerosols
C*
C*    Set max and min values for bulk surface resistances:
C*
              RSURFACE(K,LDT)=MAX(1.D0, MIN(RSURFACE(K,LDT), 9999.D0))
             end if! dodrydep
            END DO ! K loop   
            END IF ! end ice and water exclusion          
          END IF   ! IJUSE ne 0   
        END DO   ! LDT loop
C*
C* Loop through the different landuse types present in the grid square.
C* IJUSE is the fraction of the grid square occupied by surface LDT
C* in units of per mil (IJUSE=500 -> 50% of the grid square).  Add the
C* contribution of surface type LDT to the bulk surface resistance:
C*        
        TOTA = 0. ! total non-water, non-ice area
        DO LDT=1, IJREG(I,J)
          IF (IJUSE(I,J,LDT) .NE. 0) THEN
            IOLSON = IJLAND(I,J,LDT)+1  
            II = IDEP(IOLSON)
            ! exclude ice and water (except for problem points) :
            IF((II.ne.1.and.II.ne.11).or.problem_point) THEN 
            DO K = 1,ntm
             if(dodrydep(K))then
              VD(K) = VD(K) + 
     &        .001*REAL(IJUSE(I,J,LDT))/RSURFACE(K,LDT)
              TOTA(K)=TOTA(K)+.001*REAL(IJUSE(I,J,LDT))
             end if  ! dodrydep
            END DO   ! K
            END IF   ! end ice and water exclusion
          END IF     ! IJUSE ne 0
        END DO       ! LDT  
C* Calculate the deposition velocity, to be returned:
        DO K = 1,ntm
          if(dodrydep(K))dep_vel(K) = VD(K)/TOTA(K)
        END DO        
         
      CASE(1:3) ! OCEAN, OCEAN ICE, LANDICE *******************
        
C       Please see comments in land case above:     
        II = 1                  ! ice
        IF(ITYPE.eq.1) II  = 11 ! water
        LDT=1 ! just so we don't have to use all new variables
        RI(LDT)   = 1.E12 ! No stomatal deposition
        RLU(LDT)  = 1.E6  ! No cuticular surface deposition 
        RAC(LDT)  = 1.D0
        RGSS(LDT) = MAX(REAL(IRGSS(II)) + RT ,1.D0)
        IF (RGSS(LDT) .GE. 9999.) RGSS(LDT) = 1.E12
        RGSO(LDT) = MAX(REAL(IRGSO(II)) + RT ,1.D0)
        IF (RGSO(LDT) .GE. 9999.) RGSO(LDT) = 1.E12
        RCLS(LDT) = 1.E12
        RCLO(LDT) = REAL(IRCLO(II)) + RT
        IF (RCLO(LDT) .GE. 9999.) RCLO(LDT) = 1.E12
        RIX = RI(LDT)  
        IF (RIX .LT. 9999.) THEN
          IF (TEMPC.GT.0. .AND. TEMPC.LT.40.) THEN
            GFACT = 400.*byTEMPC/(40.0-TEMPC)
          ELSE
            GFACT = 100.
          END IF
          IF (RAD0.GT.0. .AND. XYLAI(I,J,LDT).GT.0.) THEN
            GFACI=1./BIOFIT
     *      (DRYCOEFF,XYLAI(I,J,LDT),SUNCOS,CFRAC(I,J))
          ELSE
            GFACI = 100.
          ENDIF
          RIX = RIX*GFACT*GFACI
        END IF
        RDC = 100.*(1.+1000./(RAD0 + 10.))
        DO K = 1,ntm
         if(dodrydep(K)) then
          tr_mm_temp = tr_mm(k)*1.D-3
#ifdef TRACERS_SPECIAL_Shindell
c         For NOx, use NO2 molecular weight to get collision diameter:
          if(trname(K).eq.'NOx') tr_mm_temp = 4.4d-2
#endif
          IF(tr_wd_TYPE(K).ne.nPART) THEN  ! NON-AEROSOLS
            RIXX = RIX*DIFFG(TEMPK,XMWH2O)/DIFFG(TEMPK,tr_mm_temp)
     &      + 1./(HSTAR(K)/3.D3+100.*F0(K))
            RLUXX = 1.E12
            RGSX = 1. / (1.D-5*HSTAR(K)/RGSS(LDT)+F0(K)/RGSO(LDT))
            RCLX = 1. / (1.D-5*HSTAR(K)/RCLS(LDT)+F0(K)/RCLO(LDT))
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Note from Greg & Gavin: it might be necessary to 
C limit some of these other resistances too:
            IF(RGSX.LT. 50.) RGSX= 50.
C           IF(RCLX.LT. 50.) RCLX= 50.
C           IF(RIXX.LT. 50.) RIXX= 50.
C           IF(RDC.LT. 50.) RDC= 50.
C           IF(RAC(LDT).LT. 50.) RAC(LDT)= 50.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            DTMP1=1./RIXX
            DTMP2=1./RLUXX   
            DTMP3=1./(RAC(LDT)+RGSX)
            DTMP4=1./(RDC+RCLX)
            RSURFACE(K,LDT) = 1./(DTMP1+DTMP2+DTMP3+DTMP4)
          ELSE IF (tr_wd_TYPE(K).eq.nPART) THEN ! AEROSOLS        
            VDS = 0.002*USTARR
            IF(OBK.LT.0.)VDS = VDS*(1.+(-300./OBK)**0.6667) 
            IF(OBK.EQ.0.) call stop_model('OBK=0 in TRDRYDEP',255)
            CZH  = ZHH/OBK
            IF(CZH.LT.-30.)VDS=0.0009*USTARR*(-CZH)**0.6667
            RSURFACE(K,LDT)=1.D0/MIN(VDS,1.D-4*REAL(IVSMAX(II)))
          END IF ! tracer type
          dep_vel(K)=1./MAX(1.D0, MIN(RSURFACE(K,LDT), 9999.D0))
         end if  ! do drydep
        END DO   ! tracer loop                  
C
      CASE DEFAULT
        call stop_model('ITYPE error in TRDRYDEP',255)
      END SELECT
C
      RETURN
      END SUBROUTINE get_dep_vel 
C      
C
      REAL*8 FUNCTION BIOFIT(COEFF1,XLAI1,SUNCOS1,CFRAC1)
!@sum BIOFIT Calculates the 'light correction' for the stomatal 
!@+   resistance?
!@auth Y.H. Wang
!@ver ? 

c
C**** GLOBAL parameters and variables:  
C
      USE tracers_DRYDEP, only : NPOLY
     
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments
C
!@param KK number of terms
!@param ND scaling factor for each variable
!@param X0 maximum for each variable
!@var COEFF1 drydep "coefficients" of fit?
!@var TERM,REALTERM ?
!@var XLAI1 local copy of leaf area index
!@var SUNCOS1 local copy of cosine of solar zenith angle
!@var CFRAC1 local copy of the cloud fraction
!@var K,K1,K2,K3,I loop index
!@var XLOW minimum for each variable
      INTEGER, PARAMETER :: KK=4
      INTEGER, PARAMETER :: ND(KK)=(/0,55,20,11/)
      REAL*8, PARAMETER  :: X0(KK)=(/0.,11.,1.,1./)
      REAL*8, DIMENSION(NPOLY), INTENT(IN) :: COEFF1
      REAL*8, DIMENSION(KK)                :: TERM
      REAL*8, DIMENSION(NPOLY)             :: REALTERM
      REAL*8, INTENT(IN) :: XLAI1,SUNCOS1,CFRAC1
      REAL*8 XLOW
      INTEGER I, K, K1, K2, K3

      TERM(1)=1. 
      TERM(2)=XLAI1 
      TERM(3)=SUNCOS1 
      TERM(4)=CFRAC1
C --- this section relaces SUNPARAM routine --- 
      DO I=2,KK ! NN
        TERM(I)=MIN(TERM(I),X0(I))
        IF (I.NE.KK) THEN ! hardcode of array position !
          XLOW=X0(I)/REAL(ND(I))
        ELSE
          XLOW= 0.
        END IF
        TERM(I)=MAX(TERM(I),XLOW)
        TERM(I)=TERM(I)/X0(I)
      END DO
C ---------------------------------------------
      K=0 
      DO K3=1,KK 
        DO K2=K3,KK 
          DO K1=K2,KK 
            K=K+1 
            REALTERM(K)=TERM(K1)*TERM(K2)*TERM(K3) 
          END DO 
        END DO 
      END DO 
      BIOFIT=0. 
      DO K=1,NPOLY 
        BIOFIT=BIOFIT+COEFF1(K)*REALTERM(K) 
      END DO 
      IF (BIOFIT.LT.0.1) BIOFIT=0.1 
      RETURN 
      END FUNCTION BIOFIT
C      
C
      REAL*8 FUNCTION DIFFG(TK,XM)
!@sum DIFFG to calculate tracer molecular diffusivity.
!@auth ? HARVARD CTM
!@ver 1.0 (based on CB436Tds3M23)
      USE CONSTANT, only : bygasc, gasc, pi, mair, avog
      IMPLICIT NONE
C=====================================================================
C  This function calculates the molecular diffusivity (m2 s-1) in air
C  for a gas X of molecular weight XM (kg) at temperature TK (K) and
C  pressure PRESS (Pa).
C  We specify the molecular weight of air (XMAIR) and the hard-sphere
C  molecular radii of air (RADAIR) & of the diffusing gas (RADX). The
C  molecular radius of air is given in a Table on p. 479 of Levine
C  [1988].  The Table also gives radii for some other molecules. Rather
C  than requesting the user to supply molecular radius we specify here
C  a generic value of 2.E-10 m for all molecules, which is good enough
C  in terms of calculating the diffusivity as long as molecule is not
C  too big.
C======================================================================
!@param XMAIR air molecular weight (KG/mole)
!@param AVOG Avogadro's number (molecules/mole)
!@param RADX hard-sphere molecular radius of the diffusing gas
!@param RADAIR hard-sphere molecular radius of air
!@param PRESS pressure (kg/s2/m) used to calculate molec. diffusivities
!@var TK passed local temperature in Kelvin
!@var XM molecular weight of diffusing gas (KG/mole)
!@var Z,DIAM ?
!@var FRPATH mean free path of the gas
!@var SPEEN average speed of the gas molecules
!@var AIRDEN local air density
      REAL*8, PARAMETER :: XMAIR = mair * 1.D-3,
     &                     RADX  = 1.5D-10,
     &                     RADAIR= 1.2D-10,
     &                     PRESS=  1.0d5
      REAL*8, INTENT(IN):: TK,XM
      REAL*8 Z,DIAM,FRPATH,SPEED,AIRDEN     
C
C* Calculate air density AIRDEN:
      AIRDEN = PRESS*avog*bygasc/TK ! can't we get this from the GCM?
C* Calculate the mean free path for gas X in air: eq. 8.5 of Seinfeld
C*  [1986]; DIAM is the collision diameter for gas X with air :
      Z = XM/XMAIR
      DIAM = RADX+RADAIR
      FRPATH = 1./(PI*SQRT(1.+Z)*AIRDEN*(DIAM**2.))
C* Calculate average speed of gas X; eq. 15.47 of Levine [1988]
      SPEED = SQRT(8.*gasc*TK/(PI*XM))
C* Calculate diffusion coefficient of gas X in air; eq. 8.9 of 
C* Seinfeld [1986]  :
      DIFFG = (3.*PI*0.03125)*(1.+Z)*FRPATH*SPEED
      RETURN
      END
C
C
      SUBROUTINE RDDRYCF
!@sum RDDRYCF Read the polynomial dry deposition coefficients from
!@+   the formatted input file 'drydep.coef'.
!@auth ? HARVARD CTM
!@ver ? 
!@calls openunit

c
C**** GLOBAL parameters and variables:  
C
      USE FILEMANAGER, only: openunit,closeunit
      USE tracers_DRYDEP, only: NPOLY,DRYCOEFF
     
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments
C
!@var DUM dummy variable for 80 character reading
!@var iu_data unit number for read
!@var I dummmy loop variable
      CHARACTER*80 DUM
      INTEGER iu_data,I
     
      call openunit('DRYCOEFF',iu_data,.false.,.true.)
      READ(iu_data,'(A80)') DUM
C--   read polynomial coefficients for drydep:   
      READ(iu_data,'(8(1PE10.2))') (DRYCOEFF(I),I=1,NPOLY)
      call closeunit(iu_data)
      RETURN
      END SUBROUTINE RDDRYCF
      

      SUBROUTINE RDLAND
!@sum RDLAND Read in land types and fractions(times 1000) from
!@+   the formatted input file:  'vegtype.global', also, Read Olson's
!@+   data from formatted input file 'drydep.table' (formerly done
!@+   in MODIN routine).
!@auth ? HARVARD CTM
!@ver ? 
!@calls openunit, MODIN
c
C      Results in non-gridded arrays:
C      IJREG(I,J), IJLAND(I,J,LDT), IJUSE(I,J,LDT)
C
C**** GLOBAL parameters and variables:  
C
      USE FILEMANAGER, only   : openunit,closeunit
      USE MODEL_COM, only    : im,jm
      USE tracers_DRYDEP, only: IJREG,IJLAND,IJUSE,IREG,NTYPE,IDEP,
     &                          IRI,IRLU,IRAC,IRGSS,IRGSO,IRCLS,
     &                          IRCLO,IVSMAX,NVEGTYPE
     
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments
C
!@var IREG # of landtypes in grid square (read in, gridded)
!@var I,J local lat lon index
!@var K,L,IDUMMY dummy loop indicies
!@var ILAND Land type ID for element LDT=1,IREG(I,J) (read in, gridded)
!@var IUSE per mil fraction of grid area with by land type element LDT
!@+   (read in, gridded)
!@var FRCLND land fraction (read in, gridded)
!@param NWAT number of olsons land types that are water
!@var NWAT2 Number of Olson's surface types that are water (read in)
!@var IWATER ID index for Olson's surface types that are water
!@var COM 70 characters of comments to read in?
!@var iu_data unit number for read
!@var iols temp read in index
!@var IZO roughness height for given surface type
      INTEGER, PARAMETER               :: NWAT=6
      CHARACTER*1 COM(70)
      REAL*8,  DIMENSION(IM,JM)        :: FRCLND
      INTEGER, DIMENSION(IM,JM,NTYPE)  :: ILAND, IUSE
      INTEGER, DIMENSION(NVEGTYPE)     :: IZO
      INTEGER I,J,K,L,NWAT2,IDUMMY,iu_data,iols
      INTEGER, DIMENSION(NWAT)         :: IWATER
      
      WRITE(6,*) 'READING land types and fractions ...'
      call openunit('VEGTYPE',iu_data,.false.,.true.)
      
100   READ(iu_data,'(20I4)',end=110) I,J,IREG(I,J),
     1                     (ILAND(I,J,K),K=1,IREG(I,J)),
     2                     (IUSE(I,J,K),K=1,IREG(I,J))
      GO TO 100
110   CONTINUE
      call closeunit(iu_data)
      DO J=1,JM
        DO I=1,IM
          FRCLND(I,J) = 1000.
          IJREG(I,J) = IREG(I,J)
          DO K=1,IJREG(I,J)
            IJLAND(I,J,K) = ILAND(I,J,K)
            IJUSE(I,J,K)  = IUSE(I,J,K)
            IF (IJLAND(I,J,K) .EQ. 0 )
     &      FRCLND(I,J) = FRCLND(I,J) - IJUSE(I,J,K)
          END DO ! K
          FRCLND(I,J) = FRCLND(I,J) * 1.E-3
        END DO   ! I
      END DO     ! J
C
C********* this section replaces call to MODIN *******************
      WRITE(6,*) 'READING drydep.table ...'
      call openunit('OLSON',iu_data,.false.,.true.)
      
      DO L = 1,5
        READ(iu_data,'(70A1)') COM
      END DO
C** Read Olson's surf. types, corresponding deposition surf. types, z0:
      DO L = 1,NVEGTYPE
         READ(iu_data,'(3I6)')  iols, IDEP(iols), IZO(iols)
      END DO
C** For the water surface types, zO is input as 1.E-4 m but is
C** recalculated elsewhere as function of wind speed.  Read the # of
C** Olson's surface types that are water (NWAT) and the corresponding
C** IDs (IWATER):
         READ(iu_data,'(70A1)') COM
         READ(iu_data,'(10I3)') NWAT2, (IWATER(I), I=1,NWAT)
         IF(NWAT2.ne.NWAT) 
     &   call stop_model('problem with NWAT vs NWAT2 in RDLAND',255)
C** Read parameters for each deposition surface type:
      DO IDUMMY=1,3
        READ(iu_data,'(70A1)') COM
      ENDDO
      DO L=1,NVEGTYPE
         READ(iu_data,'(9I5)', END=400) I,IRI(I),IRLU(I),IRAC(I),
     *   IRGSS(I),IRGSO(I),IRCLS(I),IRCLO(I),IVSMAX(I)
      END DO
 400  CONTINUE
      call closeunit(iu_data)
C******************** END MODIN SECTION **************************
      RETURN
      END SUBROUTINE RDLAND
c      
c                
      SUBROUTINE RDLAI
!@sum RDLAI Updates the Leaf Area Index (LAI) daily.
!@auth ? HARVARD CTM
!@ver ? 
!@calls READLAI
c
C**** GLOBAL parameters and variables:  
C
      USE tracers_DRYDEP, only: IJREG,XYLAI,XLAI,XLAI2,IREG
      USE MODEL_COM, only: IM,JM,JDmidOfM,JMperY,JDAY,JMON
C     
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments
C
!@var STARTDAY last (Julian) day of 1st half of month
!@var ISAVE 
!@var IMUL, ITD variables for getting right day or year
!@var M,K,I,J dummy loop variables
!@var JDAY current julian day
!@var byrITD reciprocol of real(ITD)
      INTEGER STARTDAY(JMperY),ISAVE,IMUL,ITD,M,K,I,J
      REAL*8 byrITD
      DATA ISAVE /0/
      SAVE ISAVE
c
      DO M=1,JMperY
        STARTDAY(M) = JDmidOfM(M) - 1
      END DO
c    
      IF (ISAVE.EQ.0) THEN
        ISAVE=1
        IF (JDAY.LT.STARTDAY(1)) THEN
          IMUL=365-STARTDAY(JMperY)+JDAY
          ITD = 31
        ELSE
          IMUL=JDAY-STARTDAY(JMON)
          ITD = STARTDAY(JMON+1) - STARTDAY(JMON)
        END IF
        byrITD = 1.E0/REAL(ITD)
        CALL READLAI
        DO J=1,JM
        DO I=1,IM
        DO K=1,IREG(I,J)
              XLAI2(I,J,K) = (XLAI2(I,J,K)-XLAI(I,J,K))*byrITD
              XLAI(I,J,K)=XLAI(I,J,K)+ XLAI2(I,J,K) * REAL(IMUL)
        END DO
        END DO
        END DO
      ELSE
        IF (JDAY.EQ.STARTDAY(JMON)) THEN
          ITD = STARTDAY(JMON+1) - STARTDAY(JMON)
          byrITD = 1.E0/REAL(ITD)
          CALL READLAI
          DO J=1,JM
          DO I=1,IM
          DO K=1,IREG(I,J)        
                XLAI2(I,J,K) = (XLAI2(I,J,K)-XLAI(I,J,K))*byrITD
          END DO
          END DO
          END DO
        ELSE   
          DO J=1,JM
          DO I=1,IM
          DO K=1,IREG(I,J)
                XLAI(I,J,K)=XLAI(I,J,K)+ XLAI2(I,J,K)
          END DO
          END DO
          END DO
        END IF
      END IF
      DO J=1,JM
        DO I=1,IM
          DO K=1,IJREG(I,J)
            XYLAI(I,J,K)=XLAI(I,J,K)
          END DO
        END DO
      END DO
      RETURN
      END SUBROUTINE RDLAI
c
c
      SUBROUTINE READLAI
!@sum READDLAI read in leaf area indicies from formatted file
!@+   (chem_files/lai##.global) for one month.
!@auth ? HARVARD CTM
!@ver ? 
!@calls openunit
c
C**** GLOBAL parameters and variables:  
C
      USE tracers_DRYDEP, only: XLAI,XLAI2,IREG
      USE MODEL_COM, only: IM,JM,JMON,JMperY
      USE FILEMANAGER, only: openunit,closeunit
C     
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments
C    
!@var CMONTH 2-digit of current month (character)  
!@var i,j,k dummy loop indicies
!@var IUNIT unuit number of current file being read
!@var MMM temp to hold JMON or 0
!@var INDEX ?
      CHARACTER*2, PARAMETER, DIMENSION(JMperY) :: CMONTH =
     &(/'01','02','03','04','05','06','07','08','09','10','11','12'/)
      INTEGER I,J,K,IUNIT,MMM,INDEX
C     
C lai's are in reference coordinates.  initialize them:
      DO J=1,JM
      DO I=1,IM
      DO K=1,IREG(I,J)
            XLAI(I,J,K)= 0.
            XLAI2(I,J,K)=0.
      END DO
      END DO
      END DO
      
C Read current month's lai:

      call openunit('LAI'//CMONTH(JMON),IUNIT,.false.,.true.)
10    READ(IUNIT,"(3I3,20F5.1)",END=20) I,J,INDEX,
     &(XLAI(I,J,K),K=1,INDEX)
      GOTO 10
 20   call closeunit(iunit)
      
C Read following months lai:
      IF(JMON.eq.12) THEN
        MMM=0
      ELSE
        MMM=JMON
      END IF
      call openunit('LAI'//CMONTH(MMM+1),IUNIT,.false.,.true.)
30    READ(IUNIT,"(3I3,20F5.1)",END=40) I,J,INDEX,
     &(XLAI2(I,J,K),K=1,INDEX)
      GOTO 30
40    call closeunit(iunit)
      RETURN
      END SUBROUTINE READLAI
#endif
