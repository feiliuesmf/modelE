      MODULE RADNCB
!@sum  RADNCB Model radiation arrays and parameters
!@auth Original Development Team
!@ver  1.0
      USE MODEL_COM, only : im,jm,lm
      USE DOMAIN_DECOMP, only : grid
      USE RADPAR, only : S0
!@var S0 solar 'constant' needs to be saved between calls to radiation
      IMPLICIT NONE
      SAVE

C**** DEFAULT ORBITAL PARAMETERS FOR EARTH
C**** Note PMIP runs had specified values that do not necesarily
C**** coincide with those used as the default, or the output of ORBPAR.
C****                    OMEGT          OBLIQ        ECCEN
C**** DEFAULT (2000 AD): 282.9          23.44        0.0167
C**** PMIP CONTROL:      282.04         23.446       0.016724
C**** PMIP 6kyr BP:      180.87         24.105       0.018682
C**** PMIP LGM (21k):    294.42         22.949       0.018994
!@param OMEGT_def precession angle (degrees from vernal equinox)
      real*8, parameter :: omegt_def = 282.9d0
!@param OBLIQ_def obliquity angle  (degrees)
      real*8, parameter :: obliq_def = 23.44d0
!@param ECCN_def eccentricity
      real*8, parameter :: eccn_def  = .0167d0
!@var OMEGT,OBLIQ,ECCN actual orbital parameters used
      real*8 OMEGT,OBLIQ,ECCN

C**** Database parameters to control orbital parameter calculation
C**** Note: setting calc_orb_par with paleo_orb_yr=2000 does not produce
C**** exactly the same as the default values.
!@dbparam calc_orb_par = 1 to calc orbital parameters
      integer :: calc_orb_par = 0
!@dbparam paleo_orb_yr is paleo year (BP) for orbital calc
      real*8 :: paleo_orb_yr = -50.  ! (i.e. 2000AD)

!@var LM_REQ Extra number of radiative equilibrium layers
      INTEGER, PARAMETER :: LM_REQ=3
!@var dimrad_sv dimension sum of input fields saved for radia_only runs
      INTEGER, PARAMETER :: dimrad_sv=IM*JM*(7*LM+3*LM_REQ+23)
!@var RQT Radiative equilibrium temperatures above model top
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: RQT
!@var Tchg Total temperature change in adjusted forcing runs
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: Tchg
!@var SRHR(0) Solar   raditive net flux into the ground          (W/m^2)
!@var TRHR(0) Thermal raditive net flux into ground(W/O -StB*T^4)(W/m^2)
!@*   Note: -StB*T^4 is MISSING, since T may vary a lot betw. rad. calls
!@var SRHR(1->LM) Solar   raditive heating rate (W/m^2)  (short wave)
!@var TRHR(1->LM) Thermal raditive heating rate (W/m^2)  (long wave)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: SRHR,TRHR
!@var FSF Solar Forcing over each type (W/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: FSF
!@var FSRDIR Direct beam solar incident at surface (W/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: FSRDIR
!@var SRVISSURF Incident solar direct+diffuse visible at surface (W/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: SRVISSURF
!@var SRDN Total incident solar at surface (W/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: SRDN  ! saved in rsf

!@var CFRAC Total cloud fraction as seen be radiation
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: CFRAC ! saved in rsf
!@var RCLD Total cloud optical depth as seen be radiation
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: RCLD ! saved in rsf
!@var O3_rad_save 3D ozone saved from radiation for use elsewhere
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: O3_rad_save !saved in rsf
!@var O3_tracer_save 3D ozone saved elsewhere for use in radiation
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: O3_tracer_save!saved rsf
!@var KLIQ Flag indicating dry(0)/wet(1) atmosphere (memory feature)
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:,:) :: KLIQ ! saved in rsf
!@var COSZ1 Mean Solar Zenith angle for curr. physics(not rad) time step
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: COSZ1
!@dbparam S0X solar constant multiplication factor
      REAL*8 :: S0X = 1.
!@dbparam S0_yr,S0_day obs.date of solar constant (if 0: time var)
      INTEGER :: S0_yr = 1951 , S0_day = 182
!@dbparam CO2X,... scaling factors for CO2 N2O CH4 CFC11 CFC12 XGHG
      REAL*8 :: CO2X=1.,N2OX=1.,CH4X=1., CFC11X=1.,CFC12X=1.,XGHGX=1.
!@dbparam GHG_yr,GHG_day obs.date of well-mixed GHgases (if 0: time var)
      INTEGER :: GHG_yr = 1951 , GHG_day = 182
!@dbparam Volc_yr,Volc_day obs.date of Volc.Aerosols (if 0: time var)
      INTEGER :: Volc_yr = 1951 , Volc_day = 182
!@dbparam Aero_yr obs.year of troposph.Aerosols (if 0: time var)
      INTEGER :: Aero_yr = 1951    ! always use annual cycle
!@dbparam O3_yr obs.year of Ozone (if 0: time var)
      INTEGER :: O3_yr = 1951      ! always use annual cycle
!@dbparam crops_yr obs.year of crops (if 0: time var, -1: default)
      INTEGER :: crops_yr = -1
!@dbparam H2OstratX stratospheric water vapour multiplication factor
      REAL*8 :: H2OstratX = 1.
!@dbparam H2ObyCH4 if not 0: add CH4 produced H2O into layers 1->LM
      REAL*8 :: H2ObyCH4 = 1.
!@var dH2O  zonal H2O-prod.rate in kg/m^2/ppm_CH4/second in layer L
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: dH2O
!@var RSDIST,SIND,COSD orbit related variables computed once a day
      REAL*8 :: RSDIST,SIND,COSD
!@var ALB is SRNFLB(1)/(SRDFLB(1)+1.D-20),PLAVIS,PLANIR,ALBVIS,ALBNIR,
!@+       SRRVIS,SRRNIR,SRAVIS,SRANIR (see RADIATION)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:), TARGET :: ALB

!@var SALB (1.-broadband surface albedo) - saved in rsf
      REAL*8, POINTER, DIMENSION(:,:) :: SALB   ! = ALB(:,:,1)
!      EQUIVALENCE (SALB,ALB)
!@dbparam rad_interact_tr =1 for radiatively active tracers (default=0)
      INTEGER :: rad_interact_tr = 0

C**** Local variables initialised in init_RAD
!@var COE
      REAL*8, DIMENSION(LM+LM_REQ) :: COE
!@var PLE0,QL0 global parts of local arrays (to avoid OMP-copyin)
      REAL*8, DIMENSION(LM_REQ)       :: PLB0,SHL0
!@var SINJ,COSJ sines and cosines for zenith angle calculation
      REAL*8, ALLOCATABLE, DIMENSION(:) :: SINJ,COSJ

      END MODULE RADNCB

      SUBROUTINE ALLOC_RAD_COM(grid)
!@sum  To allocate arrays who sizes now need to be determined at
!@+    run-time
!@auth Rodger Abel
!@ver  1.0

      USE DOMAIN_DECOMP, ONLY : DYN_GRID
      USE DOMAIN_DECOMP, ONLY : GET
      USE MODEL_COM, ONLY : IM, JM, LM
      USE RADNCB, ONLY : LM_REQ
      USE RADNCB, ONLY : RQT, Tchg, SRHR, TRHR, FSF, FSRDIR, SRVISSURF,
     *     SRDN, CFRAC, RCLD, O3_rad_save, O3_tracer_save, KLIQ, COSZ1,
     *     dH2O, ALB, SALB, SINJ, COSJ

      IMPLICIT NONE
      TYPE (DYN_GRID), INTENT(IN) :: grid

      INTEGER :: J_0H, J_1H
      INTEGER :: IER

      CALL GET(grid, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)

      ALLOCATE( RQT(LM_REQ, IM, J_0H:J_1H),
     *     Tchg(LM+LM_REQ, IM, J_0H:J_1H),
     *     SRHR(0:LM, IM, J_0H:J_1H),
     *     TRHR(0:LM, IM, J_0H:J_1H),
     *     FSF(4,IM, J_0H:J_1H),
     *     FSRDIR(IM, J_0H:J_1H),
     *     SRVISSURF(IM, J_0H:J_1H),
     *     SRDN(IM, J_0H:J_1H),
     *     CFRAC(IM, J_0H:J_1H),
     *     RCLD(LM, IM, J_0H:J_1H),
     *     O3_rad_save(LM, IM, J_0H:J_1H),
     *     O3_tracer_save(LM, IM, J_0H:J_1H),
     *     KLIQ(LM,4, IM, J_0H:J_1H),
     *     COSZ1(IM, J_0H:J_1H),
     *     dH2O(J_0H:J_1H, LM, 12),
     *     ALB(IM, J_0H:J_1H, 9),
     *     SINJ(J_0H:J_1H),
     *     COSJ(J_0H:J_1H),
     *     STAT=IER)

      KLIQ = 1
      dH2O = 0.
      SALB => ALB(:,:,1)
      RETURN
      END SUBROUTINE ALLOC_RAD_COM

      SUBROUTINE io_rad(kunit,iaction,ioerr)
!@sum  io_rad reads and writes radiation arrays to file
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : ioread,iowrite,irsfic,irerun,ioread_single
     *         ,lhead,Kradia,irsficnt,irsficno
      USE RADNCB
      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var HEADER Character string label for individual records
      CHARACTER*80 :: HEADER, MODULE_HEADER = "RAD05"
!@var HEADER_F Character string label for records (forcing runs)
      CHARACTER*80 :: HEADER_F, MODULE_HEADER_F = "RADF"

      if (kradia.gt.0) then
        write (MODULE_HEADER_F(lhead+1:80),'(a8,i2,a15,i2,a7)')
     *    'R8 Tchg(',lm+LM_REQ,',ijm), I4 KLIQ(',lm,',4,ijm)'
        SELECT CASE (IACTION)
        CASE (:IOWRITE)            ! output to standard restart file
          WRITE (kunit,err=10) MODULE_HEADER_F,Tchg,kliq
        CASE (IOREAD:)
          SELECT CASE  (IACTION)
          CASE (ioread,irerun,ioread_single)  ! input for restart
            READ (kunit,err=10) HEADER_F,Tchg,kliq
          CASE (IRSFIC)  ! only way to start frc. runs
            READ (kunit,err=10) HEADER,RQT,kliq ; Tchg = 0.
          END SELECT
        END SELECT
      else

      MODULE_HEADER(lhead+1:80) = 'R8:Tq(3,ijm), I:LIQ(lm,4,ijm),'//
     *  ' S0,s+tHr,fs(4),fd,sv,al,sd,cf,o3'

      SELECT CASE (IACTION)
      CASE (:IOWRITE)            ! output to standard restart file
        WRITE (kunit,err=10) MODULE_HEADER,RQT,KLIQ
  ! rest needed only if MODRAD>0 at restart
     *    ,S0,SRHR,TRHR,FSF,FSRDIR,SRVISSURF,SALB,SRDN,CFRAC,RCLD
     *    ,O3_rad_save,O3_tracer_save
      CASE (IOREAD:)
        SELECT CASE  (IACTION)
        CASE (ioread,IRERUN)  ! input for restart, rerun or extension
          READ (kunit,err=10) HEADER,RQT,KLIQ
  ! rest needed only if MODRAD>0 at restart
     *       ,S0,SRHR,TRHR,FSF,FSRDIR,SRVISSURF,SALB,SRDN,CFRAC,RCLD
     *       ,O3_rad_save,O3_tracer_save
          IF (HEADER(1:LHEAD).NE.MODULE_HEADER(1:LHEAD)) THEN
            PRINT*,"Discrepancy in module version ",HEADER,MODULE_HEADER
            GO TO 10
          END IF
        CASE (IRSFIC,irsficnt,IRSFICNO)  ! restart file of prev. run
          READ (kunit,err=10) HEADER,RQT,KLIQ
        END SELECT
      END SELECT
      end if

      RETURN
 10   IOERR=1
      RETURN
      END SUBROUTINE io_rad

