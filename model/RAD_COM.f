      MODULE RADNCB
!@sum  RADNCB Model radiation arrays and parameters
!@auth Original Development Team
!@ver  1.0
      USE MODEL_COM, only : im,jm,lm
      USE RE001, only : S0
!@var S0 solar 'constant' needs to be saved between calls to radiation
      IMPLICIT NONE
      SAVE

!@var LM_REQ Extra number of radiative equilibrium layers
      INTEGER, PARAMETER :: LM_REQ=3
!@var dimrad_sv dimension sum of input fields saved for radia_only runs
      INTEGER, PARAMETER :: dimrad_sv=IM*JM*(7*LM+3*LM_REQ+20)
!@var RQT Radiative equilibrium temperatures above model top
      REAL*8, DIMENSION(LM_REQ,IM,JM) :: RQT
!@var Tchg Total temperature change in adjusted forcing runs
      REAL*8, DIMENSION(LM+LM_REQ,IM,JM) :: Tchg
!@var SRHR(0) Solar   raditive net flux into the ground          (W/m^2)
!@var TRHR(0) Thermal raditive net flux into ground(W/O -StB*T^4)(W/m^2)
!@*   Note: -StB*T^4 is MISSING, since T may vary a lot betw. rad. calls
!@var SRHR(1->LM) Solar   raditive heating rate (W/m^2)  (short wave)
!@var TRHR(1->LM) Thermal raditive heating rate (W/m^2)  (long wave)
      REAL*8, DIMENSION(0:LM,IM,JM) :: SRHR,TRHR
!@var FSF Solar Forcing over each type (W/m^2)
      REAL*8, DIMENSION(4,IM,JM) :: FSF
!@var COSZ1 Mean Solar Zenith angle for curr. physics(not rad) time step
      REAL*8, DIMENSION(IM,JM) :: COSZ1
!@dbparam S0X solar constant multiplication factor
      REAL*8 :: S0X = 1.
!@dbparam S0_yr,S0_day obs.date of solar constant (if 0: time var)
      INTEGER :: S0_yr = 1951 , S0_day = 182
!@dbparam CO2X carbon dioxide multiplication factor
      REAL*8 :: CO2X = 1.
!@dbparam GHG_yr,GHG_day obs.date of well-mixed GHgases (if 0: time var)
      INTEGER :: GHG_yr = 1951 , GHG_day = 182
!@dbparam CH4X methane  multiplication factor
      REAL*8 :: CH4X = 1.
!@dbparam Volc_yr,Volc_day obs.date of Volc.Aerosols (if 0: time var)
      INTEGER :: Volc_yr = 1951 , Volc_day = 182
!@dbparam Aero_yr obs.year of troposph.Aerosols (if 0: time var)
      INTEGER :: Aero_yr = 1951    ! always use annual cycle
!@dbparam O3_yr obs.year of Ozone (if 0: time var)
      INTEGER :: O3_yr = 1951      ! always use annual cycle
!@dbparam H2OstratX stratospheric water vapour multiplication factor
      REAL*8 :: H2OstratX = 1.
!@dbparam H2ObyCH4 if not 0: add CH4 produced H2O into layers 1->LM
      INTEGER :: H2ObyCH4 = 0
!@var dH2O  zonal H2O-prod.rate in kg/m^2/ppm_CH4/second in layer L
      REAL*8, DIMENSION(JM,LM) :: dH2O = 0.
!@var RSDIST,SIND,COSD orbit related variables computed once a day
      REAl*8 :: RSDIST,SIND,COSD

C**** Local variables initialised in init_RAD
!@var COE
      REAL*8, DIMENSION(LM+LM_REQ) :: COE
!@var PLE0,QL0 global parts of local arrays (to avoid OMP-copyin)
      REAL*8, DIMENSION(LM_REQ)       :: PLE0,QL0
!@var LLOW,LMID,LHI max levels for low, mid and high clouds
      INTEGER LLOW,LMID,LHI
!@var SINJ,COSJ sines and cosines for zenith angle calculation
      REAL*8, DIMENSION(JM) :: SINJ,COSJ

      END MODULE RADNCB

      SUBROUTINE io_rad(kunit,iaction,ioerr)
!@sum  io_rad reads and writes radiation arrays to file
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : ioread,iowrite,irsfic,irerun,ioread_single
     *         ,lhead,Kradia
      USE RADNCB
      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var HEADER Character string label for individual records
      CHARACTER*80 :: HEADER, MODULE_HEADER = "RAD01"

      if (kradia.gt.1) then
        write (MODULE_HEADER(lhead+1:80),'(a8,i2,a5)')
     *    'R8 Tchg(',lm+LM_REQ,',ijm)'
        SELECT CASE (IACTION)
        CASE (:IOWRITE)            ! output to standard restart file
          WRITE (kunit,err=10) MODULE_HEADER,Tchg
        CASE (IOREAD:)
          SELECT CASE  (IACTION)
          CASE (ioread,irerun,ioread_single)  ! input for restart
            READ (kunit,err=10) HEADER,Tchg
          CASE (IRSFIC)  ! only way to start adj.frc. run
            Tchg = 0.
          END SELECT
        END SELECT
      else

      MODULE_HEADER(lhead+1:80) = 'R8 Teq(3,im,jm),'//
     *  ' S0, s+tHr(0:lm,im,jm,2),fs(im,jm,4)'

      SELECT CASE (IACTION)
      CASE (:IOWRITE)            ! output to standard restart file
        WRITE (kunit,err=10) MODULE_HEADER,RQT
     *    ,S0,SRHR,TRHR,FSF      ! only needed if MODRAD > 0 at restart
      CASE (IOREAD:)
        SELECT CASE  (IACTION)
        CASE (ioread,IRERUN)  ! input for restart, rerun or extension
          READ (kunit,err=10) HEADER,RQT
     *       ,S0,SRHR,TRHR,FSF   ! only needed if MODRAD > 0 at restart
        CASE (IRSFIC)            ! start from restart file of prev. run
          READ (kunit,err=10) HEADER,RQT
        END SELECT
        IF (HEADER(1:LHEAD).NE.MODULE_HEADER(1:LHEAD)) THEN
          PRINT*,"Discrepancy in module version",HEADER,MODULE_HEADER
          GO TO 10
        END IF
      END SELECT
      end if

      RETURN
 10   IOERR=1
      RETURN
      END SUBROUTINE io_rad

