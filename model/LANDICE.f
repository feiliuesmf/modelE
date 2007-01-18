#include "rundeck_opts.h"

      MODULE LANDICE
!@sum  LANDICE contains variables and routines dealing with land ice
!@auth Original Development Team
!@ver  1.0
!@cont PRECLI,LNDICE
      USE CONSTANT, only : lhm,rhoi,shi
#ifdef TRACERS_WATER
      USE TRACER_COM, only : ntm
#endif
      IMPLICIT NONE
      SAVE
!@var Z1E,Z2LI first and second layer thicknesses for land ice
      REAL*8, PARAMETER :: Z1E=.1d0, Z2LI=2.9d0
!@var ACE1LI first layer land ice mass (kg/m^2)
      REAL*8, PARAMETER :: ACE1LI = Z1E*RHOI
!@var HC1LI heat capacity of first layer land ice (J/m^2)
      REAL*8, PARAMETER :: HC1LI = ACE1LI*SHI
!@var ACE2LI second layer land ice mass (kg/m^2)
      REAL*8, PARAMETER :: ACE2LI = Z2LI*RHOI
!@var HC2LI heat capacity of second layer land ice (J/m^2)
      REAL*8, PARAMETER :: HC2LI = ACE2LI*SHI
!@dbparam glmelt_on determines whether glacial melt is used for oceans
      INTEGER :: glmelt_on = 1   ! default is 'on'
!@dbparam glmelt_fac_nh is a factor to multiply glacial melt by in NH
      REAL*8  :: glmelt_fac_nh = 1.d0
!@dbparam glmelt_fac_sh is a factor to multiply glacial melt by in SH
      REAL*8  :: glmelt_fac_sh = 1.d0

!@var FWAREA_SH, FWAREA_NH area in each hemisphere for glmelt
      REAL*8  :: FWAREA_SH, FWAREA_NH

C**** The net fluxes of glacier mass into the ocean (per hemisphere)
!@var ACCPDA total accumulation per year for Antarctica (kg/yr)
!@var ACCPDG total accumulation per year for Greenland (kg/yr)
!@var EACCPDA total energy flux per year for Antarctica (kg/yr)
!@var EACCPDG total energy flux per year for Greenland (kg/yr)
      REAL*8 :: ACCPDA, ACCPDG,  EACCPDA, EACCPDG
#ifdef TRACERS_WATER /* TNL: inserted */
#ifdef TRACERS_OCEAN
!@var TRACCPDA total tracer flux per year for Antarctica (kg/yr)
!@var TRACCPDG total tracer flux per year for Greenland (kg/yr)
      REAL*8 :: TRACCPDA(NTM), TRACCPDG(NTM)
#endif
#endif   /* TNL: inserted */

      CONTAINS

      SUBROUTINE PRECLI(SNOW,TG1,TG2,PRCP,ENRGP,
#ifdef TRACERS_WATER
     *     TRSNOW,TRLI,TRPRCP,TRDIFS,TRUN0,
#endif
     *     EDIFS,DIFS,ERUN2,RUN0)
!@sum  PRECLI apply precipitation to land ice fraction
!@auth Original Development Team
!@ver  1.0
      IMPLICIT NONE
      REAL*8, INTENT(INOUT) :: SNOW,TG1,TG2
      REAL*8, INTENT(IN) :: PRCP,ENRGP
      REAL*8, INTENT(OUT) :: EDIFS,DIFS,ERUN2,RUN0
#ifdef TRACERS_WATER
!@var TRSNOW tracer amount in snow (kg/m^2)
      REAL*8, DIMENSION(NTM), INTENT(INOUT) :: TRSNOW
!@var TRLI tracer amount in land ice (kg/m^2)
      REAL*8, DIMENSION(NTM), INTENT(INOUT) :: TRLI
!@var TRPRCP tracer amount in precip (kg/m^2)
      REAL*8, DIMENSION(NTM), INTENT(IN) :: TRPRCP
!@var TRUN0 tracer runoff from ice (kg/m^2)
      REAL*8, DIMENSION(NTM), INTENT(OUT) :: TRUN0
!@var TRDIFS implicit tracer flux at base of ice (kg/m^2)
      REAL*8, DIMENSION(NTM), INTENT(OUT) :: TRDIFS
#endif

      REAL*8 HC1,DWATER
C**** initiallize output
      EDIFS=0. ; DIFS=0. ; ERUN2=0. ; RUN0=0.
#ifdef TRACERS_WATER
      TRDIFS(:)=0. ; TRUN0(:)=0.
#endif

      HC1=HC1LI+SNOW*SHI
      IF (ENRGP.GE.0.) THEN
        IF (ENRGP.GE.-TG1*HC1) THEN
C**** RAIN HEATS UP TG1 TO FREEZING POINT AND MELTS SOME SNOW OR ICE
          DWATER=(TG1*HC1+ENRGP)/LHM
          TG1=0.
          RUN0=DWATER+PRCP
          IF (DWATER.LE.SNOW) THEN
C**** RAIN MELTS SOME SNOW
#ifdef TRACERS_WATER
            TRUN0(:)=TRPRCP(:)+DWATER*TRSNOW(:)/SNOW
            TRSNOW(:)=TRSNOW(:)*(1.-DWATER/SNOW)
#endif
            SNOW=SNOW-DWATER
          ELSE
C**** RAIN MELTS ALL SNOW AND SOME ICE, ICE MOVES UP THROUGH THE LAYERS
            DIFS=SNOW-DWATER  ! < 0
            SNOW=0.
            TG1=-TG2*DIFS/ACE1LI
            EDIFS=DIFS*(TG2*SHI-LHM)
            ERUN2=EDIFS
#ifdef TRACERS_WATER
            TRUN0(:)=TRPRCP(:)+TRSNOW(:)-DIFS*TRLI(:)/(ACE2LI+ACE1LI)
            TRDIFS(:)=DIFS*TRLI(:)/(ACE2LI+ACE1LI)
            TRSNOW(:)=0.
#endif
          END IF
        ELSE
C**** RAIN COOLS TO FREEZING POINT AND HEATS UP TG1
          TG1=TG1+ENRGP/HC1
          RUN0=PRCP
#ifdef TRACERS_WATER
          TRUN0(:)=TRPRCP(:)
#endif
        END IF
      ELSE
C**** SNOW INCREASES SNOW AMOUNT AND SNOW TEMPERATURE RECOMPUTES TG1
        TG1=(TG1*HC1+ENRGP+LHM*PRCP)/(HC1+PRCP*SHI)
        SNOW=SNOW+PRCP
#ifdef TRACERS_WATER
        TRSNOW(:)=TRSNOW(:)+TRPRCP(:)
#endif
        IF (SNOW.gt.ACE1LI) THEN
C**** SNOW IS COMPACTED INTO ICE, ICE MOVES DOWN THROUGH THE LAYERS
          DIFS=SNOW-.9d0*ACE1LI
#ifdef TRACERS_WATER
          TRDIFS(:)=TRLI(:)*DIFS/(ACE2LI+ACE1LI)
          TRLI(:)=TRLI(:)+DIFS*(TRSNOW(:)/SNOW-TRLI(:)/(ACE2LI+ACE1LI))
          TRSNOW(:)=TRSNOW(:)*(1.-DIFS/SNOW)
#endif
          SNOW=.9d0*ACE1LI
          EDIFS=DIFS*(TG1*SHI-LHM)
          ERUN2=DIFS*(TG2*SHI-LHM)
          TG2=TG2+(TG1-TG2)*DIFS/ACE2LI
        END IF
      END IF
      RETURN
      END SUBROUTINE PRECLI

      SUBROUTINE LNDICE(SNOW,TG1,TG2,F0DT,F1DT,EVAP,
#ifdef TRACERS_WATER
     *     TRSNOW,TRLI,TREVAP,TRDIFS,TRUN0,
#endif
     *     EDIFS,DIFS,RUN0)
!@sum  LNDICE apply surface fluxes to land ice fraction
!@auth Original Development team
!@ver  1.0
      IMPLICIT NONE
      REAL*8, INTENT(INOUT) :: SNOW,TG1,TG2
      REAL*8, INTENT(IN) :: F0DT,F1DT,EVAP
      REAL*8, INTENT(OUT) :: EDIFS,DIFS,RUN0
#ifdef TRACERS_WATER
C**** Tracer concentration in ice assumed constant
!@var TRSNOW tracer amount in snow (kg/m^2)
      REAL*8, DIMENSION(NTM), INTENT(INOUT) :: TRSNOW
!@var TRLI tracer amount in land ice (kg/m^2)
      REAL*8, DIMENSION(NTM), INTENT(INOUT) :: TRLI
!@var TREVAP tracer amount in evaporation (kg/m^2)
      REAL*8, DIMENSION(NTM), INTENT(IN) :: TREVAP
!@var TRUN0 tracer runoff from ice (kg/m^2)
      REAL*8, DIMENSION(NTM), INTENT(OUT) :: TRUN0
!@var TRDIFS implicit tracer flux at base of ice (kg/m^2)
      REAL*8, DIMENSION(NTM), INTENT(OUT) :: TRDIFS
      INTEGER N
#endif
      REAL*8 :: SNANDI,HC1,ENRG1

C**** initiallize output
      EDIFS=0. ; DIFS=0. ; RUN0=0.
#ifdef TRACERS_WATER
      TRDIFS(:)=0. ; TRUN0(:)=0.
#endif

C**** CALCULATE TG1
      SNANDI=SNOW+ACE1LI-EVAP
      HC1=SNANDI*SHI
      ENRG1=F0DT+EVAP*(TG1*SHI-LHM)-F1DT
      IF (ENRG1.GT.-TG1*HC1) THEN
C**** FLUXES HEAT UP TG1 TO FREEZING POINT AND MELT SOME SNOW AND ICE
        RUN0=(ENRG1+TG1*HC1)/LHM
        TG1=0.
        SNANDI=SNANDI-RUN0
#ifdef TRACERS_WATER
        IF (SNOW-EVAP.gt.RUN0) THEN  ! only snow melts
          TRUN0(:)=RUN0*(TRSNOW(:)-TREVAP(:))/(SNOW-EVAP)
          TRSNOW(:)=TRSNOW(:)-TREVAP(:)-TRUN0(:)
        ELSE ! all snow + some ice melts
          TRUN0(:)=TRSNOW(:)-TREVAP(:)+(ACE1LI-SNANDI)*TRLI(:)/(ACE2LI
     *         +ACE1LI)
c              ! done below
c          TRSNOW(:)=0.
c          TRLI(:)=TRLI(:)*(1.- ((ACE1LI-SNANDI)/(ACE2LI+ACE1LI)))
          do n=1,ntm
            TRUN0(N)=MAX(TRUN0(N),0d0)
          end do
        END IF
#endif
      ELSE
C**** FLUXES RECOMPUTE TG1 WHICH IS BELOW FREEZING POINT
        TG1=TG1+ENRG1/HC1
#ifdef TRACERS_WATER
        if (SNOW.gt.EVAP) THEN
          TRSNOW(:)=TRSNOW(:)-TREVAP(:)
c        else   ! done below
c          TRLI(:)=TRLI(:)-(TREVAP(:)-TRSNOW(:))
c          TRSNOW(:)=0.
        end if
#endif
      END IF
      IF (SNANDI.LT.ACE1LI) THEN
C**** SOME ICE HAS MELTED OR EVAPORATED, TAKE IT FROM G2
        SNOW=0.
        DIFS=SNANDI-ACE1LI
        TG1=(TG1*SNANDI-TG2*DIFS)/ACE1LI
        EDIFS=DIFS*(TG2*SHI-LHM)
#ifdef TRACERS_WATER
        TRDIFS(:)=TRLI(:)*(DIFS/(ACE2LI+ACE1LI))
        TRLI(:)=TRLI(:)-DIFS*TRLI(:)/(ACE2LI+ACE1LI)+
     *       TRSNOW(:)-TREVAP(:)-TRUN0(:)
        TRSNOW(:)=0.
#endif
      ELSE
        SNOW=SNANDI-ACE1LI
      END IF
C**** CALCULATE TG2
      TG2=TG2+F1DT/HC2LI

      RETURN
      END SUBROUTINE LNDICE

      END MODULE LANDICE

      MODULE LANDICE_COM
!@sum  LANDICE_COM contains the model arrays for land ice
!@auth Gavin Schmidt
!@ver  1.0
!@cont io_landice
      USE MODEL_COM, only : im,jm
#ifdef TRACERS_WATER
      USE TRACER_COM, only : ntm
#endif
      IMPLICIT NONE
      SAVE
!@var SNOWLI snow amount on land ice (kg/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: SNOWLI 
!@var TLANDI temperature of each land ice layer (C)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: TLANDI 
!@var MDWNIMP downward implicit ice amount accumulator (kg)
!@var EDWNIMP downward implicit energy amount accumulator (J)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: MDWNIMP, EDWNIMP 
!@var LOC_GLA, LOC_GLG mask for glmelt around Antarctica and Greeland
      LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: LOC_GLA, LOC_GLG

#ifdef TRACERS_WATER
!@var TRSNOWLI tracer amount in land ice snow (kg/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: TRSNOWLI 
!@var TRLNDI tracer amount in land ice (kg/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: TRLNDI
!@var TRLI0 default tracer conc. for land ice (kg/kg)
      REAL*8, DIMENSION(NTM) :: TRLI0
!@var TDWNIMP downward implicit tracer amount accumulator (kg)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: TRDWNIMP
#endif

      END MODULE LANDICE_COM

      SUBROUTINE ALLOC_LANDICE_COM(grid)
!@sum  To allocate arrays whose sizes now need to be determined at
!@+    run time
!@auth NCCS (Goddard) Development Team
!@ver  1.0
      USE DOMAIN_DECOMP, ONLY : DIST_GRID
      USE RESOLUTION, ONLY : IM,LM
      USE LANDICE_COM, ONLY : SNOWLI, TLANDI, MDWNIMP, EDWNIMP, LOC_GLA,
     *     LOC_GLG
#ifdef TRACERS_WATER
      USE LANDICE_COM, ONLY : TRSNOWLI, TRLNDI, TRDWNIMP
      USE TRACER_COM, only : ntm
#endif
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(IN) :: grid

      INTEGER :: J_1H, J_0H
      INTEGER :: IER

      J_0H = grid%J_STRT_HALO
      J_1H = grid%J_STOP_HALO

      ALLOCATE( SNOWLI(IM,J_0H:J_1H),
     *          TLANDI(2,IM,J_0H:J_1H),
     *          MDWNIMP(IM,J_0H:J_1H),
     *          EDWNIMP(IM,J_0H:J_1H),
     *          LOC_GLA(IM,J_0H:J_1H),
     *          LOC_GLG(IM,J_0H:J_1H),
     *          STAT=IER)
#ifdef TRACERS_WATER
      ALLOCATE( TRSNOWLI(NTM,IM,J_0H:J_1H),
     *          TRLNDI  (NTM,IM,J_0H:J_1H),
     *          TRDWNIMP(NTM,IM,J_0H:J_1H),
     *          STAT=IER)
#endif 

      RETURN
      END SUBROUTINE ALLOC_LANDICE_COM

      SUBROUTINE io_landice(kunit,iaction,ioerr)
!@sum  io_landice reads and writes landice variables to file
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : ioread,iowrite,lhead,irsfic,irsficno,irerun
      USE DOMAIN_DECOMP, only : grid, GET, AM_I_ROOT
      USE DOMAIN_DECOMP, only : PACK_DATA, UNPACK_DATA, PACK_COLUMN
      USE DOMAIN_DECOMP, only : UNPACK_COLUMN, ESMF_BCAST,
     *     BACKSPACE_PARALLEL
      USE LANDICE_COM
      USE LANDICE
      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var HEADER Character string label for individual records
      CHARACTER*80 :: HEADER, MODULE_HEADER = "GLAIC02"
!@var SNOWLI_GLOB dummy global array used for esmf i/o
      REAL*8, DIMENSION(IM,JM) :: SNOWLI_GLOB
!@var TLANDI_GLOB dummy global array used for esmf i/o
      REAL*8, DIMENSION(2,IM,JM) :: TLANDI_GLOB
      REAL*8, DIMENSION(IM,JM) :: MDWNIMP_GLOB, EDWNIMP_GLOB
      INTEGER :: J_0H,J_1H, IN
#ifdef TRACERS_WATER
!@var TRHEADER Character string label for individual records
      CHARACTER*80 :: TRHEADER, TRMODULE_HEADER = "TRGLAC02"
!@var TRSNOWLI_GLOB  dummy global arrays used for i/o
      REAL*8, DIMENSION(NTM,IM,JM) :: TRSNOWLI_GLOB
!@var TRLNDI_GLOB  dummy global arrays used for i/o
      REAL*8, DIMENSION(NTM,IM,JM) :: TRLNDI_GLOB
      REAL*8, DIMENSION(NTM,IM,JM) :: TRDWNIMP_GLOB
      write (TRMODULE_HEADER(lhead+1:80)
     *     ,'(a7,i3,a)')'R8 dim(',NTM,',im,jm):TRSNOWLI,TRLNDI,TRDWN'
#ifdef TRACERS_OCEAN
     *     //',TRACC'
#endif
#endif
      
      MODULE_HEADER(lhead+1:80) = 'R8 SNOW(im,jm),T(2,im,jm),MDWN,EDWN'
     *     //',MACC,EACC'

      SELECT CASE (IACTION)

      CASE (:IOWRITE)            ! output to standard restart file
C**** Gather into global arrays
        CALL PACK_DATA(grid,SNOWLI,SNOWLI_GLOB)
        CALL PACK_COLUMN( grid,TLANDI,TLANDI_GLOB)
        CALL PACK_COLUMN( grid,MDWNIMP,MDWNIMP_GLOB)
        CALL PACK_COLUMN( grid,EDWNIMP,EDWNIMP_GLOB)
        IF (AM_I_ROOT())
     &       WRITE (kunit,err=10) MODULE_HEADER,SNOWLI_GLOB,TLANDI_GLOB
     *       ,MDWNIMP_GLOB,EDWNIMP_GLOB,ACCPDA,ACCPDG,EACCPDA,EACCPDG
#ifdef TRACERS_WATER
C**** Gather into global arrays
          CALL PACK_COLUMN( grid,TRSNOWLI,TRSNOWLI_GLOB )
          CALL PACK_COLUMN( grid,TRLNDI,  TRLNDI_GLOB  )
          CALL PACK_COLUMN( grid,TRDWNIMP,TRDWNIMP_GLOB )
        IF (AM_I_ROOT())
     &         WRITE (kunit,err=10) TRMODULE_HEADER,TRSNOWLI_GLOB
     *         ,TRLNDI_GLOB,TRDWNIMP_GLOB
#ifdef TRACERS_OCEAN
     *         ,TRACCPDA,TRACCPDG
#endif
#endif

      CASE (IOREAD:)            ! input from restart file
        if ( AM_I_ROOT() ) then
          READ (kunit,err=10) HEADER
          CALL BACKSPACE_PARALLEL(kunit)
          IF (HEADER(1:LHEAD).EQ.MODULE_HEADER(1:LHEAD)) THEN
            READ (kunit,err=10) HEADER,SNOWLI_GLOB,TLANDI_GLOB
     *           ,MDWNIMP_GLOB,EDWNIMP_GLOB,ACCPDA,ACCPDG,EACCPDA
     *           ,EACCPDG
          ELSEIF (HEADER(1:LHEAD).EQ."GLAIC01") THEN
            READ (kunit,err=10) HEADER,SNOWLI_GLOB,TLANDI_GLOB
          ELSE
            PRINT*,"Discrepancy in module version ",HEADER,MODULE_HEADER
            GO TO 10
          END IF
        end if
C****** Get useful ESMF parameters
        CALL GET( GRID, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H )
C****** Load data into distributed arrays
        CALL UNPACK_DATA( GRID, SNOWLI_GLOB, SNOWLI)
        CALL UNPACK_COLUMN( GRID, TLANDI_GLOB, TLANDI)
        CALL UNPACK_COLUMN( GRID, MDWNIMP_GLOB, MDWNIMP)
        CALL UNPACK_COLUMN( GRID, EDWNIMP_GLOB, EDWNIMP)
        call ESMF_BCAST(grid,  ACCPDA)
        call ESMF_BCAST(grid,  ACCPDG)
        call ESMF_BCAST(grid, EACCPDA)
        call ESMF_BCAST(grid, EACCPDG)

#ifdef TRACERS_WATER
        SELECT CASE (IACTION)
        CASE (IRERUN,IOREAD,IRSFIC,IRSFICNO)    ! from reruns/restarts
          if ( AM_I_ROOT() ) then
            READ (kunit,err=10) TRHEADER
            CALL BACKSPACE_PARALLEL(kunit)
            IF (TRHEADER(1:LHEAD).EQ.TRMODULE_HEADER(1:LHEAD)) THEN
              READ (kunit,err=10) TRHEADER,TRSNOWLI_GLOB,TRLNDI_GLOB
     *             ,TRDWNIMP_GLOB
#ifdef TRACERS_OCEAN
     *             ,TRACCPDA,TRACCPDG
#endif
            ELSEIF (TRHEADER(1:LHEAD).EQ."TRGLAC01") THEN
              READ (kunit,err=10) TRHEADER,TRSNOWLI_GLOB,TRLNDI_GLOB
            ELSE
              PRINT*,"Discrepancy in module version ",TRHEADER
     *             ,TRMODULE_HEADER
              GO TO 10
            END IF
          end if                !..am_i_root
C********* Load data into distributed arrays
          CALL UNPACK_COLUMN(GRID, TRSNOWLI_GLOB, TRSNOWLI)
          CALL UNPACK_COLUMN(GRID, TRLNDI_GLOB,   TRLNDI)
          CALL UNPACK_COLUMN(GRID, TRDWNIMP_GLOB, TRDWNIMP)
#ifdef TRACERS_OCEAN
          call ESMF_BCAST(grid, TRACCPDA)
          call ESMF_BCAST(grid, TRACCPDG)
#endif
        END SELECT
#endif

      END SELECT

      RETURN
 10   IOERR=1
      RETURN
      END SUBROUTINE io_landice


