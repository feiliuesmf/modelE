#include "rundeck_opts.h"

!@sum  TRACER_COM: Exists alone to minimize the number of dependencies
!@auth Jean Lerner/Gavin Schmidt
!@ver  1.0

      MODULE TRACER_COM
!@sum  TRACER_COM tracer variables
!@auth Jean Lerner
!@ver  1.0
      USE QUSDEF, only: nmom
      USE MODEL_COM, only: im,jm,lm
      USE CONSTANT, only: mair
      IMPLICIT NONE
      SAVE
!@param NTM number of tracers
      integer, parameter :: ntm=4

C**** Each tracer has a variable name and a unique index
!@var TRNAME: Name for each tracer >>> MUST BE LEFT-JUSTIFIED <<<
      character*8, parameter :: trname(ntm)=
     * (/ 'Air     ','SF6     ','Rn222   ','CO2     '/)
!@var N_XXX: variable names of indices for tracers
      integer, parameter :: 
     *     n_air=1,     n_SF6=2, n_Rn222=3,  n_CO2=4
!@var NTM_POWER: Power of 10 associated with each tracer (for printing)
      integer, parameter, dimension(ntm) :: ntm_power=
     *     (/   -2,         -14,       -21,       -6 /)
!@var TR_MM: molecular mass of each tracer
      real*8, parameter, dimension(ntm) :: tr_mm=
     *     (/ mair,    146.01d0,    222.d0,    44.d0 /)
!@var T_QLIMIT: if t_qlimit=.true. tracer is maintained as positive
      logical, parameter, dimension (ntm) :: t_qlimit=
     *     (/ .true.,   .true.,    .true.,  .false. /)
!@var needtrs: true if surface tracer value from PBL is required
      logical, parameter, dimension (ntm) :: needtrs=
     *     (/.false.,  .false.,   .false.,  .false. /)
!@var trdecy radioactive decay constant (1/s) (=0 for stable tracers)
      real*8, parameter, dimension (ntm) :: trdecy=
     *     (/    0d0,      0d0,    2.1d-6,      0d0 /)

C**** Note units for these parameters!
C**** Example: clay dust; trpdens=2.5d3, trradius=0.73d-6 
C****          silt dust; trpdens=2.65d3, trradius=6.1d-6 
C****
!@var trpdens tracer particle density (kg/m^3) (=0 for non-particle tracers)
      real*8, parameter, dimension (ntm) :: trpdens=
     *     (/    0d0,      0d0,    0d0,      0d0 /)
!@var trradius tracer effective radius (m) (=0 for non particle tracers)
      real*8, parameter, dimension (ntm) :: trradius=
     *     (/    0d0,      0d0,    0d0,      0d0 /)

!@var ITIME_TR0: start time for each tracer
      integer, dimension(ntm) :: itime_tr0=0.
!@var TRM: Tracer array
      real*8, dimension(im,jm,lm,ntm) :: trm
!@var TRMOM: Second order moments for tracers
      real*8, dimension(nmom,im,jm,lm,ntm) :: trmom
!@var MTRACE: timing index for tracers
      integer mtrace

#ifdef TRACERS_WATER
!@var TRWM tracer in cloud liquid water amount (kg)
      real*8, dimension(im,jm,lm,ntm) :: trwm
#endif

!@var ntsurfsrcmax maximum number of surface sources/sinks
      integer, parameter :: ntsurfsrcmax=10
!@var ntsurfsrc no. of non-interactive surface sources for each tracer
      integer, parameter, dimension(ntm) :: ntsurfsrc = (/
     *     0,       1,       1,       6 /)
    
      END MODULE TRACER_COM

      SUBROUTINE io_tracer(kunit,iaction,ioerr)
!@sum  io_tracer reads and writes tracer variables to file
!@auth Jean Lerner
!@ver  1.0
      USE MODEL_COM, only: ioread,iowrite,irsfic,irerun,lhead
      USE TRACER_COM
      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var HEADER Character string label for individual records
#ifdef TRACERS_WATER
      CHARACTER*80 :: HEADER, MODULE_HEADER = "TRACERW01"

      write (MODULE_HEADER(lhead+1:80),'(a,i2,a,a,i1,a,i2,a,i2,a)')
     *     'R8 TRM(im,jm,lm,',NTM,')',
     *     ',TRmom(',NMOM,',im,jm,lm,',NTM,'),trwm(im,jm,lm,',NTM,')' 
#else
      CHARACTER*80 :: HEADER, MODULE_HEADER = "TRACER01"

      write (MODULE_HEADER(lhead+1:80),'(a,i2,a,a,i1,a,i2,a)')
     *           'R8 TRM(im,jm,lm,',NTM,')',
     *  ',TRmom(',NMOM,',im,jm,lm,',NTM,')'

#endif

      SELECT CASE (IACTION)

      CASE (:IOWRITE) ! output to end-of-month restart file
        WRITE (kunit,err=10) MODULE_HEADER,TRM,TRmom
#ifdef TRACERS_WATER
     *     ,TRWM
#endif
      CASE (IOREAD:)          ! input from restart file
        SELECT CASE (IACTION)
        CASE (IRSFIC)   ! initial conditions
          READ (kunit)
        CASE (ioread,irerun) ! restarts
          READ (kunit,err=10) HEADER,TRM,TRmom
#ifdef TRACERS_WATER
     *       ,TRWM
#endif
          IF (HEADER(1:lhead).ne.MODULE_HEADER(1:lhead)) THEN
            PRINT*,"Discrepancy in module version",HEADER,MODULE_HEADER
            GO TO 10
          END IF
        END SELECT
      END SELECT

      RETURN
 10   IOERR=1
      RETURN
      END SUBROUTINE io_tracer
