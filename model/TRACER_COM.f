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
!@var N_XXX: variable names of indices for tracers
!@var NTM_POWER: Power of 10 associated with each tracer (for printing)
!@var TR_MM: molecular mass of each tracer
!@var T_QLIMIT: if t_qlimit=.true. tracer is maintained as positive
      character*8, parameter :: trname(ntm)=
     * (/ 'Air     ','SF6     ','Rn222   ','CO2     '/)
      integer, parameter :: 
     * n_air=1,      n_SF6=2,    n_Rn222=3,n_CO2=4
      integer, parameter, dimension(ntm) :: ntm_power=
     * (/  -2,       -14,        -21,      -6/)
      real*8, parameter, dimension(ntm) :: tr_mm=
     * (/  mair,     146.01d0,   222.d0,   44.d0/)
      logical, dimension (ntm) :: t_qlimit=
     * (/  .true.,   .true.,    .true.,    .false./)
!@var ITIME_TR0: start time for each tracer
      integer, dimension(ntm) :: itime_tr0=0.
!@var TRM: Tracer array
      real*8, dimension(im,jm,lm,ntm) :: trm
!@var TRMOM: Second order moments for tracers
      real*8, dimension(nmom,im,jm,lm,ntm) :: trmom
!@var MTRACE: timing index for tracers
      integer mtrace
    
      END MODULE TRACER_COM

