



























module OldTracer_mod
  implicit none
  private

  public :: OldTracer_type
  public :: initializeOldTracers
  public :: makeNewTracers
  public :: addTracer

  public :: nGAS, nPART, nWATER

  public :: set_tr_mm, tr_mm 
public :: set_ntm_power, ntm_power 
public :: set_t_qlimit, t_qlimit 
public :: set_ntsurfsrc, ntsurfsrc 
public :: set_needtrs, needtrs 
public :: set_trdecay, trdecay 
public :: set_itime_tr0, itime_tr0 
public :: set_mass2vol, mass2vol 
public :: set_vol2mass, vol2mass 
public :: set_dodrydep, dodrydep 
public :: set_F0, F0 
public :: set_HSTAR, HSTAR 
public :: set_do_fire, do_fire 
public :: set_nBBsources, nBBsources 
public :: set_emisPerFireByVegType, emisPerFireByVegType 
public :: set_trpdens, trpdens 
public :: set_trradius, trradius 
public :: set_tr_wd_TYPE, tr_wd_TYPE 
public :: set_tr_RKD, tr_RKD 
public :: set_tr_DHD, tr_DHD 
public :: set_fq_aer, fq_aer 
public :: set_rc_washt, rc_washt 
public :: set_isDust, isDust 
public :: set_tr_H2ObyCH4, tr_H2ObyCH4 
public :: set_dowetdep, dowetdep 
public :: set_ntrocn, ntrocn 
public :: set_conc_from_fw, conc_from_fw 
public :: set_trglac, trglac 
public :: set_ntisurfsrc, ntisurfsrc 
public :: set_TRLI0, TRLI0 


  public :: set_trw0, trw0 

  public :: set_TRSI0, TRSI0 


  integer, parameter :: MAX_LEN_NAME = 8
!**** parameters for tr_wd_TYPE
!@param nGAS   index for wetdep tracer type = gas
!@param nPART  index for wetdep tracer type = particle/aerosol
!@param nWATER index for wetdep tracer type = water
      integer, parameter :: nGAS=1, nPART=2, nWATER=3

  type OldTracer_type
    private
    character(len=MAX_LEN_NAME) :: name

!@var tr_mm: !@var TR_MM: molecular mass of each tracer (g/mole)
    real*8 :: tr_mm  
!@var ntm_power: Power of 10 associated with each tracer (for printing)
    integer :: ntm_power  
!@var t_qlimit: if t_qlimit=.true. tracer is maintained as positive
    logical :: t_qlimit = .true. 
!@var ntsurfsrc: no. of non-interactive surface sources for each tracer
    integer :: ntsurfsrc = 0 
!@var needtrs: true if surface tracer value from PBL is required
    logical :: needtrs = .false. 
!@var trdecay: radioactive decay constant (1/s) (=0 for stable tracers)
    real*8 :: trdecay = 0.d0 
!@var itime_tr0: start time for each tracer (hours)
    integer :: itime_tr0  
!@var mass2vol: mass to volume ratio = mair/tr_mm
    real*8 :: mass2vol  
!@var vol2mass: volume to mass ratio = tr_mm/mair
    real*8 :: vol2mass  
!@var dodrydep: true if tracer should undergo dry deposition
    logical :: dodrydep  
!@var F0: reactivity factor for oxidation of biological substances
    real*8 :: F0 = 0.d0 
!@var HSTAR: Henry's Law const for tracer dry deposition. mole/(L atm)
!@+   Same as the tr_RKD wet dep variable
    real*8 :: HSTAR = 0.d0 
!@var do_fire: true if tracer should have emissions via flammability
    logical :: do_fire = .false. 
!@var nBBsources: number of sources attributed to biomass burning
    integer :: nBBsources = 0 
!@var emisPerFireByVegType: emisPerFireByVegType tracer emissions per fire count as a
!@+ function of 12 standard GISS (VDATA) vegetation types
    real*8, dimension(12) :: emisPerFireByVegType = 0.0d0 
!@var trpdens: tracer particle density (kg/m^3)
!@+               (=0 for non-particle tracers)
    real*8 :: trpdens = 0.d0 
!@var trradius: tracer effective radius (m) (=0 for non particle tracers)
    real*8 :: trradius = 0.d0 
!@var tr_wd_TYPE: tracer wet dep type (gas, particle, water)
    integer :: tr_wd_TYPE = nGas 
!@var tr_RKD: Henry's Law coefficient (in mole/Joule please !)
    real*8 :: tr_RKD = 0.d0 
!@var tr_DHD: coefficient of temperature-dependence term of Henry's
!@+   Law coefficient (in Joule/mole please !)
    real*8 :: tr_DHD = 0.d0 
!@var fq_aer: fraction of aerosol that condenses
    real*8 :: fq_aer = 0.d0 
!@var rc_washt: aerosol washout rate
    real*8 :: rc_washt = 1.d-1 
!@var isDust: index array for testing if tracer is a dust type
    integer :: isDust = 0 
!@var tr_H2ObyCH4: conc. of tracer in water from methane oxidation
    real*8 :: tr_H2ObyCH4 = 0.d0 
!@var dowetdep: true if tracer has some form of wet deposition
    logical :: dowetdep = .false. 
!@var ntrocn: scaling power factor for ocean/ice tracer concentrations
    integer :: ntrocn = 0 
!@var conc_from_fw: true if ocean conc is defined using fresh water
    logical :: conc_from_fw = .true. 
!@var trglac: tracer ratio in glacial runoff to ocean (kg/kg)
    real*8 :: trglac  
!@var ntisurfsrc: no. of interactive surface sources for each tracer
    integer :: ntisurfsrc  
!@var TRLI0: default tracer conc. for land ice (kg/kg)
    real*8 :: TRLI0 = 0.d0 


!@var TRSI0 default tracer conc. in sea ice (kg/m^2)
      REAL*8 :: TRSI0 = 0.0
!@var TRW0 default tracer concentration in water (kg/kg)
      real*8 :: trw0 = 0.
  end type OldTracer_type

  interface trw0
    module procedure trw0_all
    module procedure trw0_one
  end interface

  interface trsi0
    module procedure trsi0_all
    module procedure trsi0_one
  end interface

  integer, save :: numTracers = 0
  type (OldTracer_type), allocatable :: tracers(:)

contains

  subroutine initializeOldTracers()
    allocate(tracers(0))
  end subroutine initializeOldTracers

  subroutine addTracer(name)
    character(len=*), intent(in) :: name
    type (OldTracer_type), allocatable :: tmp(:)

    allocate(tmp(numTracers))
    tmp = tracers
    deallocate(tracers)

    allocate(tracers(numTracers+1))
    tracers(1:numTracers) = tmp
    deallocate(tmp)

    numTracers = numTracers + 1
    tracers(numTracers)%name = trim(name)

  end subroutine addTracer

  function makeNewTracers() result (bundle)
    use Tracers_mod, only: Tracer_type
    use Tracers_mod, only: Tracer
    use TracerBundle_mod, only: TracerBundle_type
    use TracerBundle_mod, only: TracerBundle
    use TracerBundle_mod, only: addTracer
    use Tracers_mod, only: setProperty
    use Tracers_mod, only: clean
    type (TracerBundle_type), pointer :: bundle

    type (Tracer_type) :: aTracer
    integer :: i

    bundle = TracerBundle()
    do i = 1, numTracers
      aTracer = Tracer(tracers(i)%name)
      call setProperty(aTracer, 'tr_mm', tr_mm(i))
call setProperty(aTracer, 'ntm_power', ntm_power(i))
call setProperty(aTracer, 't_qlimit', t_qlimit(i))
call setProperty(aTracer, 'ntsurfsrc', ntsurfsrc(i))
call setProperty(aTracer, 'needtrs', needtrs(i))
call setProperty(aTracer, 'trdecay', trdecay(i))
call setProperty(aTracer, 'itime_tr0', itime_tr0(i))
call setProperty(aTracer, 'mass2vol', mass2vol(i))
call setProperty(aTracer, 'vol2mass', vol2mass(i))
call setProperty(aTracer, 'dodrydep', dodrydep(i))
call setProperty(aTracer, 'F0', F0(i))
call setProperty(aTracer, 'HSTAR', HSTAR(i))
call setProperty(aTracer, 'do_fire', do_fire(i))
call setProperty(aTracer, 'nBBsources', nBBsources(i))
call setProperty(aTracer, 'emisPerFireByVegType', emisPerFireByVegType(i))
call setProperty(aTracer, 'trpdens', trpdens(i))
call setProperty(aTracer, 'trradius', trradius(i))
call setProperty(aTracer, 'tr_wd_TYPE', tr_wd_TYPE(i))
call setProperty(aTracer, 'tr_RKD', tr_RKD(i))
call setProperty(aTracer, 'tr_DHD', tr_DHD(i))
call setProperty(aTracer, 'fq_aer', fq_aer(i))
call setProperty(aTracer, 'rc_washt', rc_washt(i))
call setProperty(aTracer, 'isDust', isDust(i))
call setProperty(aTracer, 'tr_H2ObyCH4', tr_H2ObyCH4(i))
call setProperty(aTracer, 'dowetdep', dowetdep(i))
call setProperty(aTracer, 'ntrocn', ntrocn(i))
call setProperty(aTracer, 'conc_from_fw', conc_from_fw(i))
call setProperty(aTracer, 'trglac', trglac(i))
call setProperty(aTracer, 'ntisurfsrc', ntisurfsrc(i))
call setProperty(aTracer, 'TRLI0', TRLI0(i))

      call setProperty(aTracer, 'trw0', trw0(i))

      call setProperty(aTracer, 'TRSI0', TRSI0(i))


      call addTracer(bundle, aTracer)
      call clean(aTracer)

    end do
    
  end function makeNewTracers

    subroutine set_tr_mm(oldIndex, value)
    integer, intent(in) :: oldIndex
    real*8, intent(in) :: value
    tracers(oldIndex)%tr_mm = value
  end subroutine set_tr_mm
  
  function tr_mm(oldIndex)
    integer, intent(in) :: oldIndex
    real*8 :: tr_mm
    tr_mm = tracers(oldIndex)%tr_mm
  end function tr_mm

  subroutine set_ntm_power(oldIndex, value)
    integer, intent(in) :: oldIndex
    integer, intent(in) :: value
    tracers(oldIndex)%ntm_power = value
  end subroutine set_ntm_power
  
  function ntm_power(oldIndex)
    integer, intent(in) :: oldIndex
    integer :: ntm_power
    ntm_power = tracers(oldIndex)%ntm_power
  end function ntm_power

  subroutine set_t_qlimit(oldIndex, value)
    integer, intent(in) :: oldIndex
    logical, intent(in) :: value
    tracers(oldIndex)%t_qlimit = value
  end subroutine set_t_qlimit
  
  function t_qlimit(oldIndex)
    integer, intent(in) :: oldIndex
    logical :: t_qlimit
    t_qlimit = tracers(oldIndex)%t_qlimit
  end function t_qlimit

  subroutine set_ntsurfsrc(oldIndex, value)
    integer, intent(in) :: oldIndex
    integer, intent(in) :: value
    tracers(oldIndex)%ntsurfsrc = value
  end subroutine set_ntsurfsrc
  
  function ntsurfsrc(oldIndex)
    integer, intent(in) :: oldIndex
    integer :: ntsurfsrc
    ntsurfsrc = tracers(oldIndex)%ntsurfsrc
  end function ntsurfsrc

  subroutine set_needtrs(oldIndex, value)
    integer, intent(in) :: oldIndex
    logical, intent(in) :: value
    tracers(oldIndex)%needtrs = value
  end subroutine set_needtrs
  
  function needtrs(oldIndex)
    integer, intent(in) :: oldIndex
    logical :: needtrs
    needtrs = tracers(oldIndex)%needtrs
  end function needtrs

  subroutine set_trdecay(oldIndex, value)
    integer, intent(in) :: oldIndex
    real*8, intent(in) :: value
    tracers(oldIndex)%trdecay = value
  end subroutine set_trdecay
  
  function trdecay(oldIndex)
    integer, intent(in) :: oldIndex
    real*8 :: trdecay
    trdecay = tracers(oldIndex)%trdecay
  end function trdecay

  subroutine set_itime_tr0(oldIndex, value)
    integer, intent(in) :: oldIndex
    integer, intent(in) :: value
    tracers(oldIndex)%itime_tr0 = value
  end subroutine set_itime_tr0
  
  function itime_tr0(oldIndex)
    integer, intent(in) :: oldIndex
    integer :: itime_tr0
    itime_tr0 = tracers(oldIndex)%itime_tr0
  end function itime_tr0

  subroutine set_mass2vol(oldIndex, value)
    integer, intent(in) :: oldIndex
    real*8, intent(in) :: value
    tracers(oldIndex)%mass2vol = value
  end subroutine set_mass2vol
  
  function mass2vol(oldIndex)
    integer, intent(in) :: oldIndex
    real*8 :: mass2vol
    mass2vol = tracers(oldIndex)%mass2vol
  end function mass2vol

  subroutine set_vol2mass(oldIndex, value)
    integer, intent(in) :: oldIndex
    real*8, intent(in) :: value
    tracers(oldIndex)%vol2mass = value
  end subroutine set_vol2mass
  
  function vol2mass(oldIndex)
    integer, intent(in) :: oldIndex
    real*8 :: vol2mass
    vol2mass = tracers(oldIndex)%vol2mass
  end function vol2mass

  subroutine set_dodrydep(oldIndex, value)
    integer, intent(in) :: oldIndex
    logical, intent(in) :: value
    tracers(oldIndex)%dodrydep = value
  end subroutine set_dodrydep
  
  function dodrydep(oldIndex)
    integer, intent(in) :: oldIndex
    logical :: dodrydep
    dodrydep = tracers(oldIndex)%dodrydep
  end function dodrydep

  subroutine set_F0(oldIndex, value)
    integer, intent(in) :: oldIndex
    real*8, intent(in) :: value
    tracers(oldIndex)%F0 = value
  end subroutine set_F0
  
  function F0(oldIndex)
    integer, intent(in) :: oldIndex
    real*8 :: F0
    F0 = tracers(oldIndex)%F0
  end function F0

  subroutine set_HSTAR(oldIndex, value)
    integer, intent(in) :: oldIndex
    real*8, intent(in) :: value
    tracers(oldIndex)%HSTAR = value
  end subroutine set_HSTAR
  
  function HSTAR(oldIndex)
    integer, intent(in) :: oldIndex
    real*8 :: HSTAR
    HSTAR = tracers(oldIndex)%HSTAR
  end function HSTAR

  subroutine set_do_fire(oldIndex, value)
    integer, intent(in) :: oldIndex
    logical, intent(in) :: value
    tracers(oldIndex)%do_fire = value
  end subroutine set_do_fire
  
  function do_fire(oldIndex)
    integer, intent(in) :: oldIndex
    logical :: do_fire
    do_fire = tracers(oldIndex)%do_fire
  end function do_fire

  subroutine set_nBBsources(oldIndex, value)
    integer, intent(in) :: oldIndex
    integer, intent(in) :: value
    tracers(oldIndex)%nBBsources = value
  end subroutine set_nBBsources
  
  function nBBsources(oldIndex)
    integer, intent(in) :: oldIndex
    integer :: nBBsources
    nBBsources = tracers(oldIndex)%nBBsources
  end function nBBsources

  subroutine set_emisPerFireByVegType(oldIndex, value)
    integer, intent(in) :: oldIndex
    real*8, dimension(12), intent(in) :: value
    tracers(oldIndex)%emisPerFireByVegType = value
  end subroutine set_emisPerFireByVegType
  
  function emisPerFireByVegType(oldIndex)
    integer, intent(in) :: oldIndex
    real*8, dimension(12) :: emisPerFireByVegType
    emisPerFireByVegType = tracers(oldIndex)%emisPerFireByVegType
  end function emisPerFireByVegType

  subroutine set_trpdens(oldIndex, value)
    integer, intent(in) :: oldIndex
    real*8, intent(in) :: value
    tracers(oldIndex)%trpdens = value
  end subroutine set_trpdens
  
  function trpdens(oldIndex)
    integer, intent(in) :: oldIndex
    real*8 :: trpdens
    trpdens = tracers(oldIndex)%trpdens
  end function trpdens

  subroutine set_trradius(oldIndex, value)
    integer, intent(in) :: oldIndex
    real*8, intent(in) :: value
    tracers(oldIndex)%trradius = value
  end subroutine set_trradius
  
  function trradius(oldIndex)
    integer, intent(in) :: oldIndex
    real*8 :: trradius
    trradius = tracers(oldIndex)%trradius
  end function trradius

  subroutine set_tr_wd_TYPE(oldIndex, value)
    integer, intent(in) :: oldIndex
    integer, intent(in) :: value
    tracers(oldIndex)%tr_wd_TYPE = value
  end subroutine set_tr_wd_TYPE
  
  function tr_wd_TYPE(oldIndex)
    integer, intent(in) :: oldIndex
    integer :: tr_wd_TYPE
    tr_wd_TYPE = tracers(oldIndex)%tr_wd_TYPE
  end function tr_wd_TYPE

  subroutine set_tr_RKD(oldIndex, value)
    integer, intent(in) :: oldIndex
    real*8, intent(in) :: value
    tracers(oldIndex)%tr_RKD = value
  end subroutine set_tr_RKD
  
  function tr_RKD(oldIndex)
    integer, intent(in) :: oldIndex
    real*8 :: tr_RKD
    tr_RKD = tracers(oldIndex)%tr_RKD
  end function tr_RKD

  subroutine set_tr_DHD(oldIndex, value)
    integer, intent(in) :: oldIndex
    real*8, intent(in) :: value
    tracers(oldIndex)%tr_DHD = value
  end subroutine set_tr_DHD
  
  function tr_DHD(oldIndex)
    integer, intent(in) :: oldIndex
    real*8 :: tr_DHD
    tr_DHD = tracers(oldIndex)%tr_DHD
  end function tr_DHD

  subroutine set_fq_aer(oldIndex, value)
    integer, intent(in) :: oldIndex
    real*8, intent(in) :: value
    tracers(oldIndex)%fq_aer = value
  end subroutine set_fq_aer
  
  function fq_aer(oldIndex)
    integer, intent(in) :: oldIndex
    real*8 :: fq_aer
    fq_aer = tracers(oldIndex)%fq_aer
  end function fq_aer

  subroutine set_rc_washt(oldIndex, value)
    integer, intent(in) :: oldIndex
    real*8, intent(in) :: value
    tracers(oldIndex)%rc_washt = value
  end subroutine set_rc_washt
  
  function rc_washt(oldIndex)
    integer, intent(in) :: oldIndex
    real*8 :: rc_washt
    rc_washt = tracers(oldIndex)%rc_washt
  end function rc_washt

  subroutine set_isDust(oldIndex, value)
    integer, intent(in) :: oldIndex
    integer, intent(in) :: value
    tracers(oldIndex)%isDust = value
  end subroutine set_isDust
  
  function isDust(oldIndex)
    integer, intent(in) :: oldIndex
    integer :: isDust
    isDust = tracers(oldIndex)%isDust
  end function isDust

  subroutine set_tr_H2ObyCH4(oldIndex, value)
    integer, intent(in) :: oldIndex
    real*8, intent(in) :: value
    tracers(oldIndex)%tr_H2ObyCH4 = value
  end subroutine set_tr_H2ObyCH4
  
  function tr_H2ObyCH4(oldIndex)
    integer, intent(in) :: oldIndex
    real*8 :: tr_H2ObyCH4
    tr_H2ObyCH4 = tracers(oldIndex)%tr_H2ObyCH4
  end function tr_H2ObyCH4

  subroutine set_dowetdep(oldIndex, value)
    integer, intent(in) :: oldIndex
    logical, intent(in) :: value
    tracers(oldIndex)%dowetdep = value
  end subroutine set_dowetdep
  
  function dowetdep(oldIndex)
    integer, intent(in) :: oldIndex
    logical :: dowetdep
    dowetdep = tracers(oldIndex)%dowetdep
  end function dowetdep

  subroutine set_ntrocn(oldIndex, value)
    integer, intent(in) :: oldIndex
    integer, intent(in) :: value
    tracers(oldIndex)%ntrocn = value
  end subroutine set_ntrocn
  
  function ntrocn(oldIndex)
    integer, intent(in) :: oldIndex
    integer :: ntrocn
    ntrocn = tracers(oldIndex)%ntrocn
  end function ntrocn

  subroutine set_conc_from_fw(oldIndex, value)
    integer, intent(in) :: oldIndex
    logical, intent(in) :: value
    tracers(oldIndex)%conc_from_fw = value
  end subroutine set_conc_from_fw
  
  function conc_from_fw(oldIndex)
    integer, intent(in) :: oldIndex
    logical :: conc_from_fw
    conc_from_fw = tracers(oldIndex)%conc_from_fw
  end function conc_from_fw

  subroutine set_trglac(oldIndex, value)
    integer, intent(in) :: oldIndex
    real*8, intent(in) :: value
    tracers(oldIndex)%trglac = value
  end subroutine set_trglac
  
  function trglac(oldIndex)
    integer, intent(in) :: oldIndex
    real*8 :: trglac
    trglac = tracers(oldIndex)%trglac
  end function trglac

  subroutine set_ntisurfsrc(oldIndex, value)
    integer, intent(in) :: oldIndex
    integer, intent(in) :: value
    tracers(oldIndex)%ntisurfsrc = value
  end subroutine set_ntisurfsrc
  
  function ntisurfsrc(oldIndex)
    integer, intent(in) :: oldIndex
    integer :: ntisurfsrc
    ntisurfsrc = tracers(oldIndex)%ntisurfsrc
  end function ntisurfsrc

  subroutine set_TRLI0(oldIndex, value)
    integer, intent(in) :: oldIndex
    real*8, intent(in) :: value
    tracers(oldIndex)%TRLI0 = value
  end subroutine set_TRLI0
  
  function TRLI0(oldIndex)
    integer, intent(in) :: oldIndex
    real*8 :: TRLI0
    TRLI0 = tracers(oldIndex)%TRLI0
  end function TRLI0



subroutine set_trw0(oldIndex, value)
  integer, intent(in) :: oldIndex
  real*8, intent(in) :: value
  tracers(oldIndex)%trw0 = value
end subroutine set_trw0

function trw0_all()
  real*8 :: trw0_all(numTracers)
  trw0_all = tracers(:)%trw0
end function trw0_all

function trw0_one(n)
  integer, intent(in) :: n
  real*8 :: trw0_one
  trw0_one = tracers(n)%trw0
end function trw0_one

subroutine set_trsi0(oldIndex, value)
  integer, intent(in) :: oldIndex
  real*8, intent(in) :: value
  tracers(oldIndex)%trsi0 = value
end subroutine set_trsi0

function trsi0_all()
  real*8 :: trsi0_all(numTracers)
  trsi0_all = tracers(:)%trsi0
end function trsi0_all

function trsi0_one(n)
  integer, intent(in) :: n
  real*8 :: trsi0_one
  trsi0_one = tracers(n)%trsi0
end function trsi0_one

end module OldTracer_mod
