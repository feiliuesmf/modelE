




























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
public :: set_trsi0, trsi0 
public :: set_trw0, trw0 
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

  interface tr_mm
module procedure tr_mm_s
module procedure tr_mm_all
   module procedure tr_mm_m
end interface
interface ntm_power
module procedure ntm_power_s
module procedure ntm_power_all
   module procedure ntm_power_m
end interface
interface t_qlimit
module procedure t_qlimit_s
module procedure t_qlimit_all
   module procedure t_qlimit_m
end interface
interface ntsurfsrc
module procedure ntsurfsrc_s
module procedure ntsurfsrc_all
   module procedure ntsurfsrc_m
end interface
interface needtrs
module procedure needtrs_s
module procedure needtrs_all
   module procedure needtrs_m
end interface
interface trdecay
module procedure trdecay_s
module procedure trdecay_all
   module procedure trdecay_m
end interface
interface itime_tr0
module procedure itime_tr0_s
module procedure itime_tr0_all
   module procedure itime_tr0_m
end interface
interface trsi0
module procedure trsi0_s
module procedure trsi0_all
   module procedure trsi0_m
end interface
interface trw0
module procedure trw0_s
module procedure trw0_all
   module procedure trw0_m
end interface
interface mass2vol
module procedure mass2vol_s
module procedure mass2vol_all
   module procedure mass2vol_m
end interface
interface vol2mass
module procedure vol2mass_s
module procedure vol2mass_all
   module procedure vol2mass_m
end interface
interface dodrydep
module procedure dodrydep_s
module procedure dodrydep_all
   module procedure dodrydep_m
end interface
interface F0
module procedure F0_s
module procedure F0_all
   module procedure F0_m
end interface
interface HSTAR
module procedure HSTAR_s
module procedure HSTAR_all
   module procedure HSTAR_m
end interface
interface do_fire
module procedure do_fire_s
module procedure do_fire_all
   module procedure do_fire_m
end interface
interface nBBsources
module procedure nBBsources_s
module procedure nBBsources_all
   module procedure nBBsources_m
end interface
interface emisPerFireByVegType
module procedure emisPerFireByVegType_s

end interface
interface trpdens
module procedure trpdens_s
module procedure trpdens_all
   module procedure trpdens_m
end interface
interface trradius
module procedure trradius_s
module procedure trradius_all
   module procedure trradius_m
end interface
interface tr_wd_TYPE
module procedure tr_wd_TYPE_s
module procedure tr_wd_TYPE_all
   module procedure tr_wd_TYPE_m
end interface
interface tr_RKD
module procedure tr_RKD_s
module procedure tr_RKD_all
   module procedure tr_RKD_m
end interface
interface tr_DHD
module procedure tr_DHD_s
module procedure tr_DHD_all
   module procedure tr_DHD_m
end interface
interface fq_aer
module procedure fq_aer_s
module procedure fq_aer_all
   module procedure fq_aer_m
end interface
interface rc_washt
module procedure rc_washt_s
module procedure rc_washt_all
   module procedure rc_washt_m
end interface
interface isDust
module procedure isDust_s
module procedure isDust_all
   module procedure isDust_m
end interface
interface tr_H2ObyCH4
module procedure tr_H2ObyCH4_s
module procedure tr_H2ObyCH4_all
   module procedure tr_H2ObyCH4_m
end interface
interface dowetdep
module procedure dowetdep_s
module procedure dowetdep_all
   module procedure dowetdep_m
end interface
interface ntrocn
module procedure ntrocn_s
module procedure ntrocn_all
   module procedure ntrocn_m
end interface
interface conc_from_fw
module procedure conc_from_fw_s
module procedure conc_from_fw_all
   module procedure conc_from_fw_m
end interface
interface trglac
module procedure trglac_s
module procedure trglac_all
   module procedure trglac_m
end interface
interface ntisurfsrc
module procedure ntisurfsrc_s
module procedure ntisurfsrc_all
   module procedure ntisurfsrc_m
end interface
interface TRLI0
module procedure TRLI0_s
module procedure TRLI0_all
   module procedure TRLI0_m
end interface


  integer, parameter :: MAX_LEN_NAME = 8
!**** parameters for tr_wd_TYPE
!@param nGAS   index for wetdep tracer type = gas
!@param nPART  index for wetdep tracer type = particle/aerosol
!@param nWATER index for wetdep tracer type = water
      integer, parameter :: nGAS=1, nPART=2, nWATER=3

  type OldTracer_type
    private
    character(len=MAX_LEN_NAME) :: name

!@var tr_mm: !@var TR_MM: molar mass of each tracer (g/mol)
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
!@var trsi0: conc. in sea ice (kg/m^2)
    real*8 :: trsi0 = 0.d0 
!@var trw0: concentration in water (kg/kg)
    real*8 :: trw0 = 0.d0 
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


  end type OldTracer_type

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
call setProperty(aTracer, 'trsi0', trsi0(i))
call setProperty(aTracer, 'trw0', trw0(i))
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


      call addTracer(bundle, aTracer)
      call clean(aTracer)

    end do
    
  end function makeNewTracers

    subroutine set_tr_mm(oldIndex, value)
    integer, intent(in) :: oldIndex
    real*8, intent(in) :: value
    tracers(oldIndex)%tr_mm = value
  end subroutine set_tr_mm
  
  function tr_mm_s(oldIndex)
    integer, intent(in) :: oldIndex
    real*8 :: tr_mm_s
    tr_mm_s = tracers(oldIndex)%tr_mm
  end function tr_mm_s

  function tr_mm_all()
    real*8 :: tr_mm_all(size(tracers))
    tr_mm_all = tracers(:)%tr_mm
  end function tr_mm_all

  function tr_mm_m(oldIndices)
    integer, intent(in) :: oldIndices(:)
    real*8 :: tr_mm_m(size(oldIndices))
    tr_mm_m = tracers(oldIndices(:))%tr_mm
  end function tr_mm_m



  subroutine set_ntm_power(oldIndex, value)
    integer, intent(in) :: oldIndex
    integer, intent(in) :: value
    tracers(oldIndex)%ntm_power = value
  end subroutine set_ntm_power
  
  function ntm_power_s(oldIndex)
    integer, intent(in) :: oldIndex
    integer :: ntm_power_s
    ntm_power_s = tracers(oldIndex)%ntm_power
  end function ntm_power_s

  function ntm_power_all()
    integer :: ntm_power_all(size(tracers))
    ntm_power_all = tracers(:)%ntm_power
  end function ntm_power_all

  function ntm_power_m(oldIndices)
    integer, intent(in) :: oldIndices(:)
    integer :: ntm_power_m(size(oldIndices))
    ntm_power_m = tracers(oldIndices(:))%ntm_power
  end function ntm_power_m



  subroutine set_t_qlimit(oldIndex, value)
    integer, intent(in) :: oldIndex
    logical, intent(in) :: value
    tracers(oldIndex)%t_qlimit = value
  end subroutine set_t_qlimit
  
  function t_qlimit_s(oldIndex)
    integer, intent(in) :: oldIndex
    logical :: t_qlimit_s
    t_qlimit_s = tracers(oldIndex)%t_qlimit
  end function t_qlimit_s

  function t_qlimit_all()
    logical :: t_qlimit_all(size(tracers))
    t_qlimit_all = tracers(:)%t_qlimit
  end function t_qlimit_all

  function t_qlimit_m(oldIndices)
    integer, intent(in) :: oldIndices(:)
    logical :: t_qlimit_m(size(oldIndices))
    t_qlimit_m = tracers(oldIndices(:))%t_qlimit
  end function t_qlimit_m



  subroutine set_ntsurfsrc(oldIndex, value)
    integer, intent(in) :: oldIndex
    integer, intent(in) :: value
    tracers(oldIndex)%ntsurfsrc = value
  end subroutine set_ntsurfsrc
  
  function ntsurfsrc_s(oldIndex)
    integer, intent(in) :: oldIndex
    integer :: ntsurfsrc_s
    ntsurfsrc_s = tracers(oldIndex)%ntsurfsrc
  end function ntsurfsrc_s

  function ntsurfsrc_all()
    integer :: ntsurfsrc_all(size(tracers))
    ntsurfsrc_all = tracers(:)%ntsurfsrc
  end function ntsurfsrc_all

  function ntsurfsrc_m(oldIndices)
    integer, intent(in) :: oldIndices(:)
    integer :: ntsurfsrc_m(size(oldIndices))
    ntsurfsrc_m = tracers(oldIndices(:))%ntsurfsrc
  end function ntsurfsrc_m



  subroutine set_needtrs(oldIndex, value)
    integer, intent(in) :: oldIndex
    logical, intent(in) :: value
    tracers(oldIndex)%needtrs = value
  end subroutine set_needtrs
  
  function needtrs_s(oldIndex)
    integer, intent(in) :: oldIndex
    logical :: needtrs_s
    needtrs_s = tracers(oldIndex)%needtrs
  end function needtrs_s

  function needtrs_all()
    logical :: needtrs_all(size(tracers))
    needtrs_all = tracers(:)%needtrs
  end function needtrs_all

  function needtrs_m(oldIndices)
    integer, intent(in) :: oldIndices(:)
    logical :: needtrs_m(size(oldIndices))
    needtrs_m = tracers(oldIndices(:))%needtrs
  end function needtrs_m



  subroutine set_trdecay(oldIndex, value)
    integer, intent(in) :: oldIndex
    real*8, intent(in) :: value
    tracers(oldIndex)%trdecay = value
  end subroutine set_trdecay
  
  function trdecay_s(oldIndex)
    integer, intent(in) :: oldIndex
    real*8 :: trdecay_s
    trdecay_s = tracers(oldIndex)%trdecay
  end function trdecay_s

  function trdecay_all()
    real*8 :: trdecay_all(size(tracers))
    trdecay_all = tracers(:)%trdecay
  end function trdecay_all

  function trdecay_m(oldIndices)
    integer, intent(in) :: oldIndices(:)
    real*8 :: trdecay_m(size(oldIndices))
    trdecay_m = tracers(oldIndices(:))%trdecay
  end function trdecay_m



  subroutine set_itime_tr0(oldIndex, value)
    integer, intent(in) :: oldIndex
    integer, intent(in) :: value
    tracers(oldIndex)%itime_tr0 = value
  end subroutine set_itime_tr0
  
  function itime_tr0_s(oldIndex)
    integer, intent(in) :: oldIndex
    integer :: itime_tr0_s
    itime_tr0_s = tracers(oldIndex)%itime_tr0
  end function itime_tr0_s

  function itime_tr0_all()
    integer :: itime_tr0_all(size(tracers))
    itime_tr0_all = tracers(:)%itime_tr0
  end function itime_tr0_all

  function itime_tr0_m(oldIndices)
    integer, intent(in) :: oldIndices(:)
    integer :: itime_tr0_m(size(oldIndices))
    itime_tr0_m = tracers(oldIndices(:))%itime_tr0
  end function itime_tr0_m



  subroutine set_trsi0(oldIndex, value)
    integer, intent(in) :: oldIndex
    real*8, intent(in) :: value
    tracers(oldIndex)%trsi0 = value
  end subroutine set_trsi0
  
  function trsi0_s(oldIndex)
    integer, intent(in) :: oldIndex
    real*8 :: trsi0_s
    trsi0_s = tracers(oldIndex)%trsi0
  end function trsi0_s

  function trsi0_all()
    real*8 :: trsi0_all(size(tracers))
    trsi0_all = tracers(:)%trsi0
  end function trsi0_all

  function trsi0_m(oldIndices)
    integer, intent(in) :: oldIndices(:)
    real*8 :: trsi0_m(size(oldIndices))
    trsi0_m = tracers(oldIndices(:))%trsi0
  end function trsi0_m



  subroutine set_trw0(oldIndex, value)
    integer, intent(in) :: oldIndex
    real*8, intent(in) :: value
    tracers(oldIndex)%trw0 = value
  end subroutine set_trw0
  
  function trw0_s(oldIndex)
    integer, intent(in) :: oldIndex
    real*8 :: trw0_s
    trw0_s = tracers(oldIndex)%trw0
  end function trw0_s

  function trw0_all()
    real*8 :: trw0_all(size(tracers))
    trw0_all = tracers(:)%trw0
  end function trw0_all

  function trw0_m(oldIndices)
    integer, intent(in) :: oldIndices(:)
    real*8 :: trw0_m(size(oldIndices))
    trw0_m = tracers(oldIndices(:))%trw0
  end function trw0_m



  subroutine set_mass2vol(oldIndex, value)
    integer, intent(in) :: oldIndex
    real*8, intent(in) :: value
    tracers(oldIndex)%mass2vol = value
  end subroutine set_mass2vol
  
  function mass2vol_s(oldIndex)
    integer, intent(in) :: oldIndex
    real*8 :: mass2vol_s
    mass2vol_s = tracers(oldIndex)%mass2vol
  end function mass2vol_s

  function mass2vol_all()
    real*8 :: mass2vol_all(size(tracers))
    mass2vol_all = tracers(:)%mass2vol
  end function mass2vol_all

  function mass2vol_m(oldIndices)
    integer, intent(in) :: oldIndices(:)
    real*8 :: mass2vol_m(size(oldIndices))
    mass2vol_m = tracers(oldIndices(:))%mass2vol
  end function mass2vol_m



  subroutine set_vol2mass(oldIndex, value)
    integer, intent(in) :: oldIndex
    real*8, intent(in) :: value
    tracers(oldIndex)%vol2mass = value
  end subroutine set_vol2mass
  
  function vol2mass_s(oldIndex)
    integer, intent(in) :: oldIndex
    real*8 :: vol2mass_s
    vol2mass_s = tracers(oldIndex)%vol2mass
  end function vol2mass_s

  function vol2mass_all()
    real*8 :: vol2mass_all(size(tracers))
    vol2mass_all = tracers(:)%vol2mass
  end function vol2mass_all

  function vol2mass_m(oldIndices)
    integer, intent(in) :: oldIndices(:)
    real*8 :: vol2mass_m(size(oldIndices))
    vol2mass_m = tracers(oldIndices(:))%vol2mass
  end function vol2mass_m



  subroutine set_dodrydep(oldIndex, value)
    integer, intent(in) :: oldIndex
    logical, intent(in) :: value
    tracers(oldIndex)%dodrydep = value
  end subroutine set_dodrydep
  
  function dodrydep_s(oldIndex)
    integer, intent(in) :: oldIndex
    logical :: dodrydep_s
    dodrydep_s = tracers(oldIndex)%dodrydep
  end function dodrydep_s

  function dodrydep_all()
    logical :: dodrydep_all(size(tracers))
    dodrydep_all = tracers(:)%dodrydep
  end function dodrydep_all

  function dodrydep_m(oldIndices)
    integer, intent(in) :: oldIndices(:)
    logical :: dodrydep_m(size(oldIndices))
    dodrydep_m = tracers(oldIndices(:))%dodrydep
  end function dodrydep_m



  subroutine set_F0(oldIndex, value)
    integer, intent(in) :: oldIndex
    real*8, intent(in) :: value
    tracers(oldIndex)%F0 = value
  end subroutine set_F0
  
  function F0_s(oldIndex)
    integer, intent(in) :: oldIndex
    real*8 :: F0_s
    F0_s = tracers(oldIndex)%F0
  end function F0_s

  function F0_all()
    real*8 :: F0_all(size(tracers))
    F0_all = tracers(:)%F0
  end function F0_all

  function F0_m(oldIndices)
    integer, intent(in) :: oldIndices(:)
    real*8 :: F0_m(size(oldIndices))
    F0_m = tracers(oldIndices(:))%F0
  end function F0_m



  subroutine set_HSTAR(oldIndex, value)
    integer, intent(in) :: oldIndex
    real*8, intent(in) :: value
    tracers(oldIndex)%HSTAR = value
  end subroutine set_HSTAR
  
  function HSTAR_s(oldIndex)
    integer, intent(in) :: oldIndex
    real*8 :: HSTAR_s
    HSTAR_s = tracers(oldIndex)%HSTAR
  end function HSTAR_s

  function HSTAR_all()
    real*8 :: HSTAR_all(size(tracers))
    HSTAR_all = tracers(:)%HSTAR
  end function HSTAR_all

  function HSTAR_m(oldIndices)
    integer, intent(in) :: oldIndices(:)
    real*8 :: HSTAR_m(size(oldIndices))
    HSTAR_m = tracers(oldIndices(:))%HSTAR
  end function HSTAR_m



  subroutine set_do_fire(oldIndex, value)
    integer, intent(in) :: oldIndex
    logical, intent(in) :: value
    tracers(oldIndex)%do_fire = value
  end subroutine set_do_fire
  
  function do_fire_s(oldIndex)
    integer, intent(in) :: oldIndex
    logical :: do_fire_s
    do_fire_s = tracers(oldIndex)%do_fire
  end function do_fire_s

  function do_fire_all()
    logical :: do_fire_all(size(tracers))
    do_fire_all = tracers(:)%do_fire
  end function do_fire_all

  function do_fire_m(oldIndices)
    integer, intent(in) :: oldIndices(:)
    logical :: do_fire_m(size(oldIndices))
    do_fire_m = tracers(oldIndices(:))%do_fire
  end function do_fire_m



  subroutine set_nBBsources(oldIndex, value)
    integer, intent(in) :: oldIndex
    integer, intent(in) :: value
    tracers(oldIndex)%nBBsources = value
  end subroutine set_nBBsources
  
  function nBBsources_s(oldIndex)
    integer, intent(in) :: oldIndex
    integer :: nBBsources_s
    nBBsources_s = tracers(oldIndex)%nBBsources
  end function nBBsources_s

  function nBBsources_all()
    integer :: nBBsources_all(size(tracers))
    nBBsources_all = tracers(:)%nBBsources
  end function nBBsources_all

  function nBBsources_m(oldIndices)
    integer, intent(in) :: oldIndices(:)
    integer :: nBBsources_m(size(oldIndices))
    nBBsources_m = tracers(oldIndices(:))%nBBsources
  end function nBBsources_m



  subroutine set_emisPerFireByVegType(oldIndex, value)
    integer, intent(in) :: oldIndex
    real*8, dimension(12), intent(in) :: value
    tracers(oldIndex)%emisPerFireByVegType = value
  end subroutine set_emisPerFireByVegType
  
  function emisPerFireByVegType_s(oldIndex)
    integer, intent(in) :: oldIndex
    real*8, dimension(12) :: emisPerFireByVegType_s
    emisPerFireByVegType_s = tracers(oldIndex)%emisPerFireByVegType
  end function emisPerFireByVegType_s



  subroutine set_trpdens(oldIndex, value)
    integer, intent(in) :: oldIndex
    real*8, intent(in) :: value
    tracers(oldIndex)%trpdens = value
  end subroutine set_trpdens
  
  function trpdens_s(oldIndex)
    integer, intent(in) :: oldIndex
    real*8 :: trpdens_s
    trpdens_s = tracers(oldIndex)%trpdens
  end function trpdens_s

  function trpdens_all()
    real*8 :: trpdens_all(size(tracers))
    trpdens_all = tracers(:)%trpdens
  end function trpdens_all

  function trpdens_m(oldIndices)
    integer, intent(in) :: oldIndices(:)
    real*8 :: trpdens_m(size(oldIndices))
    trpdens_m = tracers(oldIndices(:))%trpdens
  end function trpdens_m



  subroutine set_trradius(oldIndex, value)
    integer, intent(in) :: oldIndex
    real*8, intent(in) :: value
    tracers(oldIndex)%trradius = value
  end subroutine set_trradius
  
  function trradius_s(oldIndex)
    integer, intent(in) :: oldIndex
    real*8 :: trradius_s
    trradius_s = tracers(oldIndex)%trradius
  end function trradius_s

  function trradius_all()
    real*8 :: trradius_all(size(tracers))
    trradius_all = tracers(:)%trradius
  end function trradius_all

  function trradius_m(oldIndices)
    integer, intent(in) :: oldIndices(:)
    real*8 :: trradius_m(size(oldIndices))
    trradius_m = tracers(oldIndices(:))%trradius
  end function trradius_m



  subroutine set_tr_wd_TYPE(oldIndex, value)
    integer, intent(in) :: oldIndex
    integer, intent(in) :: value
    tracers(oldIndex)%tr_wd_TYPE = value
  end subroutine set_tr_wd_TYPE
  
  function tr_wd_TYPE_s(oldIndex)
    integer, intent(in) :: oldIndex
    integer :: tr_wd_TYPE_s
    tr_wd_TYPE_s = tracers(oldIndex)%tr_wd_TYPE
  end function tr_wd_TYPE_s

  function tr_wd_TYPE_all()
    integer :: tr_wd_TYPE_all(size(tracers))
    tr_wd_TYPE_all = tracers(:)%tr_wd_TYPE
  end function tr_wd_TYPE_all

  function tr_wd_TYPE_m(oldIndices)
    integer, intent(in) :: oldIndices(:)
    integer :: tr_wd_TYPE_m(size(oldIndices))
    tr_wd_TYPE_m = tracers(oldIndices(:))%tr_wd_TYPE
  end function tr_wd_TYPE_m



  subroutine set_tr_RKD(oldIndex, value)
    integer, intent(in) :: oldIndex
    real*8, intent(in) :: value
    tracers(oldIndex)%tr_RKD = value
  end subroutine set_tr_RKD
  
  function tr_RKD_s(oldIndex)
    integer, intent(in) :: oldIndex
    real*8 :: tr_RKD_s
    tr_RKD_s = tracers(oldIndex)%tr_RKD
  end function tr_RKD_s

  function tr_RKD_all()
    real*8 :: tr_RKD_all(size(tracers))
    tr_RKD_all = tracers(:)%tr_RKD
  end function tr_RKD_all

  function tr_RKD_m(oldIndices)
    integer, intent(in) :: oldIndices(:)
    real*8 :: tr_RKD_m(size(oldIndices))
    tr_RKD_m = tracers(oldIndices(:))%tr_RKD
  end function tr_RKD_m



  subroutine set_tr_DHD(oldIndex, value)
    integer, intent(in) :: oldIndex
    real*8, intent(in) :: value
    tracers(oldIndex)%tr_DHD = value
  end subroutine set_tr_DHD
  
  function tr_DHD_s(oldIndex)
    integer, intent(in) :: oldIndex
    real*8 :: tr_DHD_s
    tr_DHD_s = tracers(oldIndex)%tr_DHD
  end function tr_DHD_s

  function tr_DHD_all()
    real*8 :: tr_DHD_all(size(tracers))
    tr_DHD_all = tracers(:)%tr_DHD
  end function tr_DHD_all

  function tr_DHD_m(oldIndices)
    integer, intent(in) :: oldIndices(:)
    real*8 :: tr_DHD_m(size(oldIndices))
    tr_DHD_m = tracers(oldIndices(:))%tr_DHD
  end function tr_DHD_m



  subroutine set_fq_aer(oldIndex, value)
    integer, intent(in) :: oldIndex
    real*8, intent(in) :: value
    tracers(oldIndex)%fq_aer = value
  end subroutine set_fq_aer
  
  function fq_aer_s(oldIndex)
    integer, intent(in) :: oldIndex
    real*8 :: fq_aer_s
    fq_aer_s = tracers(oldIndex)%fq_aer
  end function fq_aer_s

  function fq_aer_all()
    real*8 :: fq_aer_all(size(tracers))
    fq_aer_all = tracers(:)%fq_aer
  end function fq_aer_all

  function fq_aer_m(oldIndices)
    integer, intent(in) :: oldIndices(:)
    real*8 :: fq_aer_m(size(oldIndices))
    fq_aer_m = tracers(oldIndices(:))%fq_aer
  end function fq_aer_m



  subroutine set_rc_washt(oldIndex, value)
    integer, intent(in) :: oldIndex
    real*8, intent(in) :: value
    tracers(oldIndex)%rc_washt = value
  end subroutine set_rc_washt
  
  function rc_washt_s(oldIndex)
    integer, intent(in) :: oldIndex
    real*8 :: rc_washt_s
    rc_washt_s = tracers(oldIndex)%rc_washt
  end function rc_washt_s

  function rc_washt_all()
    real*8 :: rc_washt_all(size(tracers))
    rc_washt_all = tracers(:)%rc_washt
  end function rc_washt_all

  function rc_washt_m(oldIndices)
    integer, intent(in) :: oldIndices(:)
    real*8 :: rc_washt_m(size(oldIndices))
    rc_washt_m = tracers(oldIndices(:))%rc_washt
  end function rc_washt_m



  subroutine set_isDust(oldIndex, value)
    integer, intent(in) :: oldIndex
    integer, intent(in) :: value
    tracers(oldIndex)%isDust = value
  end subroutine set_isDust
  
  function isDust_s(oldIndex)
    integer, intent(in) :: oldIndex
    integer :: isDust_s
    isDust_s = tracers(oldIndex)%isDust
  end function isDust_s

  function isDust_all()
    integer :: isDust_all(size(tracers))
    isDust_all = tracers(:)%isDust
  end function isDust_all

  function isDust_m(oldIndices)
    integer, intent(in) :: oldIndices(:)
    integer :: isDust_m(size(oldIndices))
    isDust_m = tracers(oldIndices(:))%isDust
  end function isDust_m



  subroutine set_tr_H2ObyCH4(oldIndex, value)
    integer, intent(in) :: oldIndex
    real*8, intent(in) :: value
    tracers(oldIndex)%tr_H2ObyCH4 = value
  end subroutine set_tr_H2ObyCH4
  
  function tr_H2ObyCH4_s(oldIndex)
    integer, intent(in) :: oldIndex
    real*8 :: tr_H2ObyCH4_s
    tr_H2ObyCH4_s = tracers(oldIndex)%tr_H2ObyCH4
  end function tr_H2ObyCH4_s

  function tr_H2ObyCH4_all()
    real*8 :: tr_H2ObyCH4_all(size(tracers))
    tr_H2ObyCH4_all = tracers(:)%tr_H2ObyCH4
  end function tr_H2ObyCH4_all

  function tr_H2ObyCH4_m(oldIndices)
    integer, intent(in) :: oldIndices(:)
    real*8 :: tr_H2ObyCH4_m(size(oldIndices))
    tr_H2ObyCH4_m = tracers(oldIndices(:))%tr_H2ObyCH4
  end function tr_H2ObyCH4_m



  subroutine set_dowetdep(oldIndex, value)
    integer, intent(in) :: oldIndex
    logical, intent(in) :: value
    tracers(oldIndex)%dowetdep = value
  end subroutine set_dowetdep
  
  function dowetdep_s(oldIndex)
    integer, intent(in) :: oldIndex
    logical :: dowetdep_s
    dowetdep_s = tracers(oldIndex)%dowetdep
  end function dowetdep_s

  function dowetdep_all()
    logical :: dowetdep_all(size(tracers))
    dowetdep_all = tracers(:)%dowetdep
  end function dowetdep_all

  function dowetdep_m(oldIndices)
    integer, intent(in) :: oldIndices(:)
    logical :: dowetdep_m(size(oldIndices))
    dowetdep_m = tracers(oldIndices(:))%dowetdep
  end function dowetdep_m



  subroutine set_ntrocn(oldIndex, value)
    integer, intent(in) :: oldIndex
    integer, intent(in) :: value
    tracers(oldIndex)%ntrocn = value
  end subroutine set_ntrocn
  
  function ntrocn_s(oldIndex)
    integer, intent(in) :: oldIndex
    integer :: ntrocn_s
    ntrocn_s = tracers(oldIndex)%ntrocn
  end function ntrocn_s

  function ntrocn_all()
    integer :: ntrocn_all(size(tracers))
    ntrocn_all = tracers(:)%ntrocn
  end function ntrocn_all

  function ntrocn_m(oldIndices)
    integer, intent(in) :: oldIndices(:)
    integer :: ntrocn_m(size(oldIndices))
    ntrocn_m = tracers(oldIndices(:))%ntrocn
  end function ntrocn_m



  subroutine set_conc_from_fw(oldIndex, value)
    integer, intent(in) :: oldIndex
    logical, intent(in) :: value
    tracers(oldIndex)%conc_from_fw = value
  end subroutine set_conc_from_fw
  
  function conc_from_fw_s(oldIndex)
    integer, intent(in) :: oldIndex
    logical :: conc_from_fw_s
    conc_from_fw_s = tracers(oldIndex)%conc_from_fw
  end function conc_from_fw_s

  function conc_from_fw_all()
    logical :: conc_from_fw_all(size(tracers))
    conc_from_fw_all = tracers(:)%conc_from_fw
  end function conc_from_fw_all

  function conc_from_fw_m(oldIndices)
    integer, intent(in) :: oldIndices(:)
    logical :: conc_from_fw_m(size(oldIndices))
    conc_from_fw_m = tracers(oldIndices(:))%conc_from_fw
  end function conc_from_fw_m



  subroutine set_trglac(oldIndex, value)
    integer, intent(in) :: oldIndex
    real*8, intent(in) :: value
    tracers(oldIndex)%trglac = value
  end subroutine set_trglac
  
  function trglac_s(oldIndex)
    integer, intent(in) :: oldIndex
    real*8 :: trglac_s
    trglac_s = tracers(oldIndex)%trglac
  end function trglac_s

  function trglac_all()
    real*8 :: trglac_all(size(tracers))
    trglac_all = tracers(:)%trglac
  end function trglac_all

  function trglac_m(oldIndices)
    integer, intent(in) :: oldIndices(:)
    real*8 :: trglac_m(size(oldIndices))
    trglac_m = tracers(oldIndices(:))%trglac
  end function trglac_m



  subroutine set_ntisurfsrc(oldIndex, value)
    integer, intent(in) :: oldIndex
    integer, intent(in) :: value
    tracers(oldIndex)%ntisurfsrc = value
  end subroutine set_ntisurfsrc
  
  function ntisurfsrc_s(oldIndex)
    integer, intent(in) :: oldIndex
    integer :: ntisurfsrc_s
    ntisurfsrc_s = tracers(oldIndex)%ntisurfsrc
  end function ntisurfsrc_s

  function ntisurfsrc_all()
    integer :: ntisurfsrc_all(size(tracers))
    ntisurfsrc_all = tracers(:)%ntisurfsrc
  end function ntisurfsrc_all

  function ntisurfsrc_m(oldIndices)
    integer, intent(in) :: oldIndices(:)
    integer :: ntisurfsrc_m(size(oldIndices))
    ntisurfsrc_m = tracers(oldIndices(:))%ntisurfsrc
  end function ntisurfsrc_m



  subroutine set_TRLI0(oldIndex, value)
    integer, intent(in) :: oldIndex
    real*8, intent(in) :: value
    tracers(oldIndex)%TRLI0 = value
  end subroutine set_TRLI0
  
  function TRLI0_s(oldIndex)
    integer, intent(in) :: oldIndex
    real*8 :: TRLI0_s
    TRLI0_s = tracers(oldIndex)%TRLI0
  end function TRLI0_s

  function TRLI0_all()
    real*8 :: TRLI0_all(size(tracers))
    TRLI0_all = tracers(:)%TRLI0
  end function TRLI0_all

  function TRLI0_m(oldIndices)
    integer, intent(in) :: oldIndices(:)
    real*8 :: TRLI0_m(size(oldIndices))
    TRLI0_m = tracers(oldIndices(:))%TRLI0
  end function TRLI0_m





end module OldTracer_mod
