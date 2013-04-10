module TestTracerBundle_mod
  use pFUnit
  use Tracer_mod
  use TracerBundle_mod

  implicit none

  public :: setUp
  public :: tearDown
  public :: fixture

  type fixture
    type (Tracer), pointer :: speciesA
    type (Tracer), pointer :: speciesB
    type (TracerBundle) :: bundle
    integer :: numTracers
  end type fixture

contains

  subroutine setUp(this)
    use FileManager
    use GenericType_mod
    type (fixture) :: this

    type (Tracer), pointer :: speciesA, speciesB
    character(len=100) :: name

    this%speciesA => newTracer('speciesA')
    call this%speciesA%insert('molecularMass', 1.2d+0)
    call this%speciesA%insert('optionA', .true.)
    call this%speciesA%insert('optionC', .false.)
    call this%speciesA%insert('optionD', 1)

    this%speciesB => newTracer('speciesB')
    call this%speciesB%insert('molecularMass', 2.3d+0)
    call this%speciesB%insert('optionB', .false.)
    call this%speciesB%insert('optionD', 2)

    this%bundle = newTracerBundle()
    name = getName(this%speciesA)
    call this%bundle%insert(getName(this%speciesA), this%speciesA)
    call this%bundle%insert(getName(this%speciesB), this%speciesB)
    this%numTracers = 2

  end subroutine setUp

  subroutine tearDown(this)
    type (fixture) :: this

!!$    call clean(this%speciesA)
!!$    call clean(this%speciesB)
!!$    deallocate(this%speciesA)
!!$    deallocate(this%speciesB)

!!$    call clean(this%bundle)

  end subroutine tearDown

  subroutine testGetAttribute(this)
    use GenericType_mod
    type (fixture) :: this

    real(kind=r64) :: expected
    real(kind=r64), pointer :: found

    expected = 1.2d+0
    found = this%bundle%getAttribute('speciesA', 'molecularMass')
    call assertEqual(expected, found)

  end subroutine testGetAttribute

  subroutine testHasAttribute(this)
    use GenericType_mod
    type (fixture) :: this

    real(kind=r64) :: expected

    call assertTrue(all(this%bundle%hasAttribute('molecularMass')),'molecularMass')
    call assertTrue(all([.true.,.false.] .eqv. this%bundle%hasAttribute('optionA')),'optionA')

  end subroutine testHasAttribute

!!$  subroutine testGetSubsetA(this)
!!$    use GenericType_mod
!!$    use KeyValuePair_mod
!!$    type (fixture) :: this
!!$
!!$    type (TracerBundle) :: subset
!!$    class (Tracer), pointer :: t
!!$
!!$    print*,__LINE__,__FILE__
!!$    subset = this%bundle%makeSubset(withAttribute='optionA')
!!$    print*,__LINE__,__FILE__
!!$    call assertEqual(1, subset%size())
!!$    print*,__LINE__,__FILE__
!!$    t => subset%getReference('speciesA')
!!$    print*,__LINE__,__FILE__
!!$    call assertTrue(associated(t))
!!$    print*,__LINE__,__FILE__
!!$    call clean(subset)
!!$    print*,__LINE__,__FILE__
!!$
!!$  end subroutine testGetSubsetA
!!$
!!$  subroutine testGetSubsetB(this)
!!$    use GenericType_mod
!!$    type (fixture) :: this
!!$
!!$    type (TracerBundle) :: subset
!!$
!!$    subset = this%bundle%makeSubset(withAttribute='optionD')
!!$    call assertEqual(2, subset%size())
!!$    call clean(subset)
!!$
!!$  end subroutine testGetSubsetB

  subroutine testAddMandatoryA(this)
    type (fixture) :: this

    call this%bundle%addMandatoryAttribute('molecularMass')
    call this%bundle%insert('otherTracer', newTracer('otherTracer'))
    call assertFailedAssert( &
         & "TracerBundle_mod - species 'otherTracer' is missing mandatory attribute 'molecularMass'.", &
         & "Failed to detect missing mandatory option.")

  end subroutine testAddMandatoryA

  subroutine testAddMandatoryB(this)
    type (fixture) :: this

    call this%bundle%addMandatoryAttribute('newAttribute')

    call assertFailedAssert( &
         & "TracerBundle_mod - species 'speciesA' is missing mandatory attribute 'newAttribute'.", &
         & "Failed to detect missing mandatory option.")

  end subroutine testAddMandatoryB

  subroutine testDefaultValues(this)
    use GenericType_mod
    type (fixture) :: this
    logical :: expectedOptionC
    logical, pointer :: foundOptionC
    integer :: expectedOptionD
    integer, pointer :: foundOptionD

    call this%bundle%addDefault('optionC', .true.)
    call this%bundle%addDefault('optionD',  1)

    foundOptionC = this%bundle%getAttribute('speciesA', 'optionC')
    expectedOptionC = .false.
    call assertTrue(expectedOptionC .eqv. foundOptionC, 'first')

    foundOptionC = this%bundle%getAttribute('speciesB', 'optionC')
    expectedOptionC = .true.
    call assertTrue(expectedOptionC .eqv. foundOptionC, 'second')

    foundOptionD = this%bundle%getAttribute('speciesA', 'optionD')
    expectedOptionD = 1
    call assertTrue(expectedOptionD == foundOptionD, 'third')
    
  end subroutine testDefaultValues

! Check that a new tracers have defaults applied
  subroutine testDefaultAddNewTracer(this)
    use GenericType_mod
    type (fixture) :: this
    
    integer, parameter :: EXPECTED_VALUE = 1
    integer, parameter :: OTHER_VALUE = 2
    character(len=*), parameter :: NAME_1 = 'Tracer_1'
    character(len=*), parameter :: NAME_2 = 'Tracer_2'
    character(len=*), parameter :: ATTRIBUTE = 'newAttribute'
    integer, pointer :: foundValue

    type (Tracer), pointer :: aTracer
    type (TracerBundle) :: aBundle

!TODO - this is only needed because of unbalanced setUp/tearDown
    this%bundle = newTracerBundle()

    aBundle = newTracerBundle()
    call aBundle%addDefault(ATTRIBUTE, EXPECTED_VALUE)

    call aBundle%insert(NAME_1, newTracer(NAME_1))

    aTracer => newTracer(NAME_2)
    call aTracer%insert(ATTRIBUTE, OTHER_VALUE)
    call aBundle%insert(NAME_2, aTracer)

    foundValue = aBundle%getAttribute(NAME_1, ATTRIBUTE)
    call assertEqual(EXPECTED_VALUE, foundValue, 'default should be applied')

    foundValue = aBundle%getAttribute(NAME_2, ATTRIBUTE)
    call assertEqual(OTHER_VALUE, foundValue, 'default value should not override') 

    call clean(aTracer)
    call clean(aBundle)

  end subroutine testDefaultAddNewTracer

! Check that a new default value is applied to existing tracers.
  subroutine testDefaultApplyOld(this)
    use GenericType_mod
    type (fixture) :: this
    
    integer, parameter :: EXPECTED_VALUE = 1
    integer, parameter :: OTHER_VALUE = 2
    character(len=*), parameter :: NAME_1 = 'Tracer_1'
    character(len=*), parameter :: NAME_2 = 'Tracer_2'
    character(len=*), parameter :: ATTRIBUTE = 'newAttribute'
    integer, pointer :: foundValue

    type (Tracer), pointer :: aTracer
    type (TracerBundle) :: aBundle

!TODO - this is only needed because of unbalanced setUp/tearDown
    this%bundle = newTracerBundle()

    aBundle = newTracerBundle()
    call aBundle%insert(NAME_1, newTracer(NAME_1))

    aTracer => newTracer(NAME_2)
    call aTracer%insert(ATTRIBUTE, OTHER_VALUE)
    call aBundle%insert(NAME_2, aTracer)

    call aBundle%addDefault(ATTRIBUTE, EXPECTED_VALUE)

    foundValue = aBundle%getAttribute(NAME_1, ATTRIBUTE)
    call assertEqual(EXPECTED_VALUE, foundValue, 'should obtain default')

    foundValue = aBundle%getAttribute(NAME_2, ATTRIBUTE)
    call assertEqual(OTHER_VALUE, foundValue, 'should not override existing value')

    call clean(aTracer)
    call clean(aBundle)

  end subroutine testDefaultApplyOld

!!$  subroutine testLockSubset(this)
!!$    type (Fixture) :: this
!!$
!!$    type (TracerBundle) :: subset
!!$    class (Tracer), pointer :: t
!!$    character(len=*), parameter :: NAME_1 = 'Tracer_1'
!!$
!!$    subset = this%bundle%makeSubset(withAttribute='optionA')
!!$    call subset%insert(NAME_1, newTracer(NAME_1))
!!$
!!$    call assertFailedAssert( &
!!$         & "TracerBundle_mod - cannot insert new tracer into subset. "// &
!!$         & "Subsets are locked from modification.")
!!$
!!$  end subroutine testLockSubset

  subroutine testGetAttributeVector(this)
    use AttributeReference_mod
    use Attributes_mod
    type (fixture) :: this

    type (AttributeReference), pointer :: molecularMass(:)
    real*8, pointer :: mm1, mm2
    real*8, pointer :: mm1Found, mm2Found

    molecularMass => null()
    molecularMass => this%bundle%getAttributeVector('molecularMass')
    molecularMass => this%bundle%getAttributeVector('molecularMass')
    call assertTrue(associated(molecularMass))

    call assertEqual(2, size(molecularMass))
    mm1 = this%bundle%findAttribute('speciesA', 'molecularMass')
    mm1Found = molecularMass(1)

    mm2 = this%bundle%findAttribute('speciesB', 'molecularMass')
    mm2Found = molecularMass(2)
    
    call assertEqual(mm1, mm1Found)
    call assertEqual(mm2, mm2Found)

  end subroutine testGetAttributeVector

  ! The following test demonstrated a latent bug Intel's move_alloc used in AssociativeArray.
  subroutine test_modelErepdroducer(this)
    use Attributes_mod
    type (fixture) :: this

    class(Tracer), pointer :: t
    type (TracerBundle) :: bundle

    bundle = newTracerBundle()
    
    print*,__LINE__,__FILE__
    t=> newTracer('Ox')
    call bundle%insert(t)
    t => bundle%getReference('Ox')
    call t%insert('index', 1)

    print*,__LINE__,__FILE__
    t=> newTracer('NOx')
    call bundle%insert(t)
    t => bundle%getReference('NOx')
    call t%insert('index', 2)

    print*,__LINE__,__FILE__
    t=> newTracer('ClOx')
    call bundle%insert(t)
    t => bundle%getReference('ClOx')
    call t%insert('index', 3)

    t=> newTracer('BrOx')
    call bundle%insert(t)

    t => bundle%getReference('BrOx')
    call t%insert('index', 4)

    t=> newTracer('N2O5')
    call bundle%insert(t)

    print*,__LINE__,__FILE__
    t=> newTracer('DMS')
    call bundle%insert(t)
    t => bundle%getReference('DMS')
    call t%insert('index', 6)

    print*,__LINE__,__FILE__
    call bundle%setAttribute('N2O5', "ntm_power", newAttribute(-12))
    print*,__LINE__,__FILE__
    call bundle%setAttribute('N2O5', "tr_mm", newAttribute(108.02d0))
    print*,__LINE__,__FILE__
    call bundle%setAttribute('N2O5', "mass2vol", newAttribute(1.234d0))
    print*,__LINE__,__FILE__

  end subroutine test_modelErepdroducer

end module TestTracerBundle_mod
