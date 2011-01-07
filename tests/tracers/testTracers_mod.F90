module testTracers_mod
  use pFUnit
  use Tracers_mod
  use TracerBundle_mod
  use Dictionary_mod, only: Dictionary_type
  implicit none
  private

  public :: testGetIndex
  public :: testHasPropertySingle
  public :: testHasPropertyMulti
  public :: testSetProperty
  public :: testSetPropertyMulti
  public :: testGetPropertySingle
  public :: testGetPropertyMulti
  public :: testCountHasProperty
  public :: testGetSubsetA
  public :: testGetSubsetB
  public :: testCheckSubsetDefault
  public :: testAddMandatoryA
  public :: testAddMandatoryB
  public :: testAddMandatoryC
  public :: testDefaultValues
  public :: testDefaultAddNewTracer
  public :: testDefaultApplyOld

  public :: setUp
  public :: tearDown
  public :: fixture


  type fixture
    type (TracerBundle_type) :: bundle
    integer :: numTracers
  end type fixture

contains

  subroutine setUp(this)
    use Dictionary_mod, only: Dictionary
    use Dictionary_mod, only: insert
    use FileManager
    type (fixture) :: this

    type (Tracer_type) :: speciesA, speciesB
    
    speciesA = Tracer('speciesA')
    call setProperty(speciesA, 'molecularMass', 1.2d+0)
    call setProperty(speciesA, 'optionA', .true.)
    call setProperty(speciesA, 'optionC', .false.)
    call setProperty(speciesA, 'optionD', 1)

    speciesB = Tracer('speciesB')
    call setProperty(speciesB, 'molecularMass', 2.3d+0)
    call setProperty(speciesB, 'optionB', .false.)
    call setProperty(speciesB, 'optionD', 2)

    this%bundle = TracerBundle()
    call addTracer(this%bundle, speciesA)
    call addTracer(this%bundle, speciesB)
    this%numTracers = 2

    call clean(speciesA)
    call clean(speciesB)

  end subroutine setUp

  subroutine tearDown(this)
    use Dictionary_mod, only: clean
    type (fixture) :: this

    call clean(this%bundle)

  end subroutine tearDown

  subroutine testGetIndex(this)
    type (fixture) :: this

    call assertEqual(1, getIndex(this%bundle, 'speciesA'))
    call assertEqual(2, getIndex(this%bundle, 'speciesB'))
    call assertEqual(NOT_FOUND, getIndex(this%bundle, 'speciesZ'))

  end subroutine testGetIndex

  subroutine testGetPropertySingle(this)
    use GenericType_mod
    type (fixture) :: this

    real(kind=r64) :: expected, found
    type (Tracer_type), pointer :: tracer

    expected = 1.2d+0
    tracer => getTracer(this%bundle, 'speciesA')
    found = getProperty(tracer, 'molecularMass')
    call assertEqual(expected, found)

  end subroutine testGetPropertySingle

  subroutine testGetPropertyMulti(this)
    use GenericType_mod
    type (fixture) :: this

    real(kind=r64) :: expected, found

    expected = 1.2d+0
    found = getProperty(this%bundle, 'speciesA', 'molecularMass')
    call assertEqual(expected, found)

  end subroutine testGetPropertyMulti

  subroutine testSetProperty(this)
    use GenericType_mod
    use Dictionary_mod
    type (fixture) :: this

    real(kind=r64) :: expected, found
    type (Dictionary_type), pointer :: properties
    type (Tracer_type), pointer :: tracer

    tracer => getTracer(this%bundle, 'speciesA')

    call setProperty(tracer, 'otherProperty', expected)
    found = getProperty(tracer, 'otherProperty')
    call assertEqual(expected, found)

  end subroutine testSetProperty

  subroutine testSetPropertyMulti(this)
    use GenericType_mod
    use Dictionary_mod
    type (fixture) :: this

    real(kind=r64) :: expected, found
    type (Dictionary_type), pointer :: properties
    type (Tracer_type), pointer :: tracer


    expected = 2.1d+0
    tracer => getTracer(this%bundle, 'speciesB')
    call setProperty(tracer, 'otherProperty', expected)

    found = getProperty(this%bundle, 'speciesB', 'otherProperty')
    call assertEqual(expected, found)

  end subroutine testSetPropertyMulti

  subroutine testHasPropertySingle(this)
    use GenericType_mod
    type (fixture) :: this

    real(kind=r64) :: expected, found
    type (Tracer_type), pointer :: tracer

    tracer => getTracer(this%bundle,'speciesA')

    call assertTrue(hasProperty(tracer, 'optionA'))
    call assertFalse(hasProperty(tracer, 'optionB'))

  end subroutine testHasPropertySingle

  subroutine testHasPropertyMulti(this)
    use GenericType_mod
    type (fixture) :: this

    real(kind=r64) :: expected, found


    call assertTrue(all(hasProperty(this%bundle, 'molecularMass')))
    call assertTrue(all([.true.,.false.] .eqv. hasProperty(this%bundle, 'optionA')))

  end subroutine testHasPropertyMulti

  subroutine testCountHasProperty(this)
    use GenericType_mod
    type (fixture) :: this

    call assertEqual(2, getCount(this%bundle, 'name'))
    call assertEqual(1, getCount(this%bundle, 'optionA'))

  end subroutine testCountHasProperty

  subroutine testGetSubsetA(this)
    use GenericType_mod
    use KeyValuePair_mod
    type (fixture) :: this

    type (TracerBundle_type) :: subset

    call makeSubset(this%bundle, withProperty='optionA', subset=subset)
    call assertEqual(1, getCount(subset))
    call assertTrue(associated(getTracer(subset, 'speciesA')))
    call clean(subset)

  end subroutine testGetSubsetA

  subroutine testGetSubsetB(this)
    use GenericType_mod
    use KeyValuePair_mod
    type (fixture) :: this

    type (TracerBundle_type) :: subset

    call makeSubset(this%bundle, withProperty='optionD', subset=subset)
    call assertEqual(2, getCount(subset))
    call clean(subset)

  end subroutine testGetSubsetB

  subroutine testCheckSubsetDefault(this)
    use GenericType_mod
    use KeyValuePair_mod
    type (fixture) :: this

    type (TracerBundle_type) :: subset
    integer :: found

    call addDefault(this%bundle, 'propertyA', -1)
    call makeSubset(this%bundle, withProperty='optionD', subset=subset)

    call addTracer(subset, Tracer('speciesC'))
    call assertFalse(hasProperty(getTracer(subset, 'speciesC'), 'propertyB'))
    call assertTrue(hasProperty(getTracer(subset, 'speciesC'), 'propertyA'))
    found = getProperty(getTracer(subset, 'speciesC'), 'propertyA')
    call assertEqual(-1, found)

    call clean(subset)

  end subroutine testCheckSubsetDefault

  subroutine testAddMandatoryA(this)
    use KeyValuePair_mod
    type (fixture) :: this

#ifdef USE_PFUNIT
    call addMandatoryProperty(this%bundle,'molecularMass')
    call addTracer(this%bundle,Tracer('otherTracer'))
    call assertFailedAssert( &
         & "TracerBundle_mod - species 'otherTracer' is missing mandatory property 'molecularMass'.", &
         & "Failed to detect missing mandatory option.")
#endif

  end subroutine testAddMandatoryA

  Subroutine testAddMandatoryB(this)
    use KeyValuePair_mod
    type (fixture) :: this

#ifdef USE_PFUNIT
    call addMandatoryProperty(this%bundle,'molecularMass')
    call addTracer(this%bundle,Tracer('otherName'))
    call assertFailedAssert( &
         & "TracerBundle_mod - species 'otherName' is missing mandatory property 'molecularMass'.", &
         & "Failed to detect missing mandatory option.")
#endif

  end subroutine testAddMandatoryB

  subroutine testAddMandatoryC(this)
    use KeyValuePair_mod
    type (fixture) :: this

#ifdef USE_PFUNIT
    call addMandatoryProperty(this%bundle, 'newProperty')

    call assertFailedAssert( &
         & "TracerBundle_mod - species 'speciesA' is missing mandatory property 'newProperty'.", &
         & "Failed to detect missing mandatory option.")
#endif
  end subroutine testAddMandatoryC

  subroutine testDefaultValues(this)
    use GenericType_mod
    type (fixture) :: this
    logical :: expectedOptionC, foundOptionC
    integer :: expectedOptionD, foundOptionD
    
    call addDefault(this%bundle, 'optionC', .true.)
    call addDefault(this%bundle, 'optionD',  1)

    foundOptionC = getProperty(this%bundle, 'speciesA', 'optionC')
    expectedOptionC = .false.
    call assertTrue(expectedOptionC .eqv. foundOptionC)

    foundOptionC = getProperty(this%bundle, 'speciesB', 'optionC')
    expectedOptionC = .true.
    call assertTrue(expectedOptionC .eqv. foundOptionC)

    foundOptionD = getProperty(this%bundle, 'speciesA', 'optionD')
    expectedOptionD = 1
    call assertTrue(expectedOptionD == foundOptionD)
    
  end subroutine testDefaultValues

! Check that a new tracers have defaults applied
  subroutine testDefaultAddNewTracer(this)
    use GenericType_mod
    type (fixture) :: this
    
    integer, parameter :: EXPECTED_VALUE = 1
    integer, parameter :: OTHER_VALUE = 2
    character(len=*), parameter :: NAME_1 = 'newTracer_1'
    character(len=*), parameter :: NAME_2 = 'newTracer_2'
    character(len=*), parameter :: PROPERTY = 'newProperty'
    integer :: foundValue

    type (Tracer_type) :: aTracer
    type (TracerBundle_type) :: aBundle

!TODO - this is only needed because of unbalanced setUp/tearDown
    this%bundle = TracerBundle()

    aBundle = TracerBundle()
    call addDefault(aBundle, PROPERTY, EXPECTED_VALUE)

    call addTracer(aBundle, Tracer(NAME_1))

    aTracer = Tracer(NAME_2)
    call setProperty(aTracer, PROPERTY, OTHER_VALUE)
    call addTracer(aBundle, aTracer)

    foundValue = getProperty(aBundle, NAME_1, PROPERTY)
    call assertEqual(EXPECTED_VALUE, foundValue, 'default should be applied')

    foundValue = getProperty(aBundle, NAME_2, PROPERTY)
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
    character(len=*), parameter :: NAME_1 = 'newTracer_1'
    character(len=*), parameter :: NAME_2 = 'newTracer_2'
    character(len=*), parameter :: PROPERTY = 'newProperty'
    integer :: foundValue

    type (Tracer_type) :: aTracer
    type (TracerBundle_type) :: aBundle

!TODO - this is only needed because of unbalanced setUp/tearDown
    this%bundle = TracerBundle()

    aBundle = TracerBundle()
    call addTracer(aBundle, Tracer(NAME_1))

    aTracer = Tracer(NAME_2)
    call setProperty(aTracer, PROPERTY, OTHER_VALUE)
    call addTracer(aBundle, aTracer)

    call addDefault(aBundle, PROPERTY, EXPECTED_VALUE)

    foundValue = getProperty(aBundle, NAME_1, PROPERTY)
    call assertEqual(EXPECTED_VALUE, foundValue, 'should obtain default')

    foundValue = getProperty(aBundle, NAME_2, PROPERTY)
    call assertEqual(OTHER_VALUE, foundValue, 'should not override existing value')

    call clean(aTracer)
    call clean(aBundle)

  end subroutine testDefaultApplyOld

end module testTracers_mod
