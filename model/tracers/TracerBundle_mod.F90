module TracerBundle_mod
  use Dictionary_mod, only: Dictionary_type, Dictionary
  use KeyValuePair_mod, only: MAX_LEN_KEY
  use Tracers_mod
  implicit none
  private

  public :: TracerBundle_type ! derived type
  public :: TracerBundle      ! constructor
  public :: addTracer
  public :: addDefault
  public :: readFromText      
  public :: writeFormatted
  public :: readUnformatted
  public :: writeUnformatted
  public :: getTracer
  public :: getProperty
  public :: getProperties
  public :: setProperty
  public :: getIndex
  public :: getCount
  public :: makeSubset
  public :: hasProperty
  public :: addMandatoryProperty
  public :: operator(==)
  public :: clean

  public :: NOT_FOUND

  type TracerBundle_type
    private
    type (Dictionary_type) :: defaultValues
!TODO: Should change to allocatable when compilers work
    character(len=MAX_LEN_KEY), pointer :: mandatoryProperties(:) => null()
    type (Tracer_type), pointer :: tracers(:) => null()
  end type TracerBundle_type

  interface clean
    module procedure cleanBundle
  end interface

  integer, parameter :: NOT_FOUND = -1

  interface getProperty
    module procedure getProperty_multi
  end interface

  interface getProperties
    module procedure getProperties_multi
  end interface

  interface hasProperty
    module procedure hasProperty_multi
  end interface

  interface setProperty
    module procedure setProperty_multi_integer, setProperty_multi_integerArr
    module procedure setProperty_multi_real64,  setProperty_multi_real64Arr
    module procedure setProperty_multi_logical, setProperty_multi_logicalArr
    module procedure setProperty_multi_string,  setProperty_multi_stringArr
  end interface

  interface getCount
    module procedure getCount_
    module procedure getCount_property
  end interface

  interface operator(==)
    module procedure equals
  end interface

  interface writeFormatted
    module procedure writeFormatted_bundle
  end interface

  interface writeUnformatted
    module procedure writeUnformatted_bundle
  end interface

  interface readUnformatted
    module procedure readUnformatted_bundle
  end interface

  interface addDefault
    module procedure addDefault_integer
    module procedure addDefault_logical
    module procedure addDefault_real64
    module procedure addDefault_string
  end interface

  integer, parameter :: LEN_HEADER = 80
  integer, parameter :: VERSION = 1
  character(len=*), parameter :: DESCRIPTION = 'TracerBundle'

contains

  function TracerBundle()
!@sum Construct empty tracer    
    use Dictionary_mod, only: Dictionary
    type (TracerBundle_type) :: TracerBundle

    allocate(TracerBundle%tracers(0))
    allocate(TracerBundle%mandatoryProperties(0))
    TracerBundle%defaultValues = Dictionary()

  end function TracerBundle

  function readFromText(unit, defaultValues) result(bundle)
!@sum Populate a TracerBundle from a unit attached to a formatted
!@+ file.  Optionally apply default values.
    use Parser_mod
    use Dictionary_mod
    integer, intent(in) :: unit
    type (Dictionary_type), optional, intent(in) :: defaultValues
    type (TracerBundle_type) :: bundle
    type (Tracer_type) :: tracer
    type (Parser_type) :: parser

    integer :: status

    bundle = TracerBundle()
    if (present(defaultValues)) bundle%defaultValues = Dictionary(defaultValues)

    do
      tracer = readOneTracer(unit, status)
      if (status /= 0) exit

      call addTracer(bundle, tracer)
      call clean(tracer)
    end do

  end function readFromText

  subroutine addTracer(this, aTracer)
!@sum Insert one tracer into a bundle.  Optionally
!@+ apply default values to the tracer dictionary.
    use Dictionary_mod
    type (TracerBundle_type), intent(inout) :: this
    type (Tracer_type), intent(in) :: aTracer

    type (Tracer_type), allocatable :: oldList(:)
    integer :: count


    call assertHasProperties(aTracer, this%mandatoryProperties) 

    count = getNumTracers(this)

    allocate(oldList(count))
    oldList = this%tracers ! shallow copy ok
    deallocate(this%tracers)
    allocate(this%tracers(count+1))
    this%tracers(1:count) = oldList
    deallocate(oldList)

    this%tracers(count+1) = Tracer(aTracer)
    call merge(this%tracers(count+1), this%defaultValues)

  end subroutine addTracer

  subroutine mergeDefaults(this, pair)
    use Dictionary_mod, only: merge
    use KeyValuePair_mod, only: KeyValuePair_type
    type (TracerBundle_type), intent(inout) :: this
    type (KeyValuePair_type), intent(in) :: pair

    integer :: i

    do i = 1, getNumTracers(this)
      call merge(this%tracers(i), pair)
    end do

  end subroutine mergeDefaults

  subroutine addDefault_integer(this, property, value)
!@sum Add a default property that applies to all tracers in bundle,
!@+ unless overridden with specific value for a given tracer.
    use Dictionary_mod, only: insert, merge
    use KeyValuePair_mod, only: KeyValuePair
    use GenericType_mod, only: GenericType
    type (TracerBundle_type), intent(inout) :: this
    character(len=*), intent(in) :: property
    integer, intent(in) :: value

    call insert(this%defaultValues, property, value)
    call mergeDefaults(this, KeyValuePair(property, GenericType(value)))

  end subroutine addDefault_integer

  subroutine addDefault_logical(this, property, value)
!@sum Add a default property that applies to all tracers in bundle,
!@+ unless overridden with specific value for a given tracer.
    use Dictionary_mod, only: insert
    use KeyValuePair_mod, only: KeyValuePair
    use GenericType_mod, only: GenericType
    type (TracerBundle_type), intent(inout) :: this
    character(len=*), intent(in) :: property
    logical, intent(in) :: value

    call insert(this%defaultValues, property, value)
    call mergeDefaults(this, KeyValuePair(property, GenericType(value)))

  end subroutine addDefault_logical

  subroutine addDefault_real64(this, property, value)
!@sum Add a default property that applies to all tracers in bundle,
!@+ unless overridden with specific value for a given tracer.
    use Dictionary_mod, only: insert
    use KeyValuePair_mod, only: KeyValuePair
    use GenericType_mod, only: GenericType
    type (TracerBundle_type), intent(inout) :: this
    character(len=*), intent(in) :: property
    real*8, intent(in) :: value

    call insert(this%defaultValues, property, value)
    call mergeDefaults(this, KeyValuePair(property, GenericType(value)))

  end subroutine addDefault_real64

  subroutine addDefault_string(this, property, value)
!@sum Add a default property that applies to all tracers in bundle,
!@+ unless overridden with specific value for a given tracer.
    use Dictionary_mod, only: insert
    use KeyValuePair_mod, only: KeyValuePair
    use GenericType_mod, only: GenericType
    type (TracerBundle_type), intent(inout) :: this
    character(len=*), intent(in) :: property
    character(len=*), intent(in) :: value

    call insert(this%defaultValues, property, value)
    call mergeDefaults(this, KeyValuePair(property, GenericType(value)))

  end subroutine addDefault_string

  integer function getNumTracers(this)
!@sum Returns number of tracers in a bundle.
    type (TracerBundle_type), intent(in) :: this
    getNumTracers = size(this%tracers)
  end function getNumTracers

  subroutine writeUnformatted_bundle(this, unit)
!@sum Write a bundle to a unit attached to an unformatted sequential file.
    type (TracerBundle_type), intent(in) :: this
    integer, intent(in) :: unit

    integer :: i, n

    character(len=LEN_HEADER) :: header

    write(header, '(a," - version:",1x,i5.0)') DESCRIPTION, VERSION
    write(unit) header

    n = getNumTracers(this)
    write(unit) n
    do i = 1, n
      call writeUnformatted(this%tracers(i), unit)
    end do

  end subroutine writeUnformatted_bundle

  subroutine readUnformatted_bundle(this, unit)
!@sum Read a bundle to a unit attached to an unformatted sequential file.
    type (TracerBundle_type), intent(out) :: this
    integer, intent(in) :: unit

    integer :: i, n
    integer :: oldVersion
    character(len=len(DESCRIPTION)) :: tag
    character(len=LEN_HEADER) :: header

    read(unit) header
    read(header, '(a,11x,i10.0)') tag, oldVersion

    if (tag /= DESCRIPTION) then
      call throwException(DESCRIPTION // '::readUnformatted() - incorrect header.', 14)
    end if

    if (oldVersion /= VERSION) then
      call throwException(DESCRIPTION // '::readUnformatted() - unsupported format.', 14)
    end if
      
    read(unit) n
    this = TracerBundle()
    allocate(this%tracers(n))
    do i = 1, n
      call readUnformatted(this%tracers(i), unit)
    end do

  end subroutine readUnformatted_bundle

  subroutine writeFormatted_bundle(this, unit)
    use Parser_mod
    type (TracerBundle_type), intent(in) :: this
    integer, intent(in) :: unit

    type (Parser_type) :: parser
    type (Dictionary_type), pointer :: properties
    integer :: i

    call setBeginData(parser, '{')
    call setEndData(parser, '}')
    call setTokenSeparators(parser, '=,')

    do i = 1, getCount(this)
      properties => getProperties(this%tracers(i))
      call writeFormatted(parser, unit, properties)
    end do
    
  end subroutine writeFormatted_bundle

  integer function getIndex(this, name) result(index)
    use Dictionary_mod, only: lookup
    use KeyValuePair_mod, only: MAX_LEN_KEY
    use GenericType_mod
    type (TracerBundle_type), intent(in) :: this
    character(len=*), intent(in) :: name

    integer :: i, j
    character(len=MAX_LEN_KEY) :: aName

    index = NOT_FOUND

    j = getCount(this)

    do i = 1, j
      if (all([name] == lookup(getProperties(this%tracers(i)), 'name'))) then
        index = i
        exit
      end if
    end do

  end function getIndex

  function hasProperty_multi(this, property) result(has)
    use Dictionary_mod, only: hasKey
    type (TracerBundle_type), intent(in) :: this
    character(len=*), intent(in) :: property
    logical :: has(size(this%tracers))

    integer :: i

    has = (/ (hasProperty(this%tracers(i), property), i = 1, getCount(this)) /)
  end function hasProperty_multi

  subroutine addMandatoryProperty(this, property)
    type (TracerBundle_type), intent(inout) :: this
    character(len=*), intent(in) :: property

    character(len=MAX_LEN_KEY), allocatable :: oldProperties(:)
    integer :: n
    integer :: i

    n = size(this%mandatoryProperties)
    
    allocate(oldProperties(n))
    oldProperties = this%mandatoryProperties
    deallocate(this%mandatoryProperties)
    allocate(this%mandatoryProperties(n+1))
    this%mandatoryProperties(1:n) = oldProperties
    deallocate(oldProperties)

    this%mandatoryProperties(n+1) = property

    do i = 1, getNumTracers(this)
      if (.not. assertHasProperty(this%tracers(i), property)) return
    end do

  end subroutine addMandatoryProperty

  subroutine assertHasProperties(this, properties)
    type (Tracer_type), intent(in) :: this
    character(len=MAX_LEN_KEY), intent(in) :: properties(:)
    character(len=MAX_LEN_KEY) :: name

    integer :: i

    do i = 1, size(properties)
      if (.not. assertHasProperty(this, properties(i))) return
    end do

  end subroutine assertHasProperties

  logical function assertHasProperty(this, property)
    type (Tracer_type), intent(in) :: this
    character(len=*), intent(in) :: property

    character(len=MAX_LEN_KEY) :: name

    assertHasProperty = .true.
    if (.not. hasProperty(this, property)) then
      name = getName(this)
      call throwException("TracerBundle_mod - species '" // trim(name) // &
           & "' is missing mandatory property '" // trim(property) // "'.", 14)
      assertHasProperty = .false.
    end if

  end function assertHasProperty

  function getName(this) result(name)
    use GenericType_mod
    use KeyValuePair_mod, only: MAX_LEN_KEY
    type (Tracer_type), intent(in) :: this
    character(len=MAX_LEN_KEY) :: name
    name = getProperty(this, 'name')
  end function getName

  function getTracer(this, species) result(tracer)
    type (TracerBundle_type), target, intent(in) :: this
    character(len=*), intent(in) :: species
    type (Tracer_type), pointer :: tracer


    integer :: index

    index = getIndex(this, species)
    if (index /= NOT_FOUND) then
      tracer => this%tracers(index)
    else
      tracer => null()
    end if
    
  end function getTracer

  function getProperties_multi(this, species) result(properties)
    type (TracerBundle_type), target, intent(in) :: this
    character(len=*), intent(in) :: species
    type (Dictionary_type), pointer :: properties

    integer :: index

    index = getIndex(this, species)
    if (index /= NOT_FOUND) then
      properties => getProperties(this%tracers(index))
    else
      allocate(properties)
      properties = Dictionary()
    end if

  end function getProperties_multi

  function getProperty_multi(this, species, property) result (propertyValues)
    use GenericType_mod
    use Dictionary_mod, only: lookup
    type (TracerBundle_type), intent(in) :: this
    character(len=*), intent(in) :: species
    character(len=*), intent(in) :: property
    type (GenericType_type), pointer :: propertyValues(:)

    integer :: index

    index = getIndex(this, species)
    if (index /= NOT_FOUND) then
      propertyValues => getProperty(this%tracers(index), property)
    else
      allocate(propertyValues(0))
    end if
    
  end function getProperty_multi


  ! Select a tracer from list
  subroutine setProperty_multi_integer(this, name, property, value)
    use Dictionary_mod, only: insert
    type (TracerBundle_type), intent(inout) :: this
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: property
    integer, intent(in) :: value
    
    integer :: index

    index = getIndex(this, name)
    call setProperty(this%tracers(index), property, value)

  end subroutine setProperty_multi_integer

  subroutine setProperty_multi_integerArr(this, name, property, values)
    use Dictionary_mod, only: insert
    type (TracerBundle_type), intent(inout) :: this
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: property
    integer, intent(in) :: values(:)
    
    integer :: index

    index = getIndex(this, name)
    call setProperty(this%tracers(index), property, values)

  end subroutine setProperty_multi_integerArr

  subroutine setProperty_multi_real64(this, name, property, value)
    use Dictionary_mod, only: insert
    type (TracerBundle_type), intent(inout) :: this
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: property
    real*8, intent(in) :: value
    
    integer :: index

    index = getIndex(this, name)
    call setProperty(this%tracers(index), property, value)

  end subroutine setProperty_multi_real64

  subroutine setProperty_multi_real64Arr(this, name, property, values)
    use Dictionary_mod, only: insert
    type (TracerBundle_type), intent(inout) :: this
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: property
    real*8, intent(in) :: values(:)
    
    integer :: index

    index = getIndex(this, name)
    call setProperty(this%tracers(index), property, values)
  end subroutine setProperty_multi_real64Arr

  subroutine setProperty_multi_logical(this, name, property, value)
    use Dictionary_mod, only: insert
    type (TracerBundle_type), intent(inout) :: this
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: property
    logical, intent(in) :: value
    
    integer :: index

    index = getIndex(this, name)
    call setProperty(this%tracers(index), property, value)
  end subroutine setProperty_multi_logical

  subroutine setProperty_multi_logicalArr(this, name, property, values)
    use Dictionary_mod, only: insert
    type (TracerBundle_type), intent(inout) :: this
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: property
    logical, intent(in) :: values(:)
    
    integer :: index

    index = getIndex(this, name)
    call setProperty(this%tracers(index), property, values)
  end subroutine setProperty_multi_logicalArr

  subroutine setProperty_multi_string(this, name, property, value)
    use Dictionary_mod, only: insert
    type (TracerBundle_type), intent(inout) :: this
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: property
    character(len=*), intent(in) :: value
    
    integer :: index

    index = getIndex(this, name)
    call setProperty(this%tracers(index), property, value)
  end subroutine setProperty_multi_string

  subroutine setProperty_multi_stringArr(this, name, property, values)
    use Dictionary_mod, only: insert
    type (TracerBundle_type), intent(inout) :: this
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: property
    character(len=*), intent(in) :: values(:)
    
    integer :: index

    index = getIndex(this, name)
    call setProperty(this%tracers(index), property, values)
  end subroutine setProperty_multi_stringArr

  integer function getCount_(this) result(number)
    type (TracerBundle_type), intent(in) :: this

    number = size(this%tracers)
    
  end function getCount_

  integer function getCount_property(this, withProperty) result(number)
    type (TracerBundle_type), intent(in) :: this
    character(len=*), intent(in) :: withProperty

    number = count(hasProperty(this, withProperty))
    
  end function getCount_property

  subroutine makeSubset(this, withProperty, subset)
    use Dictionary_mod, only: Dictionary
    type (TracerBundle_type), target, intent(in) :: this
    character(len=*), intent(in) :: withProperty
    type (TracerBundle_type) :: subset

    integer :: i

    subset = TracerBundle()
    subset%defaultValues = Dictionary(this%defaultValues)

    do i = 1, getCount(this)
      if (hasProperty(this%tracers(i), withProperty)) then
        call addTracer(subset, this%tracers(i))
      end if
    end do

  end subroutine makeSubset

  logical function equals(bundleA, bundleB) result(isEqual)
    use GenericType_mod
    use Dictionary_mod, only: operator(==)
    use Parser_mod, only: MAX_LEN_TOKEN
    type (TracerBundle_type), target, intent(in) :: bundleA
    type (TracerBundle_type), target, intent(in) :: bundleB

    integer :: i
    character(len=MAX_LEN_TOKEN) :: name

    isEqual = .true.
    do i = 1, getCount(bundleA)
      name = getProperty(bundleA%tracers(1), 'name')
      if (.not. (getProperties(bundleA, name) == getProperties(bundleB, name))) then
        isEqual = .false.
        exit
      end if
    end do

  end function equals

  subroutine cleanBundle(this)
    use Dictionary_mod, only: clean
    type (TracerBundle_type), intent(inout) :: this
    integer :: i
    do i = 1, getCount(this)
      call clean(this%tracers(i))
    end do
    deallocate(this%tracers)
    call clean(this%defaultValues)
!    if (size(this%mandatoryProperties)>0) then
!       print *, 'SIZE = ',size(this%mandatoryProperties)
       deallocate(this%mandatoryProperties)
!    end if
  end subroutine cleanBundle

end module TracerBundle_mod
