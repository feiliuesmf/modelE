module Tracers_mod
  use Dictionary_mod, only: Dictionary_type
  implicit none
  private

  public :: Tracer_type       ! derived type
  public :: Tracer            ! constructor
  public :: clean

  public :: setProperty
  public :: getProperty
  public :: getProperties
  public :: hasProperty
  public :: merge
  public :: writeUnformatted
  public :: readUnformatted
  public :: readOneTracer

  type Tracer_type
    private
    type (Dictionary_type) :: properties
  end type Tracer_type

  interface Tracer
    module procedure newTracer
    module procedure newTracerName
    module procedure TracerCopy
  end interface

  interface setProperty
    module procedure setProperty_integer, setProperty_integerArr
    module procedure setProperty_real64,  setProperty_real64Arr
    module procedure setProperty_logical, setProperty_logicalArr
    module procedure setProperty_string,  setProperty_stringArr
  end interface

  interface getProperty
    module procedure getProperty_single
  end interface

  interface getProperties
    module procedure getProperties_single
  end interface

  interface hasProperty
    module procedure hasProperty_single
  end interface

  interface merge
    module procedure merge_one
    module procedure merge_multi
  end interface

  interface writeUnformatted
    module procedure writeUnformatted_tracer
  end interface

  interface readUnformatted
    module procedure readUnformatted_tracer
  end interface

  interface clean
    module procedure cleanTracer
  end interface

contains

  function newTracer() result(aTracer)
!@sum Construct empty tracer    
    use Dictionary_mod, only: Dictionary
    type (Tracer_type) :: aTracer
    aTracer%properties = Dictionary()
  end function newTracer

  function newTracerName(name) result(aTracer)
!@sum Construct named tracer
    use Dictionary_mod, only: Dictionary
    character(len=*), intent(in) :: name
    type (Tracer_type) :: aTracer

    aTracer = Tracer()
    call setProperty(aTracer, 'name', name)
  end function newTracerName

  ! Copy constructor
  function TracerCopy(original) result(copy)
    use Dictionary_mod, only: Dictionary
    type (Tracer_type), intent(in) :: original
    type (Tracer_type) :: copy

    copy%properties = Dictionary(original%properties)

  end function TracerCopy

  subroutine setProperty_integer(this, property, value)
    use Dictionary_mod, only: insert
    type (Tracer_type), intent(inout) :: this
    character(len=*), intent(in) :: property
    integer, intent(in) :: value
    call insert(this%properties, property, value)
  end subroutine setProperty_integer

  subroutine setProperty_integerArr(this, property, values)
    use Dictionary_mod, only: insert
    type (Tracer_type), intent(inout) :: this
    character(len=*), intent(in) :: property
    integer, intent(in) :: values(:)
    call insert(this%properties, property, values)
  end subroutine setProperty_integerArr

  subroutine setProperty_real64(this, property, value)
    use Dictionary_mod, only: insert
    type (Tracer_type), intent(inout) :: this
    character(len=*), intent(in) :: property
    real*8, intent(in) :: value
    call insert(this%properties, property, value)
  end subroutine setProperty_real64

  subroutine setProperty_real64Arr(this, property, values)
    use Dictionary_mod, only: insert
    type (Tracer_type), intent(inout) :: this
    character(len=*), intent(in) :: property
    real*8, intent(in) :: values(:)
    call insert(this%properties, property, values)
  end subroutine setProperty_real64Arr

  subroutine setProperty_logical(this, property, value)
    use Dictionary_mod, only: insert
    type (Tracer_type), intent(inout) :: this
    character(len=*), intent(in) :: property
    logical, intent(in) :: value
    call insert(this%properties, property, value)
  end subroutine setProperty_logical

  subroutine setProperty_logicalArr(this, property, values)
    use Dictionary_mod, only: insert
    type (Tracer_type), intent(inout) :: this
    character(len=*), intent(in) :: property
    logical, intent(in) :: values(:)
    call insert(this%properties, property, values)
  end subroutine setProperty_logicalArr

  subroutine setProperty_string(this, property, value)
    use Dictionary_mod, only: insert
    type (Tracer_type), intent(inout) :: this
    character(len=*), intent(in) :: property
    character(len=*), intent(in) :: value
    call insert(this%properties, property, value)
  end subroutine setProperty_string

  subroutine setProperty_stringArr(this, property, values)
    use Dictionary_mod, only: insert
    type (Tracer_type), intent(inout) :: this
    character(len=*), intent(in) :: property
    character(len=*), intent(in) :: values(:)
    call insert(this%properties, property, values)
  end subroutine setProperty_stringArr

  subroutine merge_one(this, pair)
    use KeyValuePair_mod, only: KeyValuePair_type
    use Dictionary_mod, only: merge
    type (Tracer_type), intent(inout) :: this
    type (KeyValuePair_type), intent(in) :: pair
    
    call merge(this%properties, pair)
  end subroutine merge_one

  function getProperty_single(this, property) result (propertyValues)
    use GenericType_mod
    use Dictionary_mod, only: lookup
    type (Tracer_type), intent(in) :: this
    character(len=*), intent(in) :: property
    type (GenericType_type), pointer :: propertyValues(:)

    propertyValues => lookup(this%properties, property)
    
  end function getProperty_single

  function getProperties_single(this) result(properties)
    type (Tracer_type), target, intent(in) :: this
    type (Dictionary_type), pointer :: properties
    properties => this%properties
  end function getProperties_single

  subroutine merge_multi(this, properties)
    use Dictionary_mod, only: merge
    type (Tracer_type), intent(inout) :: this
    type (Dictionary_type), intent(in) :: properties
    
    call merge(this%properties, properties)
  end subroutine merge_multi

  subroutine writeUnformatted_tracer(this, unit)
!@sum Write a tracer to a unit attached to an unformatted sequential file.
    use Dictionary_mod, only: writeUnformatted
    type (Tracer_type), intent(in) :: this
    integer, intent(in) :: unit

    call writeUnformatted(this%properties, unit)
  end subroutine writeUnformatted_tracer

  subroutine readUnformatted_tracer(this, unit)
!@sum Read a bundle to a unit attached to an unformatted sequential file.
    use Dictionary_mod, only: readUnformatted
    type (Tracer_type), intent(out) :: this
    integer, intent(in) :: unit
    call readUnformatted(this%properties, unit)
  end subroutine readUnformatted_tracer

  function readOneTracer(unit, status) result(aTracer)
    use Parser_mod, only: Parser_type
    use Parser_mod, only: setBeginData, setEndData
    use Parser_mod, only: setCommentCharacters, setTokenSeparators
    use Parser_mod, only: parse

    integer, intent(in) :: unit
    integer, intent(out) :: status
    type (Tracer_type) :: aTracer

    type (Parser_type) :: parser

    call setBeginData(parser, '{')
    call setEndData(parser, '}')
    call setTokenSeparators(parser, '=,')
    call setCommentCharacters(parser, '!#')

    aTracer = Tracer()
    aTracer%properties = parse(parser, unit, status)
    if (status /= 0) return

  end function readOneTracer

  logical function hasProperty_single(this, property) result(has)
    use Dictionary_mod, only: hasKey
    type (Tracer_type), intent(in) :: this
    character(len=*), intent(in) :: property
    has = hasKey(this%properties, property)
  end function hasProperty_single

  subroutine cleanTracer(this)
    use Dictionary_mod, only: clean
    type (Tracer_type), intent(inout) :: this
    call clean(this%properties)
  end subroutine cleanTracer

end module Tracers_mod
