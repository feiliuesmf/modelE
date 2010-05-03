module GenericType_mod
  implicit none
  private

  public :: GenericType
  public :: GenericType_type
  public :: operator(==)
  public :: assignment(=)
  public :: toString
  public :: fromString
  public :: getType

  public :: MAX_LEN_VALUE
  public :: INTEGER_TYPE
  public :: REAL64_TYPE
  public :: LOGICAL_TYPE
  public :: STRING_TYPE

  integer, parameter :: MAX_LEN_VALUE = 128
  integer, parameter :: MAX_LEN_TYPE_STRING = 10

  type GenericType_type
    integer :: type = 0
    integer :: integerValue
    real*8 :: real64Value
    logical :: logicalvalue
    character(len=MAX_LEN_VALUE) :: stringValue
  end type GenericType_type

  interface GenericType ! constructors
    module procedure GenericType_copy
    module procedure GenericType_integer
    module procedure GenericType_real64
    module procedure GenericType_logical
    module procedure GenericType_string
  end interface

  interface operator(==)
    module procedure equals_generic
    module procedure equals_integer
!!$    module procedure equals_integerArray
    module procedure equals_real64
    module procedure equals_logical
    module procedure equals_string
  end interface

  interface assignment(=)
    module procedure assignToInteger
    module procedure assignToIntegerArray
    module procedure assignToReal64
    module procedure assignToReal64Array
    module procedure assignToLogical
    module procedure assignToLogicalArray
    module procedure assignToString
    module procedure assignToStringArray
  end interface

  interface getType
    module procedure getType_generic
    module procedure getType_integer
    module procedure getType_real64
    module procedure getType_logical
    module procedure getType_string
  end interface

  integer, parameter :: INTEGER_TYPE = 1
  integer, parameter :: REAL64_TYPE  = 2
  integer, parameter :: LOGICAL_TYPE = 3
  integer, parameter :: STRING_TYPE  = 4

contains

  function GenericType_copy(original) result(generic)
    type (GenericType_type), intent(in) :: original
    type (GenericType_type) :: generic

    generic = original

  end function GenericType_copy

  elemental function GenericType_integer(value) result(generic)
    integer, intent(in) :: value
    type (GenericType_type) :: generic

    generic%type = INTEGER_TYPE
    generic%integerValue = value
  end function GenericType_integer

  elemental function GenericType_real64(value) result(generic)
    real*8, intent(in) :: value
    type (GenericType_type) :: generic

    generic%type = REAL64_TYPE
    generic%real64Value = value
  end function GenericType_real64

  elemental function GenericType_logical(value) result(generic)
    logical, intent(in) :: value
    type (GenericType_type) :: generic

    generic%type = LOGICAL_TYPE
    generic%logicalValue = value
  end function GenericType_logical

  elemental function GenericType_string(value) result(generic)
    character(len=*), intent(in) :: value
    type (GenericType_type) :: generic

    generic%type = STRING_TYPE
    generic%stringValue = value
  end function GenericType_string

  elemental logical function equals_generic(expected, generic) result(same)
    type (GenericType_type), intent(in) :: expected
    type (GenericType_type), intent(in) :: generic
    same = (expected%type == generic%type)
    if (same) then
      select case (generic%type)
      case (INTEGER_TYPE)
        same = (expected%integerValue == generic%integerValue)
      case (REAL64_TYPE)
        same = (expected%real64Value == generic%real64Value)
      case (LOGICAL_TYPE)
        same = (expected%logicalValue == generic%logicalValue)
      case (STRING_TYPE)
        same = (expected%stringValue == generic%stringValue)
      end select
    end if
  end function equals_generic
  
  elemental logical function equals_integer(expected, generic) result(same)
    integer, intent(in) :: expected
    type (GenericType_type), intent(in) :: generic
    same = (generic%type==INTEGER_TYPE) .and. (generic%integerValue==expected)
  end function equals_integer

  elemental logical function equals_real64(expected, generic) result(same)
    real*8, intent(in) :: expected
    type (GenericType_type), intent(in) :: generic
    same = (generic%type == REAL64_TYPE) .and. (generic%real64Value==expected)
  end function equals_real64

  elemental logical function equals_logical(expected, generic) result(same)
    logical, intent(in) :: expected
    type (GenericType_type), intent(in) :: generic
    same = (generic%type==LOGICAL_TYPE) .and. (generic%logicalValue==expected)
  end function equals_logical

  elemental logical function equals_string(expected, generic) result(same)
    character(len=*), intent(in) :: expected
    type (GenericType_type), intent(in) :: generic
    same = (generic%type==STRING_TYPE) .and. (generic%stringValue == expected)
  end function equals_string

  subroutine assignToInteger(value, this)
    integer, intent(out) :: value
    type (GenericType_type), intent(in) :: this
    value = this%integerValue
  end subroutine assignToInteger

  subroutine assignToIntegerArray(values, this)
    integer, intent(out) :: values(:)
    type (GenericType_type), intent(in) :: this(:)
    values(:) = this(:)%integerValue
  end subroutine assignToIntegerArray
  
  subroutine assignToReal64(value, this)
    real*8, intent(out) :: value
    type (GenericType_type), intent(in) :: this
    value = this%real64Value
  end subroutine assignToReal64

  subroutine assignToReal64Array(values, this)
    real*8, intent(out) :: values(:)
    type (GenericType_type), intent(in) :: this(:)
    values = this(:)%real64Value
  end subroutine assignToReal64Array

  subroutine assignToLogical(value, this)
    logical, intent(out) :: value
    type (GenericType_type), intent(in) :: this
    value = this%logicalValue
  end subroutine assignToLogical

  subroutine assignToLogicalArray(values, this)
    logical, intent(out) :: values(:)
    type (GenericType_type), intent(in) :: this(:)
    values = this(:)%logicalValue
  end subroutine assignToLogicalArray

  subroutine assignToString(value, this)
    character(len=*), intent(out) :: value
    type (GenericType_type), intent(in) :: this
    value = this%stringValue
  end subroutine assignToString

  subroutine assignToStringArray(values, this)
    character(len=*), intent(out) :: values(:)
    type (GenericType_type), intent(in) :: this(:)
    values = this%stringValue
  end subroutine assignToStringArray

  function toString(this) result(string)
    type (GenericType_type), intent(in) :: this
    character(len=MAX_LEN_VALUE) :: string
    select case (this%type)
    case (INTEGER_TYPE)
      write(string, '(i0)') this%integerValue
    case (REAL64_TYPE)
      write(string, '(g25.17)') this%real64Value
    end select
  end function toString
  
  function fromString(string, type) result(generic)
    use StringUtilities_mod, only: toLowerCase
    character(len=*), intent(in) :: string
    integer, intent(in) :: type
    type (GenericType_type) :: generic
    integer :: status

    generic%type = type
    select case (type)
    case (INTEGER_TYPE)
      generic%type = INTEGER_TYPE
      read(string,'(i)',iostat=status) generic%integerValue
    case (REAL64_TYPE)
      generic%type = REAL64_TYPE
      read(string,'(g)',iostat=status) generic%real64Value
    case (LOGICAL_TYPE)
      generic%type = LOGICAL_TYPE
      status = 0
      select case(trim(toLowerCase(string)))
      case ('true','t','.true.')
        generic%logicalValue = .true.
      case ('false','f','.false.')
        generic%logicalValue = .false.
      case default
        status = -1 ! cannot convert
      end select
    case (STRING_TYPE)
      generic%type = STRING_TYPE
      generic%stringValue = string
    case default
      call throwException('GenericType::fromString() - no such type.',14)
    end select

    if (status /= 0) then
      call throwException('GenericType::fromString() - cannot convert string "' // &
           & trim(string) // '" to ' // trim(typeString(type)) // '.', 14)
    end if
  contains

    function typeString(type)
      integer, intent(in) :: type
      character(len=MAX_LEN_TYPE_STRING) :: typeString

      select case (type)
      case (INTEGER_TYPE)
        typeString = 'integer'
      case (REAL64_TYPE)
        typeString = 'real*8'
      case (LOGICAL_TYPE)
        typeString = 'logical'
      case (STRING_TYPE)
        typeString = 'string'
      end select
    end function typeString

  end function fromString

  elemental integer function getType_generic(this) result(type)
    type (GenericType_type), intent(in) :: this
    type = this%type
  end function getType_generic

  integer function getType_integer(value)
    integer, intent(in) :: value
    getType_integer = INTEGER_TYPE
  end function getType_integer

  integer function getType_real64(value)
    real*8, intent(in) :: value
    getType_real64 = REAL64_TYPE
  end function getType_real64

  integer function getType_logical(value)
    logical, intent(in) :: value
    getType_logical = LOGICAL_TYPE
  end function getType_logical

  integer function getType_string(value)
    character(len=*), intent(in) :: value
    getType_string = STRING_TYPE
  end function getType_string

end module GenericType_mod
