module KeyValuePair_mod
  use GenericType_mod
  implicit none
  private

  public :: KeyValuePair
  public :: KeyValuePair_type
  public :: getKey
  public :: getNumValues
  public :: getValue, getValues
  public :: writeUnformatted, readUnformatted
  public :: operator(==)
  public :: clean

  public :: MAX_LEN_KEY
  integer, parameter :: MAX_LEN_KEY = 128

  type KeyValuePair_type
    character(len=MAX_LEN_KEY) :: key
    type (GenericType_type), pointer :: values(:) => null()
  end type KeyValuePair_type

  interface KeyValuePair
    module procedure KeyValuePair_scalar
    module procedure KeyValuePair_array
    module procedure KeyValuePair_copy
  end interface

  interface clean
    module procedure cleanKeyValuePair
  end interface

  interface getValue
    module procedure getValue_1
    module procedure getValue_i
  end interface

  interface operator(==)
    module procedure equals
  end interface

  interface readUnformatted
    module procedure readUnformatted_pair
  end interface

  interface writeUnformatted
    module procedure writeUnformatted_pair
  end interface

Contains

  function KeyValuePair_scalar(key, value) result(pair)
    character(len=*), intent(in) :: key
    type (GenericType_type), intent(in) :: value
    type (KeyValuePair_type) :: pair
    
    pair%key = key
    allocate(pair%values(1))
    pair%values(1) = value
  end function KeyValuePair_scalar

  function KeyValuePair_array(key, values) result(pair)
    character(len=*), intent(in) :: key
    type (GenericType_type), intent(in) :: values(:)
    type (KeyValuePair_type) :: pair
    
    pair%key = key
    allocate(pair%values(size(values)))
    pair%values = values
  end function KeyValuePair_array

  function KeyValuePair_copy(original) result(pair)
    type (KeyValuePair_type), intent(in) :: original
    type (KeyValuePair_type) :: pair
    
    pair%key = original%key
    allocate(pair%values(size(original%values)))
    pair%values = original%values
  end function KeyValuePair_copy

  function getKey(this) result(key)
    type (KeyValuePair_type) :: this
    character(len=MAX_LEN_KEY) :: key
    key = trim(this%key)
  end function getKey

  function getNumValues(this) result(numValues)
    type (KeyValuePair_type), intent(in) :: this
    integer :: numValues
    numValues = size(this%values)
  end function getNumValues

  function getValue_1(this) result(value)
    type (KeyValuePair_type), intent(in) :: this
    type (GenericType_type) :: value
    
    value = this%values(1)
  end function getValue_1

  function getValues(this) result(values)
    type (KeyValuePair_type), intent(in) :: this
    type (GenericType_type) :: values(size(this%values))
    
    values = this%values(:)
  end function getValues

  function getValue_i(this, i) result(value)
    type (KeyValuePair_type), intent(in) :: this
    integer, intent(in) :: i
    type (GenericType_type) :: value
    
    
    value = this%values(i)
  end function getValue_i

  subroutine readUnformatted_pair(this, unit)
    type (KeyValuePair_type), intent(out) :: this
    integer, intent(in) :: unit

    character (len=MAX_LEN_KEY) :: key
    type (GenericType_type), allocatable :: values(:)
    integer :: i, n

    read(unit) key, n
    allocate(values(n))
    do i = 1, n
      call readUnformatted(values(i), unit)
    end do
    this = KeyValuePair(key, values)
  end subroutine readUnformatted_pair

  subroutine writeUnformatted_pair(this, unit)
    type (KeyValuePair_type), intent(in) :: this
    integer, intent(in) :: unit

    integer :: i, n

    n = getNumValues(this)
    write(unit) this%key, n
    do i = 1, n
      call writeUnformatted(this%values(i), unit)
    end do

  end subroutine writeUnformatted_pair

  logical function check(this, valueType, numValues)
    type (KeyValuePair_type), intent(in) :: this
    integer, intent(in) :: valueType
    integer, intent(in) :: numValues

    check = all(getType(this%values) == valueType)
    if (.not. check) then
      call throwException('Incorrect type for specified key: <' &
           & // trim(this%key) // '>', 14)
      return
    end if

    check = (numValues == size(this%values))
    if (.not. check) then
      call throwException('Incorrect number of elements for specified key: <' &
           & // trim(this%key) // '>', 14)
      return
    end if

  end function check

  logical function equals(pairA, pairB) result(isEqual)
    type (KeyValuePair_type), intent(in) :: pairA
    type (KeyValuePair_type), intent(in) :: pairB
    
    isEqual = .true.

    if (trim(getKey(pairA)) /= trim(getKey(pairB))) then
      isEqual = .false.
      return
    end if

    if (getNumValues(pairA) /= getNumValues(pairB)) then
      isEqual = .false.
      return
    end if

    if (.not. all(getValues(pairA)  == getValues(pairB))) then
      isEqual = .false.
      return
    end if

  end function equals

  subroutine cleanKeyValuePair(this)
    type (KeyValuePair_type), intent(in) :: this

    deallocate(this%values)

  end subroutine cleanKeyValuePair

end module KeyValuePair_mod
