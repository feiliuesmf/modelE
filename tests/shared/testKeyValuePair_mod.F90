module testKeyValuePair_mod
  use pFUnit
  use Dictionary_mod
  implicit none
  private

  public :: testGetKey
  public :: testGetType
  public :: testGetNumValues

  public :: testGetInteger
  public :: testGetIntegerB
  public :: testGetReal64
  public :: testGetLogical
  public :: testGetString
  public :: testGetDictionary

  public :: testGetIntegerArray
  public :: testGetReal64Array
  public :: testGetLogicalArray
  public :: testGetStringArray

#ifdef USE_PFUNIT
  public :: testWrongType
#endif
  
contains

  subroutine testGetKey()
    type (KeyValuePair_type) :: pair

    pair = KeyValuePair('aKey', 1)
    call assertEqual('aKey', getKey(pair)) 

  end subroutine testGetKey

  subroutine testGetType()
    type (KeyValuePair_type) :: pair

    pair = KeyValuePair('aKey', 1)
    call assertEqual(INTEGER_TYPE, getType(pair))

    pair = KeyValuePair('aKey', 1.234d+9)
    call assertEqual(REAL64_TYPE, getType(pair))

    pair = KeyValuePair('aKey', .true.)
    call assertEqual(LOGICAL_TYPE, getType(pair))

    pair = KeyValuePair('aKey', 'hello')
    call assertEqual(STRING_TYPE, getType(pair))

    pair = KeyValuePair('aKey', Dictionary())
    call assertEqual(DICT_TYPE, getType(pair))

    pair = KeyValuePair('aKey', [1])
    call assertEqual(INTEGER_TYPE, getType(pair))

    pair = KeyValuePair('aKey', [1.234d+9])
    call assertEqual(REAL64_TYPE, getType(pair))

    pair = KeyValuePair('aKey', [.true.])
    call assertEqual(LOGICAL_TYPE, getType(pair))

    pair = KeyValuePair('aKey', ['hello'])
    call assertEqual(STRING_TYPE, getType(pair))

    pair = KeyValuePair('aKey', [Dictionary()])
    call assertEqual(DICT_TYPE, getType(pair))

  end subroutine testGetType

  subroutine testGetNumValues()
    type (KeyValuePair_type) :: pair

    pair = KeyValuePair('aKey', 1)
    call assertEqual(1, getNumValues(pair))

    pair = KeyValuePair('aKey', 1.234d+0)
    call assertEqual(1, getNumValues(pair))

    pair = KeyValuePair('aKey', .true.)
    call assertEqual(1, getNumValues(pair))

    pair = KeyValuePair('aKey', 'hello')
    call assertEqual(1, getNumValues(pair))

    ! lists
    pair = KeyValuePair('aKey', [1,2,3])
    call assertEqual(3, getNumValues(pair))

    pair = KeyValuePair('aKey', [1.234d+0])
    call assertEqual(1, getNumValues(pair))

    pair = KeyValuePair('aKey', [.true.,.false.,.true.])
    call assertEqual(3, getNumValues(pair))

    pair = KeyValuePair('aKey', ['hello', 'goodbye'])
    call assertEqual(2, getNumValues(pair))


  end subroutine testGetNumValues

  subroutine testGetInteger()
    type (KeyValuePair_type) :: pair
    integer :: expectedValue, foundValue

    expectedValue = 1
    pair = KeyValuePair('aKey', expectedValue)
    call getValue(pair, foundValue)
    call assertEqual(expectedValue, foundValue)

  end subroutine testGetInteger

  subroutine testGetIntegerB()
    type (KeyValuePair_type) :: pair
    integer :: expectedValue, foundValue

    expectedValue = 2
    pair = KeyValuePair('aKey', expectedValue)
    call getValue(pair, foundValue)
    call assertEqual(expectedValue, foundValue)

  end subroutine testGetIntegerB

  subroutine testGetReal64()
    type (KeyValuePair_type) :: pair
    real(kind=r64) :: expectedValue, foundValue

    expectedValue = 1.234
    pair = KeyValuePair('aKey', expectedValue)
    call getValue(pair, foundValue)
    call assertEqual(expectedValue, foundValue)

  end subroutine testGetReal64

  subroutine testGetLogical()
    type (KeyValuePair_type) :: pair
    logical :: expectedValue, foundValue

    expectedValue = .true.
    foundValue = .false.
    pair = KeyValuePair('aKey', expectedValue)
    call getValue(pair, foundValue)
    call assertTrue(expectedValue .eqv. foundValue)

  end subroutine testGetLogical

  subroutine testGetString()
    type (KeyValuePair_type) :: pair
    character(len=10) :: expectedValue, foundValue

    expectedValue = 'happy'
    pair = KeyValuePair('aKey', expectedValue)
    call getValue(pair, foundValue)
    call assertEqual(expectedValue, foundValue)

  end subroutine testGetString

  subroutine testGetDictionary()
    type (KeyValuePair_type) :: pair
    type (Dictionary_type) :: expectedValue, foundValue
    character(len=MAX_LEN_KEY), pointer :: expectedKeys(:), foundKeys(:)
    integer :: m

    expectedValue = Dictionary()
    call insert(expectedValue, 'wordA', 1)
    call insert(expectedValue, 'wordB', 2)

    pair = KeyValuePair('name', expectedValue)
    call getValue(pair, foundValue)

    call lookup(foundValue, 'wordA', m)
    call assertEqual(1, m)
    call lookup(foundValue, 'wordB', m)
    call assertEqual(2, m)

    call clean(foundValue)
    call clean(expectedValue)

  end subroutine testGetDictionary

  subroutine testGetIntegerArray()
    type (KeyValuePair_type) :: pair
    integer, dimension(2) :: expectedValues, foundValues

    expectedValues = [1,2]
    pair = KeyValuePair('aKey', expectedValues)
    call getValue(pair, foundValues)
    call assertEqual(expectedValues, foundValues)

  end subroutine testGetIntegerArray

  subroutine testGetReal64Array()
    type (KeyValuePair_type) :: pair
    real(kind=r64), dimension(2) :: expectedValues, foundValues

    expectedValues = [1.2,2.3]
    pair = KeyValuePair('aKey', expectedValues)
    call getValue(pair, foundValues)
    call assertEqual(expectedValues, foundValues)

  end subroutine testGetReal64Array

  subroutine testGetLogicalArray()
    type (KeyValuePair_type) :: pair
    integer, parameter :: NUM_ENTRIES = 2
    logical, dimension(NUM_ENTRIES) :: expectedValues, foundValues

    integer :: i

    expectedValues = [ .true., .false.]
    pair = KeyValuePair('aKey', expectedValues)
    call getValue(pair, foundValues)
    do i = 1, NUM_ENTRIES
      call assertTrue(expectedValues(i) .eqv. foundValues(i))
    end do

  end subroutine testGetLogicalArray

  subroutine testGetStringArray()
    type (KeyValuePair_type) :: pair
    integer, parameter :: NUM_ENTRIES = 2
    character(len=20), dimension(NUM_ENTRIES) :: expectedValues, foundValues
    
    integer :: i

    expectedValues = ['hello', 'goodbye']
    pair = KeyValuePair('aKey', expectedValues)
    call getValue(pair, foundValues)
    do i = 1, NUM_ENTRIES
      call assertEqual(expectedValues(i), foundValues(i))
    end do

  end subroutine testGetStringArray

  subroutine testWrongType()
    type (KeyValuePair_type) :: pair
    integer :: integerValue
    real(kind=r64) :: realValue
    logical :: logicalValue
    character(len=10) :: stringValue

    integerValue = 1
    pair = KeyValuePair('aKey', integerValue)
    call getValue(pair, realValue)
    if (.not. catch('Incorrect type for specified key: <aKey>')) then
      call throw(Exception('Failed to detect incorrect type in dictionary.'))
    end if

    call getValue(pair, logicalValue) ! real
    if (.not. catch('Incorrect type for specified key: <aKey>')) then
      call throw(Exception('Failed to detect incorrect type in dictionary.'))
    end if

    call getValue(pair, stringValue) ! real
    if (.not. catch('Incorrect type for specified key: <aKey>')) then
      call throw(Exception('Failed to detect incorrect type in dictionary.'))
    end if

    realValue = 1.2
    pair = KeyValuePair('aKey', realValue)
    call getValue(pair, integerValue) ! real
    if (.not. catch('Incorrect type for specified key: <aKey>')) then
      call throw(Exception('Failed to detect incorrect type in dictionary.'))
    end if


  end subroutine testWrongType

end module testKeyValuePair_mod
