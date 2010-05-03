module testGenericType_mod
  use pFUnit
  use GenericType_mod
  implicit none

contains

  subroutine testGetInteger()
    integer :: found

    found = GenericType(1)
    call assertEqual(1, found)

    call assertTrue(1 == GenericType(1))
    call assertTrue(2 == GenericType(2))
    call assertFalse(1 == GenericType(2))
  end subroutine testGetInteger

  subroutine testGetIntegerArray()
    integer, parameter :: NUM_VALUES=2
    integer :: expected(NUM_VALUES), found(NUM_VALUES)

    expected = [1,2]
    found = GenericType(expected)
    call assertEqual(expected, found)
    call assertTrue(all(expected == GenericType(expected)))

  end subroutine testGetIntegerArray

  subroutine testGetReal64()
    real(kind=r64) :: expected, found

    expected = 1.234d+0
    found = GenericType(expected)
    call assertEqual(expected, found)
    call assertTrue(expected == GenericType(expected))
    call assertFalse(expected-1 == GenericType(expected))

  end subroutine testGetReal64

  subroutine testToStringInteger()
    call assertEqual('1', toString(GenericType(1)))
    call assertEqual('2', toString(GenericType(2)))
  end subroutine testToStringInteger

  subroutine testToStringReal64()
    character(len=MAX_LEN_VALUE) :: string
    real(kind=r64) :: expectedValue, foundValue

    expectedValue = 321.1234567890123456789d+0
    string = toString(GenericType(expectedValue))
    read(string,'(g)') foundValue
    call assertEqual(expectedValue, foundValue)

  end subroutine testToStringReal64

  subroutine testFromStringInteger()
    call assertTrue(1 == fromString('1', INTEGER_TYPE))
    call assertTrue(1234 == fromString(' 1234', INTEGER_TYPE))
  end subroutine testFromStringInteger
    
  subroutine testFromStringReal64()
    real(kind=r64) :: expected
    type (GenericType_type) :: generic

    expected = 1.23d+0
    generic = GenericType(expected)
    call assertTrue(expected == fromString(toString(generic), REAL64_TYPE))

    expected = -123.789d+0
    generic = GenericType(expected)
    call assertTrue(expected == fromString(toString(generic), REAL64_TYPE))

  end subroutine testFromStringReal64

  subroutine testFailFromString()
    integer :: i

    i = fromString('1.234', INTEGER_TYPE)
    call assertFailedAssert('GenericType::fromString() - cannot convert string "1.234" to integer.', &
         & 'Failed to detect conversion error.')

    i = fromString('1.234', -1)
    call assertFailedAssert('GenericType::fromString() - no such type.')

  end subroutine testFailFromString
    
end module testGenericType_mod
