module testKeyValuePair_mod
  use pFUnit
  use GenericType_mod
  use KeyValuePair_mod
  implicit none
  private

  public :: testGetKey
  public :: testGetNumValues
  public :: testGetValue

contains

  subroutine testGetKey()
    type (KeyValuePair_type) :: pair

    pair = KeyValuePair('aKey', GenericType(1))
    call assertEqual('aKey', getKey(pair)) 

  end subroutine testGetKey

  subroutine testGetNumValues()
    type (KeyValuePair_type) :: pair

    pair = KeyValuePair('aKey', GenericType(1))
    call assertEqual(1, getNumValues(pair))

    pair = KeyValuePair('aKey', GenericType(1.234d+0))
    call assertEqual(1, getNumValues(pair))

    pair = KeyValuePair('aKey', GenericType(.true.))
    call assertEqual(1, getNumValues(pair))

    pair = KeyValuePair('aKey', GenericType('hello'))
    call assertEqual(1, getNumValues(pair))

    ! lists
    pair = KeyValuePair('aKey', GenericType([1,2,3]))
    call assertEqual(3, getNumValues(pair))

    pair = KeyValuePair('aKey', GenericType([1.234d+0]))
    call assertEqual(1, getNumValues(pair))

    pair = KeyValuePair('aKey', GenericType([.true.,.false.,.true.]))
    call assertEqual(3, getNumValues(pair))

    pair = KeyValuePair('aKey', GenericType(['hello', 'goodbye']))
    call assertEqual(2, getNumValues(pair))

  end subroutine testGetNumValues

  subroutine testGetValue()
    integer :: expected(2)

    call assertTrue(1 == getValue(KeyValuePair('aKey', GenericType(1))))

    expected = [1,2]
    call assertTrue(all(expected == getValues(KeyValuePair('aKey', GenericType(expected)))))

  end subroutine testGetValue

end module testKeyValuePair_mod
