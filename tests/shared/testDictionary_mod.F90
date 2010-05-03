module testDictionary_mod
  use pfunit
  use Dictionary_mod
  use KeyValuePair_mod
  implicit none
  private

  ! tests of legacy functionality
  public :: testSetInteger
  public :: testSetReal
  public :: testSetString

  public :: testSetIntegerList
  public :: testSetRealList
  public :: testSetStringList

  public :: testQuery
  public :: testReadWrite

  ! tests of new implementation
  public :: testKeyNotFound
  public :: testGetValueA
  public :: testGetValueB
  public :: testGetValueReal64

  public :: testGetNumEntries
  public :: testHasKey
  public :: testGetKeys
  
  public :: testFailDuplicateKey

contains

  subroutine testSetInteger()
    integer :: expected
    integer :: found

    call reset()
    expected = 1
    call set_param('testParam_Int', expected)
    call get_param('testParam_Int', found)
    call assertEqual(expected, found)

  end subroutine testSetInteger

  subroutine testSetReal()
    real(kind=r64) :: expected
    real(kind=r64) :: found

    call reset()
    expected = 1.234
    call set_param('testParam_Real', expected)
    call get_param('testParam_Real', found)
    call assertEqual(expected, found)
    
  end subroutine testSetReal

  subroutine testSetString()
    character(len=20) :: expected
    character(len=20) :: found

    call reset()
    expected = 'testValue'
    call set_param('testParam_String', expected)
    call get_param('testParam_String', found)
    call assertEqual(expected, found)
    
  end subroutine testSetString

  subroutine testSetIntegerList()
    integer, parameter :: NUM_ELEMENTS = 2
    integer :: expected(NUM_ELEMENTS)
    integer :: found(NUM_ELEMENTS)

    call reset()
    expected = [1,2]
    call set_param('testParam_IntList', expected, NUM_ELEMENTS)
    call get_param('testParam_IntList', found, NUM_ELEMENTS)
    call assertEqual(expected, found)

  end subroutine testSetIntegerList

  subroutine testSetRealList()
    integer, parameter :: NUM_ELEMENTS = 2
    real(kind=r64) :: expected(NUM_ELEMENTS)
    real(kind=r64) :: found(NUM_ELEMENTS)

    expected = [1.1, 2.2]
    call set_param('testParam_RealList', expected, NUM_ELEMENTS)
    call get_param('testParam_RealList', found, NUM_ELEMENTS)
    call assertEqual(expected, found)

  end subroutine testSetRealList

  subroutine testSetStringList()
    integer, parameter :: NUM_ELEMENTS = 2
    character*(128) :: expected(NUM_ELEMENTS)
    character*(128) :: found(NUM_ELEMENTS)
    integer :: i

    call reset()
    expected = ['hello ','goodbye']
    call set_param('testParam_StringList', expected, NUM_ELEMENTS)
    call get_param('testParam_StringList', found, NUM_ELEMENTS)
    do i = 1, NUM_ELEMENTS
      call assertEqual(expected(i), found(i))
    end  do

  end subroutine testSetStringList

  subroutine testQuery()
    integer, parameter :: NUM_ELEMENTS = 2
    integer :: values(NUM_ELEMENTS)

    character(len=20) :: expectedName, foundName
    character(len=1) :: ptype
    integer :: dim

    call reset()
    values = [1,2]
    expectedName = 'testquery'

    call set_param(expectedName, values, NUM_ELEMENTS)
    call query_param(1, foundName, dim, ptype)

    call assertEqual(expectedName, foundName, 'check order that tests are run')
    call assertEqual(NUM_ELEMENTS, dim)
    call assertEqual('i', ptype)

  end subroutine testQuery

  subroutine testReadWrite()
    use FileManager
    integer :: unit
    character(len=20) :: foundName
    integer :: dim
    character(len=1) :: pType

    call reset()
    call set_param('testReadWriteA', 2)
    call set_param('testReadWriteB', [3,4], 2)

    call openUnit('dictionary.dat',unit,qold=.false.,qbin=.true.)
    call write_param(unit)
    close(unit)

    call reset()
    ! start with one parameter just to make sure things were cleared
    ! before reading
    call set_param('testReadWriteC', 'aString')

    call openUnit('dictionary.dat',unit,qold=.true.,qbin=.true.)
    call read_param(unit, ovrwrt = .false.)
    call query_param(3, foundName, dim, pType)
    close(unit)

    call assertEqual('testreadwriteb', foundname)
    call assertEqual(2, dim)
    call assertEqual('i', pType)
    call deleteFile('dictionary.dat')

  end subroutine testReadWrite


  subroutine testKeyNotFound()
    type (Dictionary_type) :: aDictionary
    integer :: i

    aDictionary = Dictionary()
    call lookup(aDictionary, 'alpha', i)
    call assertFailedAssert('Key not found: <alpha>.', &
         & 'Failed to detect missing key.')

  end subroutine testKeyNotFound

  subroutine testGetValueA()
    type (Dictionary_type) :: aDictionary
    integer :: expected, found

    expected = 1
    aDictionary = Dictionary()
    call insert(aDictionary, 'key', expected)
    call lookup(aDictionary, 'key', found)
    call assertEqual(expected, found)

    call clean(aDictionary)

  end subroutine testGetValueA

  subroutine testGetValueB()
    type (Dictionary_type) :: aDictionary
    integer :: expectedA, foundA
    integer :: expectedB, foundB

    expectedA = 1
    expectedB = 2
    aDictionary = Dictionary()
    call insert(aDictionary, 'keyA', expectedA)
    call insert(aDictionary, 'keyB', expectedB)

    call lookup(aDictionary, 'keyB', foundB)
    call assertEqual(expectedB, foundB)

    call lookup(aDictionary, 'keyA', foundA)
    call assertEqual(expectedA, foundA)

    call clean(aDictionary)

  end subroutine testGetValueB

  subroutine testGetValueReal64()
    type (Dictionary_type) :: aDictionary
    real(kind=r64) :: expected, found

    expected = 1
    aDictionary = Dictionary()
    call insert(aDictionary, 'key', expected)
    call lookup(aDictionary, 'key', found)
    call assertEqual(expected, found)

    call clean(aDictionary)

  end subroutine testGetValueReal64

  subroutine testGetNumEntries()
    type (Dictionary_type) :: aDictionary
    real(kind=r64) :: expected, found

    expected = 1
    aDictionary = Dictionary()

    call assertEqual(0, getNumEntries(aDictionary))

    call insert(aDictionary, 'key1', expected)
    call assertEqual(1, getNumEntries(aDictionary))

    call insert(aDictionary, 'key2', expected)
    call assertEqual(2, getNumEntries(aDictionary))

  end subroutine testGetNumEntries

  subroutine testHasKey()
    type (Dictionary_type) :: aDictionary

    aDictionary = Dictionary()
    call assertFalse(hasKey(aDictionary, 'a'))

    call insert(aDictionary, 'a', 1)
    call assertTrue(hasKey(aDictionary, 'a'))
    call assertFalse(hasKey(aDictionary, 'b'))

    call insert(aDictionary, 'b', 1)
    call assertTrue(hasKey(aDictionary, 'a'))
    call assertTrue(hasKey(aDictionary, 'b'))

    call clean(aDictionary)

  end subroutine testHasKey

  subroutine testGetKeys()
    type (Dictionary_type) :: aDictionary
    character(len=MAX_LEN_KEY), pointer :: keys(:)

    aDictionary = Dictionary()
    call insert(aDictionary, 'a', 1)
    call insert(aDictionary, 'b', 1)

    keys => getKeys(aDictionary)
    call assertEqual(2, size(keys))
    call assertTrue(any(keys == 'b'))
    call assertTrue(any(keys == 'a'))

    call clean(aDictionary)
  end subroutine testGetKeys

  subroutine testFailDuplicateKey()
    type (Dictionary_type) :: aDictionary
    aDictionary = Dictionary()
    call insert(aDictionary, 'a', 1)
    call insert(aDictionary, 'a', 2)
    call assertFailedAssert('Dictionary: duplicate key - <a>.', &
         & 'Failed to diagnose duplicate key.')
    call clean(aDictionary)
  end subroutine testFailDuplicateKey

end module testDictionary_mod
