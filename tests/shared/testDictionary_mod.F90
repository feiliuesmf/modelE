module testDictionary_mod
  use pfunit
  use Dictionary_mod
  implicit none
  private

  public :: testSetInteger
  public :: testSetReal
  public :: testSetString

contains

  subroutine testSetInteger()
    integer :: expected
    integer :: found

    expected = 1
    call set_param('testParam_Int', expected)
    call get_param('testParam_Int', found)
    call assertEqual(expected, found)

  end subroutine testSetInteger

  subroutine testSetReal()
    real(kind=r64) :: expected
    real(kind=r64) :: found

    expected = 1.234
    call set_param('testParam_Real', expected)
    call get_param('testParam_Real', found)
    call assertEqual(expected, found)
    
  end subroutine testSetReal

  subroutine testSetString()
    character(len=20) :: expected
    character(len=20) :: found

    expected = 'testValue'
    call set_param('testParam_String', expected)
    call get_param('testParam_String', found)
    call assertEqual(expected, found)
    
  end subroutine testSetString

end module testDictionary_mod
