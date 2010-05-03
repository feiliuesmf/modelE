module testStringUtilities_mod
  use pFUnit
  use StringUtilities_mod
  implicit none

contains

  subroutine testToLowerCase()
    call assertEqual('abcdef123abcdef',toLowerCase('AbCdEf123aBcDeF'))
  end subroutine testToLowerCase

end module testStringUtilities_mod

