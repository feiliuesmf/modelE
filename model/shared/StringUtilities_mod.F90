module StringUtilities_mod
  implicit none
  private

  public :: toLowerCase

contains

  function toLowerCase(string) result(newString)
    character(len=*), intent(in) :: string
    character(len=len_trim(string)) :: newString

    integer n, i
    integer A, Z, shift, c

    A = iachar( 'A' )
    Z = iachar( 'Z' )
    shift = iachar( 'a' ) - iachar( 'A' )

    newString = trim(string)
    n = len(newString)
    do i=1,n
      c = iachar( newString(i:i) )
      if ( c>=A .and. c<=Z ) newString(i:i) = achar( c + shift )
    enddo

  end function toLowerCase

end module StringUtilities_mod
