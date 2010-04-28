module Parser_mod
!@sum procedures to read parameters from the rundeck into the database
!@auth I. Aleinov
!@ver 1.0     
  implicit none

  public :: Parser_type
  public :: parse
  public :: stripComment
  public :: skipHeader
  public :: isEndOfList
  public :: setCommentCharacters
  public :: setEndHeader
  public :: setEndOfList
  public :: getValueType
  public :: splitTokens

  public :: MAX_LEN_LINE
  public :: MAX_LEN_TOKEN
  integer, parameter :: MAX_COMMENT_CHARACTERS = 2
  integer, parameter :: MAX_TOKEN_SEPARATORS   = 2
  integer, parameter :: MAX_LEN_LINE  = 256
  integer, parameter :: MAX_LEN_TOKEN = 32
  character(len=*), parameter :: ENTIRE_LINE = '(a257)' ! MAX_LEN_LINE + 1 char
  character(len=*), parameter :: ENTIRE_TOKEN = '(a33)' ! MAX_LEN_TOKEN + 1 char

  type Parser_type
    integer :: numCommentCharacters = 2
    integer :: numTokenSeparators = 2
    character(len=MAX_COMMENT_CHARACTERS) :: commentCharacters = '!#' ! legacy default
    character(len=MAX_TOKEN_SEPARATORS) :: tokenSeparators = '=,'     ! legacy default
    character(len=MAX_LEN_LINE) :: endHeader = '&&PARAMETERS'         ! legacy default
    character(len=MAX_LEN_LINE) :: endOfList = '&&END_PARAMETERS'     ! legacy default
  end type Parser_type

  type (Parser_type) :: globalParser

  interface getValueType
    module procedure getValueType_single
    module procedure getValueType_multi
  end interface

contains

  function strip_comment( str ) result(newStr)
    ! remove comment at the end of the line. comments symbols: !#
    character*(*), intent(in) :: str
    character(len=len(str)) :: newStr
    
    newStr = stripComment(globalParser, str)

  end function strip_comment


  subroutine skip_junk( str )
    character*(*) str

    do while ( len_trim( str ) > 0 .and. scan( str, ' =,' ) == 1 )
      str = str(2:)
    enddo
    return
  end subroutine skip_junk


  subroutine sread_int( str, value )
    character*(*) str
    integer value
    integer n

    read ( str, * ) value

    ! remove chars till the next [ =,]
    n = scan( str, ' =,' )
    str = str(n+1:)

    call skip_junk( str )

    return
  end subroutine sread_int


  subroutine sread_real( str, value )
    character*(*) str
    real*8 value
    integer n

    read ( str, * ) value

    ! remove chars till the next [ =,]
    n = scan( str, ' =,' )
    str = str(n+1:)

    call skip_junk( str )

    return
  end subroutine sread_real


  subroutine sread_char( str, value )
    character*(*) str
    character*(*) value
    !character*256 tstr
    integer n, n1

    ! replace '=' with space if not quoted
    n1 = scan( str, '''' )
    n  = scan( str, '=' )
    if ( n>0 .and. ( n1==0 .or. n<n1 ) ) str(n:n) = ' '

    read ( str, * ) value

    if ( scan( str, '''' ) == 1 ) then  ! quoted string
      str = str(2:)
      n = scan( str, '''' )
      str = str(n+1:)
    else  ! remove chars till the next [ =,]
      n = scan( str, ' =,' ) 
      str = str(n+1:)
    endif

    call skip_junk( str )

    return
  end subroutine sread_char


  subroutine parse_params( kunit )
    use PARAM
    integer, parameter :: MAXDIM=64
    integer, intent(in) :: kunit
    character*256 bufs
    character*32 name
    character*1 type
    integer np
    integer ivars(MAXDIM)
    real*8 rvars(MAXDIM)
    character*64 cvars(MAXDIM)

    ! skip unrelated stuff
    do
      read( kunit, '(a256)', err=666, end=667 ) bufs
      if ( len_trim(bufs) < 1 ) cycle
      read( bufs, * ) name
      if ( name == '&&PARAMETERS' ) exit
    enddo

    do
      read( kunit, '(a256)', err=666, end=666 ) bufs

      if ( len_trim(bufs) < 1 ) cycle

      bufs = strip_comment( bufs )
      call skip_junk( bufs )

      if ( len_trim(bufs) < 1 ) cycle

      !read the name of the variable
      call sread_char( bufs, name )

      if ( name == '&&END_PARAMETERS' ) exit  ! end of list 

      if ( len_trim(bufs) < 1 ) then
        print *,'PARSER: no values were given to param: ', name
        call stop_model('PARSER error',255)
      endif

      ! now check the type of variables
      if ( scan( bufs, '''' ) > 0 ) then
        type = 'c'
      else if ( scan( bufs, '.' ) > 0 ) then
        type = 'r'
      else
        type = 'i'
      endif

      select case ( type )
      case ('i')
        np = 0
        do while ( len_trim(bufs) > 0 )
          np = np+1
          call sread_int( bufs, ivars(np) )
        end do
        call set_param( name, ivars, np, 'or' )
      case ('r')
        np = 0
        do while ( len_trim(bufs) > 0 )
          np = np+1
          call sread_real( bufs, rvars(np) )
        end do
        call set_param( name, rvars, np, 'or' )
      case ('c')
        np = 0
        do while ( len_trim(bufs) > 0 )
          np = np+1
          call sread_char( bufs, cvars(np) )
        end do
        call set_param( name, cvars, np, 'or' )
      end select

    enddo

    return
666 print *, 'PARSER: Error reading params'
    call stop_model( 'PARSER: Error reading params', 255 )
667 print *, 'PARSER: No &&PARAMETERS or &&END_PARAMETERS found'
    call stop_model( &
    &     'PARSER: No &&PARAMETERS or &&END_PARAMETERS found',255)
  end subroutine parse_params

  function parse(this, unit) result(aDictionary)
    use Dictionary_mod
    type (Parser_type), intent(in) :: this
    integer, intent(in) :: unit

    type (Dictionary_type) :: aDictionary
    integer :: status
    character(len=MAX_LEN_TOKEN), pointer :: tokens(:)
    character(len=MAX_LEN_LINE) :: line

    aDictionary = Dictionary()

    call skipHeader(this, unit)
    
    do
      read(unit,fmt=ENTIRE_LINE,iostat=status) line
      if (status /= 0) exit
      if (isEndOfList(this, line)) exit

      tokens => splitTokens(this, line)
      select case (getValueType(tokens(2:)))
      case (INTEGER_TYPE)
        call insert(aDictionary, tokens(1), readIntegers(tokens(2:)))
      case (REAL64_TYPE)
        call insert(aDictionary, tokens(1), readReal64s(tokens(2:)))
      case (LOGICAL_TYPE)
        call insert(aDictionary, tokens(1), readLogicals(tokens(2:)))
      case (STRING_TYPE)
        call insert(aDictionary, tokens(1), readStrings(tokens(2:)))
      end select

      deallocate(tokens)
    end do

  contains

    function readIntegers(tokens) result(values)
      character(len=*), intent(in) :: tokens(:)
      integer, pointer :: values(:)
      
      integer :: i, n

      n = size(tokens)
      allocate(values(n))
      do i = 1, n
        read(tokens(i),'(i)') values(i)
      end do
    end function readIntegers

    function readReal64s(tokens) result(values)
      character(len=*), intent(in) :: tokens(:)
      real*8, pointer :: values(:)
      
      integer :: i, n

      n = size(tokens)
      allocate(values(n))
      do i = 1, n
        read(tokens(i),'(g)') values(i)
      end do
    end function readReal64s

    function readLogicals(tokens) result(values)
      character(len=*), intent(in) :: tokens(:)
      logical, pointer :: values(:)
      
      integer :: i, n
      character(len=MAX_LEN_TOKEN) :: token

      n = size(tokens)
      allocate(values(n))
      do i = 1, n
        token = tokens(i)
        select case(toLowerCase(token))
        case ('t', 'true', '.true.')
          values(i) = .true.
        case ('f', 'false', '.false.')
          values(i) = .false.
        end select
      end do
    end function readLogicals

    function readStrings(tokens) result(values)
      character(len=*), intent(in) :: tokens(:)
      character(len=MAX_LEN_TOKEN), pointer :: values(:)
      
      integer :: i, n

      n = size(tokens)
      allocate(values(n))
      do i = 1, n
        read(tokens(i),fmt=ENTIRE_TOKEN) values(i)
      end do
    end function readStrings

  end function parse
  
  subroutine setCommentCharacters(this, commentCharacters)
    type (Parser_type), intent(inout) :: this
    character(len=*), intent(in) :: commentCharacters

    this%commentCharacters = commentCharacters
    this%numCommentCharacters = len(commentCharacters) ! do NOT use trim - might want blanks

  end subroutine setCommentCharacters

  function stripComment(this, str) result(newStr)
    ! remove comment at the end of the line.
    type (Parser_type), intent(in) :: this
    character*(*), intent(in) :: str
    character(len=len(str)) :: newStr

    integer :: n

    n = scan( str, this%commentCharacters(1:this%numCommentCharacters))
    select case (n)
    case (0)
      newStr = trim(str)
    case (1:)
      newStr = str(:n-1)
    end select

  end function stripComment

  subroutine setEndHeader(this, endHeader)
    type (Parser_type), intent(inout) :: this
    character(len=*), intent(in) :: endHeader
    this%endHeader = endHeader
  end subroutine setEndHeader

  subroutine setEndOfList(this, endOfList)
    type (Parser_type), intent(inout) :: this
    character(len=*), intent(in) :: endOfList
    this%endOfList = endOfList
  end subroutine setEndOfList

  subroutine skipHeader(this, unit)
    use pFUnit
    type (Parser_type), intent(in) :: this
    integer, intent(in) :: unit
    character(len=MAX_LEN_LINE) :: line
    integer :: status
    
    do
      read(unit,fmt=ENTIRE_LINE,iostat=status) line
#ifdef USE_PFUNIT
      if (status < 0) then
        call throw(Exception('end of file before reaching header'))
        return
      end if
#endif
      if (trim(line) == this%endHeader) exit
    end do

  end subroutine skipHeader

  logical function isInteger(string)
    character(len=*), intent(in) :: string
    integer :: integerValue
    integer :: status
    
    read(string,'(i)',iostat=status) integerValue
    isInteger = (status == 0)

  end function isInteger

  logical function isReal64(string)
    character(len=*), intent(in) :: string
    real*8 :: real64Value
    integer :: status
    
    read(string,'(g)',iostat=status) real64Value
    isReal64 = (status == 0)

  end function isReal64

  logical function isLogical(string)
!@sum Allow a variety of convenient formats for true/false values
    use Dictionary_mod, only: toLowerCase
    character(len=*), intent(in) :: string
    logical :: logicalValue
    integer :: status
    
    select case (trim(toLowerCase(string)))
    case ('t','f','true','false','.true.','.false.')
      isLogical = .true.
    case default
      isLogical = .false.
    end select

  end function isLogical

  integer function getValueType_single(string) result(type)
    use Dictionary_mod
    character(len=*), intent(in) :: string

    real*8  :: real64Value
    logical :: logicalValue
    character(len=MAX_LEN_VALUE)  :: stringValue
    integer :: status
  
    if (isInteger(string)) then
      type = INTEGER_TYPE
      return
    else if (isReal64(string)) then
      type = REAL64_TYPE
      return
    else if (isLogical(string)) then
      type = LOGICAL_TYPE
      return
    else
      type = STRING_TYPE
    end if

  end function getValueType_single

  integer function getValueType_multi(tokens) result(type)
    use Dictionary_mod
    character(len=*), intent(in) :: tokens(:)

    integer :: numTokens
    integer :: i
    integer, allocatable :: types(:)

    numTokens = size(tokens)
    allocate(types(numTokens))
    do i = 1, numTokens
      types(i) = getValueType(tokens(i))
    end do

    if (all(types == types(1))) then ! simple case
      type = types(1)
    else ! mixed type
      if (all((types == INTEGER_TYPE) .or. (types == REAL64_TYPE))) then
        type = REAL64_TYPE ! integers treated as subset of reals
      else
        type = STRING_TYPE ! most general category
      end if
    end if

    deallocate(types)

  end function getValueType_multi

  logical function isEndOfList(this, string)
    type (Parser_type), intent(in) :: this
    character(len=*), intent(in) :: string
    
    isEndOfList = (trim(this%endOfList) == trim(string))

  end function isEndOfList

  subroutine setTokenSeparators(this, tokenSeparators)
    type (Parser_type), intent(inout) :: this
    character(len=*), intent(in) :: tokenSeparators
    this%tokenSeparators = tokenSeparators
    this%numTokenSeparators = len(tokenSeparators)
  end subroutine setTokenSeparators

  function splitTokens(this, string) result(tokens)
    type (Parser_type), intent(in) :: this
    character(len=*), intent(in) :: string
    character(len=MAX_LEN_TOKEN), pointer :: tokens(:)

    character(len=len(string)) :: buffer
    integer :: numTokens
    integer :: i, idxStart
    integer :: idxNextSeparator

    numTokens = countTokens(this, string)
    allocate(tokens(numTokens))
    
    buffer = adjustl(string)
    i = 0
    do
      if (len_trim(buffer) == 0) exit ! done
      i = i + 1
      ! skip strings which might have embedded separator characters
      if (scan(buffer, '"') == 1) then
        idxStart = scan(buffer(2:),'"') + 1 + 1 ! include 1st char
      else if (scan(buffer, "'") == 1) then
        idxStart = scan(buffer(2:),"'") + 1 + 1 ! include 1st char
      else
        idxStart = 1
      end if
      idxNextSeparator = scan(buffer(idxStart:), this%tokenSeparators)
      if (idxNextSeparator > 0) then
        tokens(i) = trim(buffer(:idxStart+idxNextSeparator-2))
        buffer = adjustl(buffer(idxStart+idxNextSeparator:))
      else
        tokens(i) = trim(buffer) ! take the rest
        exit
      end if
    end do
    
  end function splitTokens

  integer function countTokens(this, string)
    type (Parser_type), intent(in) :: this
    character(len=*), intent(in) :: string

    character(len=len(string)) :: buffer
    integer :: numTokens
    integer :: idxNextSeparator, idxStart

    buffer = adjustl(string)
    numTokens = 0
    do
      if (len_trim(buffer) == 0) exit ! done
      numTokens = numTokens + 1
      ! skip strings which might have embedded separator characters
      if (scan(buffer, '"') == 1) then
        idxStart = scan(buffer(2:),'"') + 1 + 1 ! include 1st char
      else if (scan(buffer, "'") == 1) then
        idxStart = scan(buffer(2:),"'") + 1 + 1 ! include 1st char
      else
        idxStart = 1
      end if
      idxNextSeparator = scan(buffer(idxStart:), this%tokenSeparators)
      if (idxNextSeparator > 0) then
        buffer = adjustl(buffer(idxStart + idxNextSeparator:))
      else
        exit
      end if
    end do

    countTokens = numTokens

  end function countTokens

end module PARSER_MOD
