module testParser_mod
  use pFUnit
  use Parser_mod
  implicit none
  private

  ! legacy
  public :: testStripComment_noComment
  public :: testStripComment
  public :: testStripComment_1stChar

  public :: testStripComment_new
  public :: testSkipHeader
  public :: testNoEndHeader
  public :: testIsEndOfList

  public :: testGetValueType
  public :: testGetCommonValueType

  public :: testSplitTokensA
  public :: testSplitTokensB
  public :: testSplitTokensC
  public :: testSplitTokensD
  public :: testSplitEmbeddedComma

  public :: testParse

contains

  subroutine testStripComment_noComment()
    character(len=30) :: expectedString
    character(len=30) :: string

    string = 'has no comment characters'
    expectedString = string
    call assertEqual(expectedString, strip_comment(string))

  end subroutine testStripComment_noComment

  subroutine testStripComment()
    character(len=30) :: expectedString
    character(len=30) :: string

    string = 'comment starts here ! ignore this'
    expectedString = 'comment starts here'
    
    call assertEqual(expectedString, strip_comment(string))

    string = 'comment starts #'
    expectedString = 'comment starts'
    call assertEqual(expectedString, strip_comment(string))

  end subroutine testStripComment

  subroutine testStripComment_1stChar()
    character(len=30) :: expectedString
    character(len=30) :: string

    string = '! all comment'
    expectedString = ''
    call assertEqual(expectedString, strip_comment(string))

  end subroutine testStripComment_1stChar

  subroutine testStripComment_new()
    type (Parser_type) :: parser
    character(len=40) :: expectedString
    character(len=40) :: string

    call setCommentCharacters(parser, ':')
    string = 'has no comment characters? ! comment'
    expectedString = string
    call assertEqual(expectedString, stripComment(parser, string))

    call setCommentCharacters(parser, ' ')
    expectedString = 'has'
    call assertEqual(expectedString, stripComment(parser, string))

    call setCommentCharacters(parser, '#!')
    expectedString = 'has no comment characters?'
    call assertEqual(expectedString, stripComment(parser, string))

  end subroutine testStripComment_new

  subroutine testSkipHeader()
    use FileManager
    type (Parser_type) :: parser
    integer :: unit
    character(len=*), parameter :: END_HEADER = '*** end header ***'
    character(len=256) :: line

    call openUnit('testParser.txt', unit, qold=.false., qbin=.false.)
    write(unit,*) 'a'
    write(unit,*) 'b'
    write(unit,*) END_HEADER
    write(unit,*) 'first line after header'
    rewind(unit)

    call setEndHeader(parser, END_HEADER)
    call skipHeader(parser, unit)
    read(unit,'(a257)') line
    call assertEqual('first line after header', line)

    close(unit, status='delete')

  end subroutine testSkipHeader

  subroutine testNoEndHeader()
    use FileManager
    type (Parser_type) :: parser
    integer :: unit
    character(len=*), parameter :: END_HEADER = '*** end header ***'
    character(len=256) :: line

    call openUnit('testParser.txt', unit, qold=.false., qbin=.false.)
    write(unit,*) 'a'
    write(unit,*) 'b'
    rewind(unit)

    call setEndHeader(parser, END_HEADER)
    call skipHeader(parser, unit)
    if (.not. catch('end of file before reaching header')) then
      call throw(Exception('did not notice missing end of header'))
    end if

    close(unit, status='delete')
  end subroutine testNoEndHeader

  subroutine testGetValueType()
    use Dictionary_mod

    call assertEqual(INTEGER_CASE, getValueType('3'))
    call assertEqual(INTEGER_CASE, getValueType('123456'))

    call assertEqual(REAL64_CASE, getValueType('1.'))
    call assertEqual(REAL64_CASE, getValueType('1.2'))
    call assertEqual(REAL64_CASE, getValueType('1.2e+10'))
    call assertEqual(REAL64_CASE, getValueType('234.6d-12'))

    call assertEqual(LOGICAL_CASE, getValueType('.true.'))
    call assertEqual(LOGICAL_CASE, getValueType('.false.'))
    call assertEqual(LOGICAL_CASE, getValueType('T'))
    call assertEqual(LOGICAL_CASE, getValueType('F'))
    call assertEqual(LOGICAL_CASE, getValueType('t'))
    call assertEqual(LOGICAL_CASE, getValueType('f'))
    call assertEqual(LOGICAL_CASE, getValueType('t'))
    call assertEqual(LOGICAL_CASE, getValueType('false'))
    call assertEqual(LOGICAL_CASE, getValueType('true'))

    call assertEqual(STRING_CASE, getValueType('hello'))
    call assertEqual(STRING_CASE, getValueType('truely'))
    call assertEqual(STRING_CASE, getValueType('Fa'))
    call assertEqual(STRING_CASE, getValueType('''.true.''')) ! has quotes => string
    call assertEqual(STRING_CASE, getValueType('".true."')) ! has quotes => string

  end subroutine testGetValueType

  subroutine testGetCommonValueType()
!@sum Test that a mixed set of tokens are treated as a category that
!@+ sufficiently broad.  E.g mix of integer and string should be
!@+ treated as string.
    use Dictionary_mod

    call assertEqual(INTEGER_CASE, getValueType(['3','4']))
    call assertEqual(REAL64_CASE, getValueType(['3.','4 ']))
    call assertEqual(REAL64_CASE, getValueType(['3 ','4.', '5 ']))
    call assertEqual(REAL64_CASE, getValueType(['3 ','4.']))
    call assertEqual(REAL64_CASE, getValueType(['3.','4.']))

    call assertEqual(LOGICAL_CASE, getValueType(['T     ','F     ','.true.']))
    call assertEqual(STRING_CASE,  getValueType(['3     ','F     ','.true.']))

  end subroutine testGetCommonValueType



  subroutine testIsEndOfList()
    use FileManager
    type (Parser_type) :: parser
    character(len=*), parameter :: END_LIST = '*** end params ***'

    call setEndOfList(parser, END_LIST)
    call assertTrue(isEndOfList(parser, END_LIST))
    call assertFalse(isEndOfList(parser, 'some other string'))
  end subroutine testIsEndOfList

  subroutine testSplitTokensA()
    type (Parser_type) :: parser
    character(len=*), parameter :: separators = '=,'
    character(len=128) :: string
    character(len=MAX_LEN_TOKEN), pointer :: tokens(:)

    call setTokenSeparators(parser, separators)

    tokens => splitTokens(parser, '  ')
    call assertEqual(0, size(tokens))
    deallocate(tokens)

  end subroutine testSplitTokensA

  subroutine testSplitTokensB()
    type (Parser_type) :: parser
    character(len=*), parameter :: separators = '=,'
    character(len=128) :: string
    character(len=MAX_LEN_TOKEN), pointer :: tokens(:)

    call setTokenSeparators(parser, separators)

    tokens => splitTokens(parser, '  hello  ')
    call assertEqual(1, size(tokens))
    if (size(tokens) == 1) call assertEqual('hello', tokens(1))
    deallocate(tokens)

  end subroutine testSplitTokensB

  subroutine testSplitTokensC()
    type (Parser_type) :: parser
    character(len=*), parameter :: separators = '=,'
    character(len=128) :: string
    character(len=MAX_LEN_TOKEN), pointer :: tokens(:)

    call setTokenSeparators(parser, separators)

    tokens => splitTokens(parser, ' var = 3  ')
    call assertEqual(2, size(tokens))
    if (size(tokens) == 2) then
      call assertEqual('var', tokens(1))
      call assertEqual('3', tokens(2))
    end if
    deallocate(tokens)

  end subroutine testSplitTokensC

  subroutine testSplitTokensD()
    type (Parser_type) :: parser
    character(len=*), parameter :: separators = '=,'
    character(len=128) :: string
    character(len=MAX_LEN_TOKEN), pointer :: tokens(:)

    call setTokenSeparators(parser, separators)

    tokens => splitTokens(parser, ' var = "hello", 3, 7.891e+10, .true., 123456  ')
    call assertEqual(6, size(tokens))
    if (size(tokens) == 6) then
      call assertEqual('var', tokens(1))
      call assertEqual('"hello"', tokens(2))
      call assertEqual('3', tokens(3))
      call assertEqual('7.891e+10', tokens(4))
      call assertEqual('.true.', tokens(5))
      call assertEqual('123456', tokens(6))
    end if
    deallocate(tokens)

  end subroutine testSplitTokensD

  subroutine testSplitEmbeddedComma()
    type (Parser_type) :: parser
    character(len=*), parameter :: separators = '=,'
    character(len=128) :: string
    character(len=MAX_LEN_TOKEN), pointer :: tokens(:)

    call setTokenSeparators(parser, separators)

    tokens => splitTokens(parser, " var = 'short name', 'name with , in it'")
    call assertEqual(3, size(tokens))
    if (size(tokens) == 3) then
      call assertEqual('var', tokens(1))
      call assertEqual("'short name'", tokens(2))
      call assertEqual("'name with , in it'", tokens(3))
    end if
    deallocate(tokens)

  end subroutine testSplitEmbeddedComma

  subroutine testParse()
    use FileManager
    use Dictionary_mod
    type (Parser_type) :: parser

    type (Dictionary_type) :: aDictionary
    integer :: unit
    character(len=*), parameter :: END_HEADER = '*** end header ***'
    character(len=*), parameter :: END_LIST = '*** end params ***'
    character(len=256) :: line

    integer :: betaValue
    real(kind=r64) :: gammaValue
    logical :: deltaValues(4)
    call openUnit('testParser.txt', unit, qold=.false., qbin=.false.)
    write(unit,*) 'alhap = 0' ! ignore - in header
    write(unit,*) 'b'
    write(unit,*) END_HEADER
    write(unit,*) '! alpha = 1' ! ignore - only a comment
    write(unit,*) 'beta = 2'
    write(unit,*) 'gamma = 3.'
    write(unit,*) '! alpha = 2' ! ignore - only a comment
    write(unit,*) 'delta = T, F, T, T'
    write(unit,*) END_LIST
    write(unit,*) '! alpha = 3' ! ignore - after header
    rewind(unit)

    call setEndHeader(parser, END_HEADER)
    call setEndOfList(parser, END_LIST)
    call setTokenSeparators(parser, '=,')
    call setCommentCharacters(parser, '!#')

    aDictionary = parse(parser, unit)
    close(unit, status='delete')

    call assertFalse(hasKey(aDictionary, 'alpha'), 'There should not be an "alpha".')
    call assertTrue(hasKey(aDictionary, 'beta'), 'There should be a "beta".')
    call assertTrue(hasKey(aDictionary, 'gamma'), 'There should be a "gamma".')
    call assertTrue(hasKey(aDictionary, 'delta'), 'There should be a "delta".')

    call lookup(aDictionary, 'beta', betaValue)
    call assertEqual(2, betaValue, 'Wrong value for "beta".')

    call lookup(aDictionary, 'gamma', gammaValue)
    call assertEqual(3.d+0, gammaValue, 'Wrong value for "gamma".')

    call lookup(aDictionary, 'delta', deltaValues)
    call assertTrue(all([.true.,.false.,.true.,.true.] .eqv. deltaValues), 'Wrong value for "delta".')

    call clean(aDictionary)

  end subroutine testParse

end module testParser_mod
