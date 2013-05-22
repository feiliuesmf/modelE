module testParser_mod
  use pFUnit
  implicit none
  private

  ! legacy
  public :: testStripComment_noComment
  public :: testStripComment
  public :: testStripComment_1stChar

  public :: testStripComment_new
  public :: testSkipHeader
  public :: testNoBeginData
  public :: testIsEndData

  public :: testSplitTokensA
  public :: testSplitTokensB
  public :: testSplitTokensC
  public :: testSplitTokensD
  public :: testSplitEmbeddedComma
  public :: testBadFirstSeparator

  public :: testParseEOF
  public :: testParseSimple
  public :: testParse
  public :: testParseAttribute
  public :: testParseNoValue

  public :: testWriteFormatted
  public :: testWriteTextExample


!TO DO - create tests for these:
!!$  public :: testMissingEndQuote
!!$  public :: testFirstSep ! '='
!!$  public :: testIgnoreEmbeddedComments

contains
  subroutine testStripComment_noComment()
    use Parser_mod
    character(len=30) :: expectedString
    character(len=30) :: string

    string = 'has no comment characters'
    expectedString = string
    call assertEqual(expectedString, strip_comment(string))

  end subroutine testStripComment_noComment

  subroutine testStripComment()
    use Parser_mod
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
    use Parser_mod
    character(len=30) :: expectedString
    character(len=30) :: string

    string = '! all comment'
    expectedString = ''
    call assertEqual(expectedString, strip_comment(string))

  end subroutine testStripComment_1stChar

  subroutine testStripComment_new()
    use Parser_mod
    type (Parser_type) :: parser
    character(len=40) :: expectedString
    character(len=40) :: string

    call setCommentCharacters(parser, ':')
    string = 'has no comment characters? ! comment'
    expectedString = string
    call assertEqual(expectedString, stripComment(parser, string))

    call setCommentCharacters(parser, ' ,')
    expectedString = 'has'
    call assertEqual(expectedString, stripComment(parser, string))

    call setCommentCharacters(parser, '#!')
    expectedString = 'has no comment characters?'
    call assertEqual(expectedString, stripComment(parser, string))

  end subroutine testStripComment_new

  subroutine testSkipHeader()
    use Parser_mod
    use FileManager
    type (Parser_type) :: parser
    integer :: unit
    integer :: status
    character(len=*), parameter :: BEGIN_DATA = '*** end header ***'
    character(len=MAX_LEN_LINE) :: line

    call openUnit('testParser.txt', unit, qold=.false., qbin=.false.)
    write(unit,'(a)') 'a'
    write(unit,'(a)') 'b'
    write(unit,'(a)') BEGIN_DATA
    write(unit,'(a)') 'first line after header'
    rewind(unit)

    call setBeginData(parser, BEGIN_DATA)
    call skipHeader(parser, unit, status)
    read(unit,'(a256)') line
    call assertEqual('first line after header', line)

    close(unit, status='delete')

  end subroutine testSkipHeader

  subroutine testNoBeginData()
    use Parser_mod
    use FileManager
    type (Parser_type) :: parser
    integer :: unit
    integer :: status
    character(len=*), parameter :: BEGIN_DATA = '*** end header ***'
    character(len=256) :: line

    call openUnit('testParser.txt', unit, qold=.false., qbin=.false.)
    write(unit,'(a)') 'a'
    write(unit,'(a)') 'b'
    rewind(unit)

    call setBeginData(parser, BEGIN_DATA)
    call skipHeader(parser, unit, status)
    call assertTrue(status /= 0)

    close(unit, status='delete')
  end subroutine testNoBeginData

  subroutine testIsEndData()
    use Parser_mod
    use FileManager
    type (Parser_type) :: parser
    character(len=*), parameter :: END_DATA = '*** end params ***'

    call setEndData(parser, END_DATA)
    call assertTrue(isEndData(parser, END_DATA))
    call assertFalse(isEndData(parser, 'some other string'))
  end subroutine testIsEndData

  subroutine testSplitTokensA()
    use Parser_mod
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
    use Parser_mod
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
    use Parser_mod
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
    use Parser_mod
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
    use Parser_mod
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

  subroutine testBadFirstSeparator()
    use Parser_mod
    type (Parser_type) :: parser
    character(len=*), parameter :: separators = '= ,'
    character(len=128) :: string
    character(len=MAX_LEN_TOKEN), pointer :: tokens(:)

#ifdef USE_PFUNIT
    call setTokenSeparators(parser, separators)

    tokens => splitTokens(parser, " my var = 3")
    call assertFailedAssert('Parser_mod: Illegal syntax.  "=" not first separator.', &
         & 'Failed to detect illegal syntax.')
    deallocate(tokens)
#endif

  end subroutine testBadFirstSeparator

  subroutine testParseEOF()
    use Parser_mod
    use FileManager
    use AttributeDictionary_mod
    type (Parser_type) :: parser
    type (AttributeDictionary) :: aDictionary

    integer :: status
    integer :: unit

    call openUnit('testParser.txt', unit, qold=.false., qbin=.false.)
    rewind(unit)

    aDictionary = parse(parser, unit, status)
    call assertTrue(status /= 0) ! EOF
    close(unit, status='delete')

  end subroutine testParseEOF

  subroutine testParseSimple()
    use Parser_mod
    use FileManager
    use AttributeDictionary_mod
    use Attributes_mod
    type (Parser_type) :: parser

    type (AttributeDictionary) :: aDictionary
    integer :: unit
    integer, pointer :: i(:)
    character(len=*), parameter :: BEGIN_DATA = '*** end header ***'
    character(len=*), parameter :: END_DATA = '*** end params ***'
    integer :: status

    call openUnit('testParser.txt', unit, qold=.false., qbin=.false.)
    write(unit,'(a)') BEGIN_DATA
    write(unit,'(a)') 'key = 1, 2, 3'
    write(unit,'(a)') END_DATA
    rewind(unit)

    call setBeginData(parser, BEGIN_DATA)
    call setEndData(parser, END_DATA)
    call setTokenSeparators(parser, '=,')
    call setCommentCharacters(parser, '!#')

    aDictionary = parse(parser, unit, status)
    close(unit, status='delete')

    call assertEqual(1, aDictionary%size())
    call assertTrue(aDictionary%has('key'), 'There should not be a "key".')
    i = aDictionary%getReference('key')
    call assertEqual([1,2,3], i, 'Wrong value for "key".')

    call clean(aDictionary)

  end subroutine testParseSimple

  subroutine testParse()
    use Parser_mod
    use FileManager
    use AttributeDictionary_mod
    use Attributes_mod
    type (Parser_type) :: parser

    type (AttributeDictionary) :: aDictionary
    integer :: unit
    character(len=*), parameter :: BEGIN_DATA = '*** end header ***'
    character(len=*), parameter :: END_DATA = '*** end params ***'

    integer :: status

    integer, pointer :: betaValue
    real(kind=r64), pointer :: gammaValue
    logical, pointer :: deltaValues(:)

    call openUnit('testParser.txt', unit, qold=.false., qbin=.false.)
    write(unit,'(a)') 'alpha = 0' ! ignore - in header
    write(unit,'(a)') 'b'
    write(unit,'(a)') BEGIN_DATA
    write(unit,'(a)') '! alpha = 1' ! ignore - only a comment
    write(unit,'(a)') 'beta = 2'
    write(unit,'(a)') ! empty line
    write(unit,'(a)') '!' ! effectively empty line
    write(unit,'(a)') 'gamma = 3.'
    write(unit,'(a)') '! alpha = 2' ! ignore - only a comment
    write(unit,'(a)') 'delta = T, F, T, T'
    write(unit,'(a)') END_DATA
    write(unit,'(a)') '! alpha = 3' ! ignore - after header
    rewind(unit)

    call setBeginData(parser, BEGIN_DATA)
    call setEndData(parser, END_DATA)
    call setTokenSeparators(parser, '=,')
    call setCommentCharacters(parser, '!#')

    aDictionary = parse(parser, unit, status)
    close(unit, status='delete')

    call assertEqual(3, aDictionary%size())
    call assertFalse(aDictionary%has('alpha'), 'There should not be an "alpha".')
    call assertTrue(aDictionary%has('beta'), 'There should be a "beta".')
    call assertTrue(aDictionary%has('gamma'), 'There should be a "gamma".')
    call assertTrue(aDictionary%has('delta'), 'There should be a "delta".')

    betaValue = aDictionary%getReference('beta')
    call assertEqual(2, betaValue, 'Wrong value for "beta".')

    gammaValue = aDictionary%getReference('gamma')
    call assertEqual(3.d+0, gammaValue, 'Wrong value for "gamma".')

    gammaValue = aDictionary%getReference('Gamma')
    call assertEqual(3.d+0, gammaValue, 'Wrong value for "Gamma".')

    deltaValues = aDictionary%getReference('delta')
    call assertTrue(all([.true.,.false.,.true.,.true.] .eqv. deltaValues), 'Wrong value for "delta".')

    call clean(aDictionary)

  end subroutine testParse

  subroutine testParseAttribute()
    use Parser_mod, only: Parser_type, parse
    use Parser_mod, only: setBeginData, setCommentCharacters, setEndData
    use Parser_mod, only: setTokenSeparators
    use FileManager
    use AttributeDictionary_mod
    use Attributes_mod
    type (Parser_type) :: parser

    type (AttributeDictionary) :: aDictionary
    integer :: unit
    character(len=*), parameter :: BEGIN_DATA = '*** end header ***'
    character(len=*), parameter :: END_DATA = '*** end params ***'

    integer :: status

    integer, pointer :: betaValue
    real(kind=r64), pointer :: gammaValue
    logical, pointer :: deltaValues(:)
    class (AbstractAttribute), pointer :: p

    call openUnit('testParser.txt', unit, qold=.false., qbin=.false.)
    write(unit,'(a)') 'alpha = 0' ! ignore - in header
    write(unit,'(a)') 'b'
    write(unit,'(a)') BEGIN_DATA
    write(unit,'(a)') '! alpha = 1' ! ignore - only a comment
    write(unit,'(a)') 'beta = 2'
    write(unit,'(a)') ! empty line
    write(unit,'(a)') '!' ! effectively empty line
    write(unit,'(a)') 'gamma = 3.'
    write(unit,'(a)') '! alpha = 2' ! ignore - only a comment
    write(unit,'(a)') 'delta = T, F, T, T'
    write(unit,'(a)') END_DATA
    write(unit,'(a)') '! alpha = 3' ! ignore - after header
    rewind(unit)

    call setBeginData(parser, BEGIN_DATA)
    call setEndData(parser, END_DATA)
    call setTokenSeparators(parser, '=,')
    call setCommentCharacters(parser, '!#')

    aDictionary = parse(parser, unit, status)

    close(unit, status='delete')


    call assertEqual(3, aDictionary%size())
    call assertFalse(aDictionary%has('alpha'), 'There should not be an "alpha".')
    call assertTrue(aDictionary%has('beta'), 'There should be a "beta".')
    call assertTrue(aDictionary%has('gamma'), 'There should be a "gamma".')
    call assertTrue(aDictionary%has('delta'), 'There should be a "delta".')

    betaValue = aDictionary%getReference('beta')
    call assertEqual(2, betaValue, 'Wrong value for "beta".')

    gammaValue = aDictionary%getReference('gamma')
    call assertEqual(3.d+0, gammaValue, 'Wrong value for "gamma".')

    deltaValues = aDictionary%getReference('delta')
    call assertTrue(all([.true.,.false.,.true.,.true.] .eqv. deltaValues), 'Wrong value for "delta".')

  end subroutine testParseAttribute

  subroutine testParseNoValue()
!@sum Tests that an exception is thrown if there is no value to associate with a key.
    use Parser_mod
    use FileManager
    use Dictionary_mod
    type (Parser_type) :: parser

    type (Dictionary) :: aDictionary
    integer :: unit
    character(len=*), parameter :: BEGIN_DATA = '*** end header ***'
    character(len=*), parameter :: END_DATA = '*** end params ***'
    character(len=256) :: line

    integer :: status

#ifdef USE_PFUNIT
    call openUnit('testParser.txt', unit, qold=.false., qbin=.false.)
    write(unit,'(a)') BEGIN_DATA
    write(unit,'(a)') 'beta'
    write(unit,'(a)') END_DATA
    rewind(unit)

    call setBeginData(parser, BEGIN_DATA)
    call setEndData(parser, END_DATA)
    call setTokenSeparators(parser, ' =,')
    call setCommentCharacters(parser, '!#')
    aDictionary = parse(parser, unit, status)
    call assertFailedAssert('Parser_mod: syntax error in input unit.', &
         & 'Failed to detect syntax error in parse')

    close(unit, status='delete')
#endif
  end subroutine testParseNoValue

  subroutine testWriteFormatted()
    use Parser_mod
    use FileManager
    use AttributeDictionary_mod
    type (Parser_type) :: parser

    type (AttributeDictionary) :: aDictionary
    type (AttributeDictionary) :: bDictionary
    integer :: unit
    character(len=*), parameter :: BEGIN_DATA = '*** end header ***'
    character(len=*), parameter :: END_DATA = '*** end params ***'
    character(len=256) :: line

    integer :: status

    integer, pointer :: i(:)

    aDictionary = newAttributeDictionary()
    call aDictionary%insert('key',[1,2,3])

    call setBeginData(parser, BEGIN_DATA)
    call setEndData(parser, END_DATA)
    call setTokenSeparators(parser, '=,')
    call setCommentCharacters(parser, '!#')
    call openUnit('testParser.txt', unit, qold=.false., qbin=.false.)
    call writeFormatted(parser, unit, aDictionary)
    close(unit)

    call openUnit('testParser.txt', unit, qold=.true., qbin=.false.)

    bDictionary = parse(parser, unit, status)

    i = bDictionary%getReference('key')

    call assertTrue(aDictionary%equals(bDictionary))

    call clean(bDictionary)
    call clean(aDictionary)

    close(unit, status='delete')

  end subroutine testWriteFormatted

  subroutine testWriteTextExample()
    use Parser_mod
    use FileManager
    use AttributeDictionary_mod
    type (Parser_type) :: parser

    type (AttributeDictionary) :: aDictionary
    type (AttributeDictionary) :: bDictionary
    integer :: unit
    character(len=*), parameter :: BEGIN_DATA = '{'
    character(len=*), parameter :: END_DATA = '}'
    character(len=256) :: line

    integer :: status

    integer :: scale
    real*8 :: mass
    character(len=MAX_LEN_TOKEN) :: name

    call setBeginData(parser, BEGIN_DATA)
    call setEndData(parser, END_DATA)
    call setTokenSeparators(parser, '=,')
    call setCommentCharacters(parser, '!#')
    call openUnit('testParser.txt', unit, qold=.false., qbin=.false.)

    aDictionary = newAttributeDictionary()
    call aDictionary%insert('name', 'Air')
    call aDictionary%insert('logScale', -2)
    call aDictionary%insert('molecularMass', 28.9655d0)
    call writeFormatted(parser, unit, aDictionary)
    
    rewind(unit)
    bDictionary = parse(parser, unit, status)

    call assertTrue(aDictionary%equals(bDictionary), 'Dictionaries differ.')

    call clean(bDictionary)
    call clean(aDictionary)

    close(unit, status='delete')

  end subroutine testWriteTextExample

end module testParser_mod
