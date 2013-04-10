module testTracerIO_mod
  use pFUnit
  use Tracer_mod
  use TracerBundle_mod
  use Dictionary_mod, only: Dictionary
  use AttributeDictionary_mod, only: AttributeDictionary
  implicit none
  private

  public :: testReadFromText
  public :: testWriteFormatted
  public :: testWriteUnformatted

  public :: setUp
  public :: tearDown
  public :: fixture


  type fixture
    type (TracerBundle) :: bundle
    type (Dictionary) :: defaultValues
    type (AttributeDictionary) :: newDefaultValues
    integer :: unit
    integer :: numTracers
  end type fixture

contains

  subroutine setUp(this)
    use Dictionary_mod, only: Dictionary
    use AttributeDictionary_mod
    use FileManager
    type (fixture) :: this
    integer :: unit

    call openUnit('testTracers.txt', unit, qold=.false., qbin=.false.)
    this%unit = unit

    write(unit,'(a)')'{'
    write(unit,'(a)')' name = speciesA'
    write(unit,'(a)')' molecularMass = 1.2 ! g'
    write(unit,'(a)')' optionA = True'
    write(unit,'(a)')' optionC = False'
    write(unit,'(a)')'}'
    write(unit,'(a)')' '
    write(unit,'(a)')'{'
    write(unit,'(a)')' name = speciesB'
    write(unit,'(a)')' molecularMass = 2.3 ! g'
    write(unit,'(a)')' optionB = False'
    write(unit,'(a)')'}'

    rewind(unit)

    this%numTracers = 2

    this%defaultValues = Dictionary()
    call this%defaultValues%insert('optionC', .true.)
    call this%defaultValues%insert('optionD',  1)

    this%newDefaultValues = newAttributeDictionary()
    call this%newDefaultValues%insert('optionC', .true.)
    call this%newDefaultValues%insert('optionD', 1)

    this%bundle = newTracerBundle()

  end subroutine setUp

  subroutine tearDown(this)
    use Dictionary_mod, only: clean
    type (fixture) :: this

    close(this%unit, status='delete')
    this%unit = -1
    call clean(this%defaultValues)
    call clean(this%newDefaultValues)
    call clean(this%bundle)

  end subroutine tearDown

  subroutine testReadFromText(this)
    type (fixture) :: this

    integer :: unit

    unit = this%unit
    this%bundle = readFromText(unit, this%newDefaultValues)

    call assertEqual(this%numTracers, this%bundle%size(), &
         & 'Incorrect number of tracers found.')

  end subroutine testReadFromText

  ! support routine
  subroutine readTracers(this)
    type (fixture) :: this
    this%bundle = readFromText(this%unit, this%newDefaultValues)
  end subroutine readTracers

  subroutine testWriteFormatted(this)
    use FileManager
    type (fixture) :: this
    type (TracerBundle) :: bundle
    integer :: unit

    call readTracers(this)

    call openUnit('testTracersOut.txt', unit, qold=.false., qbin=.false.)
    call this%bundle%writeFormatted(unit)
    rewind(unit)

    bundle = readFromText(unit)
    call assertTrue(bundle == this%bundle)

    close(unit, status ='delete')
    call clean(bundle)
    
  end subroutine testWriteFormatted

  subroutine testWriteUnformatted(this)
     use FileManager
    type (fixture) :: this
    type (TracerBundle) :: bundle
    integer :: unit

    call readTracers(this)

    call openUnit('testTracersOut.bin', unit, qold=.false., qbin=.true.)
    call this%bundle%writeUnformatted(unit)
    rewind(unit)

    bundle = readUnformattedBundle(unit)
    call assertTrue(bundle == this%bundle, 'should be ==')

    close(unit, status ='delete')
    call clean(bundle)
    
  end subroutine testWriteUnformatted

end module testTracerIO_mod
