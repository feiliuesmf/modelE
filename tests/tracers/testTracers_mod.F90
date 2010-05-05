module testTracers_mod
  use pFUnit
  use Tracers_mod
  implicit none

  public :: setUp
  public :: tearDown
  public :: fixture

  type fixture
    integer :: unit
    integer :: numTracers
  end type fixture

contains

  subroutine setUp(this)
    use FileManager
    type (fixture) :: this
    integer :: unit

    call openUnit('testTracers.txt', unit, qold=.false., qbin=.false.)
    this%unit = unit

    write(unit,*)'{'
    write(unit,*)' name = speciesA'
    write(unit,*)' molecularMass = 1.2 ! g'
    write(unit,*)' optionA = True'
    write(unit,*)'}'

    write(unit,*)'{'
    write(unit,*)' name = speciesB'
    write(unit,*)' molecularMass = 2.3 ! g'
    write(unit,*)' optionB = False'
    write(unit,*)'}'

    rewind(unit)

    this%numTracers = 2
  end subroutine setUp

  subroutine tearDown(this)
    type (fixture) :: this

    close(this%unit, status='delete')
    this%unit = -1
  end subroutine tearDown

  subroutine testReadFromText(this)
    type (fixture) :: this

    type (Tracer_type), pointer :: tracers(:)
    integer :: unit

    unit = this%unit
    tracers => readFromText(unit)

    call assertEqual(this%numTracers, size(tracers), &
         & 'Incorrect number of tracers found.')

  end subroutine testReadFromText

  subroutine testGetIndex(this)
    type (fixture) :: this

    type (Tracer_type), pointer :: tracers(:) => null()
    integer :: unit

    unit = this%unit
    tracers => readFromText(unit)

    call assertEqual(1, getIndex(tracers, 'speciesA'))
    call assertEqual(2, getIndex(tracers, 'speciesB'))
    call assertEqual(NOT_FOUND, getIndex(tracers, 'speciesZ'))

  end subroutine testGetIndex

  subroutine testSetProperty(this)
    use GenericType_mod
    use Dictionary_mod
    type (fixture) :: this

    type (Tracer_type), pointer :: tracers(:) => null()
    integer :: unit
    real(kind=r64) :: expected, found
    type (Dictionary_type), pointer :: properties

    unit = this%unit
    tracers => readFromText(unit)

    expected = 2.1d+0
    properties => getProperties(tracers(1))
    call insert(properties, 'otherProperty', expected)

    found = getProperty(tracers(1), 'otherProperty')
    call assertEqual(expected, found)

  end subroutine testSetProperty

  subroutine testSetPropertyMulti(this)
    use GenericType_mod
    use Dictionary_mod
    type (fixture) :: this

    type (Tracer_type), pointer :: tracers(:) => null()
    integer :: unit
    real(kind=r64) :: expected, found
    type (Dictionary_type), pointer :: properties

    unit = this%unit
    tracers => readFromText(unit)

    expected = 2.1d+0
    properties => getProperties(tracers, 'speciesB')
    call insert(properties, 'otherProperty', expected)

    found = getProperty(tracers, 'speciesB', 'otherProperty')
    call assertEqual(expected, found)

  end subroutine testSetPropertyMulti

  subroutine testGetPropertySingle(this)
    use GenericType_mod
    type (fixture) :: this

    type (Tracer_type), pointer :: tracers(:) => null()
    integer :: unit
    real(kind=r64) :: expected, found

    unit = this%unit
    tracers => readFromText(unit)

    expected = 1.2d+0
    found = getProperty(tracers(1), 'molecularMass')
    call assertEqual(expected, found)

  end subroutine testGetPropertySingle

  subroutine testGetPropertyMulti(this)
    use GenericType_mod
    type (fixture) :: this

    type (Tracer_type), pointer :: tracers(:) => null()
    integer :: unit
    real(kind=r64) :: expected, found

    unit = this%unit
    tracers => readFromText(unit)

    expected = 1.2d+0
    found = getProperty(tracers, 'speciesA', 'molecularMass')
    call assertEqual(expected, found)

  end subroutine testGetPropertyMulti


end module testTracers_mod
