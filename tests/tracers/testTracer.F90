module TestTracer_mod
  use pFUnit
  use Tracer_mod
  use Dictionary_mod, only: Dictionary
  implicit none
  private

  public :: testHasProperty
  public :: testSetProperty

contains

  subroutine testSetProperty()
    use GenericType_mod
    use Dictionary_mod

    real(kind=r64) :: expected
    real(kind=r64), pointer :: found
    type (Dictionary), pointer :: properties
    type (Tracer), pointer :: aTracer

    aTracer => newTracer('speciesA')
    expected = 1.0
    call aTracer%insert('otherProperty', expected)
    allocate(found)
    found = aTracer%getReference('otherProperty')
    call assertEqual(expected, found)
    call clean(aTracer)

  end subroutine testSetProperty

  subroutine testHasProperty()
    use GenericType_mod

    real(kind=r64) :: expected, found
    type (Tracer), pointer :: aTracer

    aTracer => newTracer('speciesA')
    call aTracer%insert('optionA', 'description')
    call assertFalse(aTracer%has('optionB'))

  end subroutine testHasProperty

end module TestTracer_mod
