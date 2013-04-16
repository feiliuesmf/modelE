module TestAttributeDictionary_mod
  use pfunit
  use AttributeDictionary_mod
  use Attributes_mod

  implicit none

  type Fixture
    type (AttributeDictionary) :: dictionary
  end type Fixture

contains

  subroutine setUp(this)
    type (Fixture), intent(inout) :: this
    this%dictionary = newAttributeDictionary()
  end subroutine setUp

  subroutine tearDown(this)
    type (Fixture), intent(inout) :: this
    call clean(this%dictionary)
  end subroutine tearDown

  ! Test that class can insert intrinsic integer and cast to an integer pointer.
  subroutine testInsertIntegerAttribute(this)
    type (Fixture), intent(inout) :: this
    
    integer, pointer :: pInteger

    call this%dictionary%insert('A', 1)

    pInteger = this%dictionary%getReference('A')
    call assertEqual(1, pInteger)

  end subroutine testInsertIntegerAttribute

  subroutine testAssignAttributeReference(this)
    use AttributeReference_mod
    use IntegerAttribute_mod, only: assignment(=)
    type (Fixture), intent(inout) :: this
    
    type (AttributeReference) :: scalar
    integer, pointer :: i => null()

    type (IntegerAttribute) :: iAttr

    iAttr = newAttribute(1)
    
    call scalar%set(iattr)

    i = scalar
    call assertEqual(1, i)

  end subroutine testAssignAttributeReference

  subroutine testAssignAttributeReferenceMulti(this)
    use AttributeReference_mod
    type (Fixture), intent(inout) :: this
    
    type (AttributeReference), target :: vector(3)
    integer, pointer :: i(:)
    integer, pointer :: j

    type (IntegerAttribute) :: iAttr(3)

    iAttr(1) = newAttribute(1)
    iAttr(2) = newAttribute(2)
    iAttr(3) = newAttribute(5)

    call vector(1)%set(iAttr(1))
    call vector(2)%set(iAttr(2))
    call vector(3)%set(iAttr(3))

    i = vector

    call assertEqual([1,2,5], i)

  end subroutine testAssignAttributeReferenceMulti
  
  ! Test that class can insert/cast intrinsic 1D integer arrays
  subroutine testInsertInteger1DAttribute(this)
    type (Fixture), intent(inout) :: this
    
    integer, pointer :: pInteger(:)
    call this%dictionary%insert('A', [1,2])

    pInteger = this%dictionary%getReference('A')
    call assertEqual([1,2], pInteger)

  end subroutine testInsertInteger1DAttribute

  ! Test that class can insert intrinsic double and cast to an double pointer.
  subroutine testInsertRealDPAttribute(this)
    type (Fixture), intent(inout) :: this
    
    real*8, pointer :: pRealDP
    real*8, parameter :: VALUE = 1.23d0
    call this%dictionary%insert('A', VALUE)

    pRealDP = this%dictionary%getReference('A')
    call assertEqual(VALUE, pRealDP)

  end subroutine testInsertRealDPAttribute

  ! Test that class can insert/cast intrinsic 1D integer arrays
  subroutine testInsertRealDP1DAttribute(this)
    type (Fixture), intent(inout) :: this
    
    real*8, pointer :: pRealDP(:)
    real*8, parameter :: VALUES(3) = [1.2,3.4,5.6]
    call this%dictionary%insert('A', VALUES)

    pRealDP = this%dictionary%getReference('A')
    call assertEqual(VALUES, pRealDP)

  end subroutine testInsertRealDP1DAttribute

  ! Test that class can insert intrinsic integer and cast to an integer pointer.
  subroutine testInsertStringAttribute(this)
    type (Fixture), intent(inout) :: this
    
    character(len=MAX_LEN_ATTRIBUTE_STRING), pointer :: pString

    call this%dictionary%insert('A', 'foo')

    pString = this%dictionary%getReference('A')
    call assertEqual('foo', pString)

  end subroutine testInsertStringAttribute

  ! Test that class can insert intrinsic logical and cast to a logical pointer.
  subroutine testInsertLogicalAttribute(this)
    type (Fixture), intent(inout) :: this
    
    logical, pointer :: pLogical
    logical, parameter :: TRUE = .true.
    logical, parameter :: FALSE = .false.

    call this%dictionary%insert('A', TRUE)
    call this%dictionary%insert('B', FALSE)

    pLogical = this%dictionary%getReference('A')
    call assertTrue(pLogical)
 
    pLogical = this%dictionary%getReference('B')
    call assertFalse(pLogical)

  end subroutine testInsertLogicalAttribute

  subroutine testComplexCase(this)
    type (Fixture), intent(inout) :: this
    integer, pointer :: betaValue
    real(kind=r64), pointer :: gammaValue
    logical, pointer :: deltaValues(:)

    call this%dictionary%insert('beta', 2)
    call this%dictionary%insert('gamma', 3.d0)
    call this%dictionary%insert('delta', [.true., .false., .true., .true.])

    call assertEqual(3, this%dictionary%size())
    call assertFalse(this%dictionary%has('alpha'), 'There should not be an "alpha".')
    call assertTrue(this%dictionary%has('beta'), 'There should be a "beta".')
    call assertTrue(this%dictionary%has('gamma'), 'There should be a "gamma".')
    call assertTrue(this%dictionary%has('delta'), 'There should be a "delta".')

    betaValue = this%dictionary%getReference('beta')
    call assertEqual(2, betaValue, 'Wrong value for "beta".')

    gammaValue = this%dictionary%getReference('gamma')
    call assertEqual(3.d+0, gammaValue, 'Wrong value for "gamma".')

    gammaValue = this%dictionary%getReference('Gamma')
    call assertEqual(3.d+0, gammaValue, 'Wrong value for "Gamma".')

    deltaValues = this%dictionary%getReference('delta')
    call assertTrue(all([.true.,.false.,.true.,.true.] .eqv. deltaValues), 'Wrong value for "delta".')

  end subroutine testComplexCase
  
end module TestAttributeDictionary_mod
