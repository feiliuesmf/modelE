module testAssociativeArray_mod
  use pfunit
  use Foo_mod
  use FooAssociativeArray_mod
  implicit none

  ! using a fixture to help enforce cleanup
  type Fixture
    type (FooAssociativeArray) :: dictionary
  end type Fixture

contains

  subroutine setUp(this)
    type (Fixture), intent(inout) :: this

    this%dictionary = newFooAssociativeArray()

  end subroutine setUp

  subroutine tearDown(this)
    type (Fixture), intent(inout) :: this
    call clean(this%dictionary)
  end subroutine tearDown

  ! An empty dictionary has 0 size
  subroutine testSize_empty(this)
    type (Fixture), intent(inout) :: this

    call assertEqual(0, this%dictionary%size())

  end subroutine testSize_empty

  ! Each insert() should increase the size by 1
  subroutine testSizeAfterInsert(this)
    type (Fixture), intent(inout) :: this

    call this%dictionary%insert('A', Foo(2))
    call assertEqual(1, this%dictionary%size())

    call this%dictionary%insert('B', Foo(3))
    call assertEqual(2, this%dictionary%size())

  end subroutine testSizeAfterInsert

  ! has() should return FALSE if key has not been used.
  subroutine testHas_missing(this)
    type (Fixture), intent(inout) :: this

    call assertFalse(this%dictionary%has('A'))

  end subroutine testHas_missing

  ! has() should return true only for keys that have been used.
  subroutine testHas(this)
    type (Fixture), intent(inout) :: this

    call this%dictionary%insert('A', Foo(1))

    call assertTrue(this%dictionary%has('A'))  ! A is in dictionary
    call assertFalse(this%dictionary%has('B')) ! B is not

  end subroutine testHas

  ! getReference() returns nonassociated pointer for keys that have not been used.
  subroutine testGetReference_notFound(this)
    type (Fixture), intent(inout) :: this
    class (Foo), pointer :: p

    call this%dictionary%insert('A', Foo(1))
    p => this%dictionary%getReference('B')
    call assertTrue(.not. associated(p))

  end subroutine testGetReference_notFound

  ! getReference() returns pointer to entry.
  subroutine testGetReference(this)
    type (Fixture), intent(inout) :: this
    class (Foo), pointer :: p

    call this%dictionary%insert('A', Foo(2))

    p => this%dictionary%getReference('A')
    if (.not. catch(preserve=.true.)) then
      call assertEqual(2, p%value)
    end if

  end subroutine testGetReference

  ! This test ensures that the correct entry is selected when there
  ! are multiple entries.
  subroutine testGetReferenceB(this)
    type (Fixture), intent(inout) :: this
    class (Foo), pointer :: p

    print*,__LINE__,__FILE__
    call this%dictionary%insert('A', Foo(1))
    call this%dictionary%insert('B', Foo(2))
    call this%dictionary%insert('C', Foo(3))
    print*,__LINE__,__FILE__

    p => this%dictionary%getReference('B')
    if (.not. catch(preserve=.true.)) then
      call assertEqual(2, p%value)
    end if

    p => this%dictionary%getReference('A')
    p => this%dictionary%getReference('C')

  end subroutine testGetReferenceB

  subroutine testMerge(this)
    type (Fixture), intent(inout) :: this
    type (FooAssociativeArray) :: otherDictionary

    call this%dictionary%insert('A', Foo(1))
    otherDictionary = newFooAssociativeArray()
    call otherDictionary%insert('B', Foo(2))
    
    call otherDictionary%merge(this%dictionary)
    call assertTrue(otherDictionary%has('A'))

  end subroutine testMerge

  subroutine testMergeFail(this)
    type (Fixture), intent(inout) :: this
    type (FooAssociativeArray) :: otherDictionary

    call this%dictionary%insert('A', Foo(1))
    call this%dictionary%insert('B', Foo(2))
    otherDictionary = newFooAssociativeArray()
    call otherDictionary%insert('B', Foo(2))
    
    call otherDictionary%merge(this%dictionary)
    call assertFailedAssert('AssociativeArray::merge() failed due to duplicate keys: <b>.')

 end subroutine testMergeFail

  subroutine testStress(this)
    type (Fixture), intent(inout) :: this

    integer :: i, j
    integer, parameter :: NUM_ENTRIES = 10
    integer, parameter :: NUM_CHAR = 16
    real :: x
    class (Foo), pointer :: f
    integer :: c0, c1, cr

    character(len=NUM_CHAR) :: key

    call system_clock(c0,cr)
    do i = 1, NUM_ENTRIES

      key = ' '
      
      do j = 1, NUM_CHAR
        call random_number(x)
        key(j:j) = char(ichar('a') + floor(x*26))
      end do

      call this%dictionary%insert(key, Foo(i))
      f => this%dictionary%getReference(key)

      call assertEqual(i, f%value)

    end do
    call system_clock(c1)

  end subroutine testStress

end module testAssociativeArray_mod
