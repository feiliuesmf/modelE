module TestHashMapIterator_mod
  use pfunit
  use Foo_mod
  use FooHashMap_mod
  implicit none

  type Fixture
    type (FooHashMap), pointer :: dictionary
  end type Fixture

contains

  subroutine setUp(this)
    type (Fixture) :: this
    this%dictionary => newFooHashMap()
    call this%dictionary%insert('A', Foo(1))
    call this%dictionary%insert('B', Foo(2))
    call this%dictionary%insert('C', Foo(3))
  end subroutine setUp

  subroutine tearDown(this)
    type (Fixture) :: this
    call clean(this%dictionary)
    deallocate(this%dictionary)
  end subroutine tearDown
  
  ! hasNext() returns false for empty associative array
  subroutine testIteratorHasNext_empty(this)
    type (Fixture) :: this
    type (FooHashMap), pointer :: emptyDictionary
    type (FooHashMapIterator) :: iterator

    emptyDictionary => newFooHashMap()
    iterator = emptyDictionary%begin()

    call assertFalse(iterator%hasNext())

    call clean(iterator)
    call clean(emptyDictionary)
    
  end subroutine testIteratorHasNext_empty

  ! iter /= end() returns false for empty array
  subroutine testIteratorNotEqual_empty(this)
    type (Fixture) :: this
    type (FooHashMap), pointer :: emptyDictionary
    type (FooHashMapIterator) :: iterator

    emptyDictionary => newFooHashMap()
    iterator = emptyDictionary%begin()

    call assertFalse(iterator /= emptyDictionary%last())

    call clean(iterator)
    call clean(emptyDictionary)
    
  end subroutine testIteratorNotEqual_empty

  ! hasNext() returns true for non-empty associative array
  subroutine testIteratorHasNext(this)
    type (Fixture) :: this
    type (FooHashMapIterator) :: iterator

    iterator = this%dictionary%begin()
    call assertTrue(iterator%hasNext())

    call clean(iterator)
    
  end subroutine testIteratorHasNext

  ! hasNext() returns false for non-empty associative array after
  ! multiple calls to next()
  subroutine testIteratorHasNext_multi(this)
    type (Fixture) :: this
    type (FooHashMapIterator) :: iterator

    iterator = this%dictionary%begin()
    call iterator%next() ! 2nd element of 3
    call assertTrue(iterator%hasNext())
    call iterator%next() ! 3rd element of 3
    call assertFalse(iterator%hasNext())

    call clean(iterator)
    
  end subroutine testIteratorHasNext_multi

  ! Use iterator in a loop context.
  subroutine testLoop(this)
    type (Fixture) :: this
    type (FooHashMapIterator) :: iterator
    type (FooHashMapIterator) :: iterEnd

    integer :: counter

    counter = 0
    iterator = this%dictionary%begin()
    iterEnd  = this%dictionary%last()
    do while (iterator /= iterEnd)
      call iterator%next()
      counter = counter + 1
    end do

    call assertEqual(this%dictionary%size(), counter)

    call clean(iterator)
    
  end subroutine testLoop

  subroutine testKeyValue(this)
    type (Fixture) :: this
    type (FooHashMapIterator) :: iterator
    
    logical :: foundA
    logical :: foundB
    logical :: foundC

    class (Foo), pointer :: valueA
    class (Foo), pointer :: valueB
    class (Foo), pointer :: valueC

    iterator = this%dictionary%begin()

    foundA = .false.
    foundB = .false.
    foundC = .false.

    do while (iterator .ne. this%dictionary%last())
      select case (iterator%key())
      case ('a')
        foundA = .true.
        valueA => iterator%value()
      case ('b')
        foundB = .true.
        valueB => iterator%value()
      case ('c')
        foundC = .true.
        valueC => iterator%value()
      end select
      call iterator%next()
    end do

    call assertTrue( all([ foundA, foundB, foundC ]) )
    call assertEqual([1,2,3], [ valueA%value,valueB%value,valueC%value ])

    call clean(iterator)
    nullify(valueA)
    nullify(valueB)
    nullify(valueC)

  end subroutine testKeyValue

end module TestHashMapIterator_mod
