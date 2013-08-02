.PHONY: clean all

SRCS = $(wildcard *.pf)
OBJS = $(SRCS:.pf=.o FooAssociativeArray.o FooHashMap.o)

LIB = libshared.a

$(LIB): $(OBJS)
	ar -r $@ $?
	-@ranlib $@

F90FLAGS += -I../../model/include -I../../model/mod

all: FooAssociativeArray.o FooHashMap.o $(OBJS) $(LIB)

testSuites.inc: $(SRCS)

%.F90: %.pf
	python $(PFUNIT)/bin/pFUnitParser.py $<  $@

%.o: %.F90
	$(F90) -c $(F90FLAGS) $<

clean:
	$(RM) *.o *.mod *.a

FooAssociativeArray.o: ../../model/include/AssociativeArrayTemplate.h
FooHashMap.o: ../../model/include/HashMapTemplate.h
Test_AssociativeArray_mod.o : FooAssociativeArray.o
Test_HashMap_mod.o : FooHashMap.o FooAssociativeArray.o
Test_AssociativeArrayIterator_mod.o : FooAssociativeArray.o

