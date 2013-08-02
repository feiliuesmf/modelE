.PHONY: clean tests

SRCS = $(wildcard *.pf)
OBJS = $(SRCS:.pf=.o)

LIB = libtracers.a

$(LIB): $(OBJS)
	ar -r $@ $?
	-@ranlib $@

all: $(OBJS) $(LIB)

testSuites.inc: $(SRCS)

%.F90: %.pf
	python $(PFUNIT)/bin/pFUnitParser.py $<  $@

%.o: %.F90
	$(F90) -c $(F90FLAGS) $<

clean:
	$(RM) *.o *.mod *.a

