.PHONY: clean all

ALLFILES:=$(wildcard *.pf)
# We need to exclude MPI test files whe running serial tests to 
# avoid compilation errors.
MPIFILES:=Test_TimerParallel.pf

ifeq ($(USE_MPI),YES)
  SRCS := $(ALLFILES)
else
  SRCS := $(filter-out $(MPIFILES),$(ALLFILES))
endif
OBJS := $(SRCS:.pf=.o)

LIB = libprofiler.a

$(LIB): $(OBJS)
	ar -r $@ $?
	-@ranlib $@

all: createTestSuites $(OBJS) $(LIB)

createTestSuites:
	rm -f testSuites.inc
	@for f in ${SRCS}; do \
	  echo ADD_TEST_SUITE\($${f%.*}\_mod\_suite\) >> testSuites.inc; \
	done;

testSuites.inc: $(SRCS)

%.F90: %.pf
	python $(PFUNIT)/bin/pFUnitParser.py $<  $@

%.o: %.F90
	$(F90) -c $(F90FLAGS) $<

clean:
	$(RM) *.o *.mod *.a

