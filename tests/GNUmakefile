COMPILER=intel
EXTERNAL_CPP=NO
CONFIG_DIR=../config
#sinclude $(CONFIG_DIR)/rules.mk

PFUNIT_TEST_DIRECTORIES =./profiler ./shared ./tracers

SOURCE_DIR = ../model

PFUNIT_LFLAGS+= -L$(SOURCE_DIR)/profiler -lprofiler
PFUNIT_LFLAGS+= -L$(SOURCE_DIR)/shared -lshared
PFUNIT_LFLAGS+= -L$(SOURCE_DIR)/tracers -ltracers

EXTRA_FFLAGS += -$I$(SOURCE_DIR)/profiler
EXTRA_FFLAGS += -$I$(SOURCE_DIR)/shared
EXTRA_FFLAGS += -$I$(SOURCE_DIR)/tracers

MPI_FC=mpif90
CSTRING = $(shell mpif90 --version)
COMPILER = $(findstring ifort, $(CSTRING))

ifeq ($(COMPILER),ifort)
  FC =ifort
  EXTRA_FFLAGS += -g -O0 -traceback -fpp -O2
else
  FC =gfortran
  EXTRA_FFLAGS += -g -O0 
endif

EXTRA_FFLAGS += $(CPPFLAGS)

MPIRUN=mpirun
PFUNIT_USE_MPI=YES
PFUNIT_NPES=10

clean: 
	make pfunitdistclean
	$(RM) *.o *_wrap.F90 *.mod suite_list

distclean:
	make clean
	$(RM) *.a *.x 


ifdef PFUNIT
  include $(PFUNIT)/include/pFUnit.makefile
endif

$(PFUNIT_TEST_EXECUTABLE) : ../model/profiler/libprofiler.a ../model/shared/libshared.a ../model/tracers/libtracers.a

%.o : %.F90
	$(FC) -c $(EXTRA_FFLAGS) $($(PFUNIT_FFLAGS_VARIABLE)) $<

