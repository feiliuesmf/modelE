COMPILER=intel
EXTERNAL_CPP=NO
CONFIG_DIR=../config
#include $(CONFIG_DIR)/rules.mk

PFUNIT_TEST_DIRECTORIES =./profiler ./shared

SOURCE_DIR = ../model

PFUNIT_LFLAGS+= -L$(SOURCE_DIR)/profiler -lprofiler -L$(SOURCE_DIR)/shared -lshared
EXTRA_FFLAGS += -$I$(SOURCE_DIR)/profiler -$I$(SOURCE_DIR)/shared -DUSE_PFUNIT
PFUNIT_TEST_PATTERNS=test%
MPI_FC=mpif90
FC=ifort
MPIRUN=mpirun
PFUNIT_USE_MPI=YES
PFUNIT_NPES=10


ifdef PFUNIT
  include $(PFUNIT)/include/pFUnit.makefile
endif

$(PFUNIT_TEST_EXECUTABLE) : ../model/profiler/libprofiler.a ../model/shared/libshared.a

%.o : %.F90
	$(FC) -c $(EXTRA_FFLAGS) $($(PFUNIT_FFLAGS_VARIABLE)) $<

