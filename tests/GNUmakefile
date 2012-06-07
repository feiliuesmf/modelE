COMPILER=intel
EXTERNAL_CPP=NO
CONFIG_DIR=../config
#sinclude $(CONFIG_DIR)/rules.mk

COMPONENTS := profiler shared tracers MPI_Support

PFUNIT_TEST_DIRECTORIES = $(foreach ITEM,$(COMPONENTS), ./$(ITEM))

SOURCE_DIR := ../model

PFUNIT_LFLAGS += $(foreach ITEM,$(COMPONENTS), -L$(SOURCE_DIR)/$(ITEM) -l$(ITEM))

EXTRA_FFLAGS += $(foreach ITEM,$(COMPONENTS), -I$(SOURCE_DIR)/$(ITEM))

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

EXECS = $(foreach ITEM,$(COMPONENTS), $(SOURCE_DIR)/$(ITEM)/lib$(ITEM).a)

$(PFUNIT_TEST_EXECUTABLE) : $(EXECS)

%.o : %.F90
	$(FC) -c $(EXTRA_FFLAGS) $($(PFUNIT_FFLAGS_VARIABLE)) $<

