COMPILER = intel
EXTERNAL_CPP = NO
CONFIG_DIR = ../config
SOURCE_DIR = ../model
DECKS_DIR = ../decks

# List of all testable modelE components

include $(DECKS_DIR)/components

INCS := $(foreach ITEM, $(COMPONENTS), $(SOURCE_DIR)/$(ITEM))
LIBS := $(foreach ITEM, $(COMPONENTS), $(SOURCE_DIR)/$(ITEM)/lib$(ITEM).a)

# Determine modelE components to be tested based on library availability

LIBSAVAIL := $(wildcard $(LIBS))
PFUNIT_LFLAGS += $(foreach ITEM, $(LIBSAVAIL), $(ITEM))
PFUNIT_TEST_DIRECTORIES := $(foreach ITEM, $(COMPONENTS), ./$(ITEM))

INCSAVAIL = $(wildcard $(INCS))
EXTRA_FFLAGS += $(foreach ITEM, $(INCSAVAIL), -I$(ITEM))

# Defaults

MPI_FC=mpif90
CSTRING = $(shell mpif90 --version)
COMPILER = $(findstring ifort, $(CSTRING))
MPIRUN=mpirun
PFUNIT_NPES=10
PFUNIT_USE_MPI=YES

ifeq ($(COMPILER),ifort)
  FC =ifort
  EXTRA_FFLAGS += -g -O0 -traceback -fpp -O2
else
  FC =gfortran
  EXTRA_FFLAGS += -g -O0 
endif

EXTRA_FFLAGS += $(CPPFLAGS)

# Targets:

clean: 
	make pfunitdistclean
	$(RM) *.o *_wrap.F90 *.mod suite_list

distclean:
	make clean
	$(RM) *.a *.x 

ifdef PFUNIT
  include $(PFUNIT)/include/pFUnit.makefile
endif

$(PFUNIT_TEST_EXECUTABLE) : $(LIBSAVAIL)

# Rules

%.o : %.F90
	$(FC) -c $(EXTRA_FFLAGS) $($(PFUNIT_FFLAGS_VARIABLE)) $<

