EXTERNAL_CPP = NO
CONFIG_DIR = ../config
SOURCE_DIR = ../model
DECKS_DIR = ../decks

# List of all testable modelE components
INCS := $(foreach ITEM, $(COMPONENTS), $(SOURCE_DIR)/$(ITEM))
COMPLIBS = $(shell perl -e 'print reverse <>;' < $(SOURCE_DIR)/.liblist)
#LIBS := $(foreach ITEM, $(COMPONENTS), $(SOURCE_DIR)/$(ITEM)/lib$(ITEM).a)
LIBS := $(foreach ITEM, $(COMPLIBS), $(SOURCE_DIR)/$(ITEM))

# Determine modelE components to be tested based on library availability

LIBSAVAIL := $(wildcard $(LIBS))
PFUNIT_LFLAGS += $(foreach ITEM, $(LIBSAVAIL), $(ITEM))
PFUNIT_LFLAGS += -L$(PNETCDFHOME)/lib -lpnetcdf
PFUNIT_LFLAGS += -L$(NETCDFHOME)/lib -lnetcdf
PFUNIT_TEST_DIRECTORIES := $(foreach ITEM, $(COMPONENTS), ./$(ITEM))

INCSAVAIL := $(wildcard $(INCS))
EXTRA_FFLAGS += $(foreach ITEM, $(INCSAVAIL), -I$(ITEM)) -I$(SOURCE_DIR)/include

# Defaults

USE_MPI = NO
PFUNIT_USE_MPI = NO

ifeq ($(COMPILER),intel)
  FC = ifort
  EXTRA_FFLAGS += -g -O0 -traceback -fpp -O2 -check bounds
else
  FC = gfortran
  EXTRA_FFLAGS += -g -O0 
endif

EXTRA_FFLAGS += $(CPPFLAGS)

ifeq ($(MPI),YES)
  USE_MPI = YES
  MPI_FC := mpif90
  MPIRUN := mpirun
  PFUNIT_NPES := 10
  PFUNIT_USE_MPI := YES
  PFUNIT_MPI_HOSTFILE=--hostfile my-hostfile
  export PFUNIT_MPI_HOSTFILE
endif

# Targets:

clean: 
	$(RM) *.o *_wrap.F90 *.mod suite_list 

distclean:
	make clean
	$(RM) *.a *.x $(DECKS_DIR)/components

$(PFUNIT_TEST_EXECUTABLE) : $(LIBSAVAIL)

PFUNIT_TEST_PATTERNS = test% Test%
ifeq ($(MPI),YES)
  PFUNIT_TEST_PATTERNS += mpiTest% 
endif

ifdef PFUNIT
  include $(PFUNIT)/include/pFUnit.makefile
endif

# Rules

%.o : %.F90
	$(FC) -c $(EXTRA_FFLAGS) $($(PFUNIT_FFLAGS_VARIABLE)) $<

