
F90 = ifort
FMAKEDEP = $(SCRIPTS_DIR)/sfmakedepend
CMP_MOD = $(SCRIPTS_DIR)/compare_module_file.pl -compiler INTEL-ifort-9-0-on-LINUX
FFLAGS = -fpp -O2 -ftz         -convert big_endian 
F90FLAGS = -fpp -O2 -ftz        -convert big_endian -free 
LFLAGS = -O2 -ftz
CPPFLAGS += -DCOMPILER_Intel8 -DCONVERT_BIGENDIAN
F90_VERSION = $(shell $(F90) -v 2>&1)
ifeq ($(MP),YES)
FFLAGS += -openmp
F90FLAGS += -openmp
LFLAGS += -openmp
endif
ifeq ($(COMPILE_WITH_TRAPS),YES)
FFLAGS += -CB -fpe0
LFLAGS += -CB -fpe0
F90FLAGS += -CB -fpe0
LFLAGSF += -CB -fpe0
endif
