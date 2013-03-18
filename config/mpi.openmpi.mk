
# if MPI*DIR are not defined try to get them from MPIDIR
ifneq ($(MPIDIR),)
  ifneq ($(wildcard $(MPIDIR)/include/mpi.h),)
    MPIINCLUDEDIR ?= $(MPIDIR)/include
  else
    ifneq ($(wildcard $(MPIDIR)/include/openmpi/mpi.h),)
      MPIINCLUDEDIR ?= $(MPIDIR)/include/openmpi
    endif
  endif
  ifeq ($(MPIINCLUDEDIR),)
    $(error MPI distribution not found. Check settings in ~/.modelErc)
  endif
  MPILIBDIR ?= $(MPIDIR)/lib
endif

# if MPI*DIR are not yet defined try to get them from a module
MPIINCLUDEDIR ?= $(MPI_INCLUDE)
MPILIBDIR ?= $(MPI_LIB)

ifneq ($(MPIINCLUDEDIR),)
  FFLAGS += -I$(MPIINCLUDEDIR)
  F90FLAGS += -I$(MPIINCLUDEDIR)
  CPPFLAGS += -I$(MPIINCLUDEDIR)
endif
ifneq ($(MPILIBDIR),)
  LIBS += -L$(MPILIBDIR)
endif


# try to work around memory leak
CPPFLAGS += -DMPITYPE_LOOKUP_HACK

LIBS += -lmpi_f77 -lmpi -lmpi_cxx -lstdc++
ifneq ($(shell uname),Darwin)
  LIBS += -lrt
endif

