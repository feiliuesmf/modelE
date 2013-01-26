
ifneq (${MPIDIR},)
ifneq ($(wildcard $(MPIDIR)/include/mpi.h),)
  FFLAGS += -I${MPIDIR}/include
  F90FLAGS += -I${MPIDIR}/include
  CPPFLAGS += -I${MPIDIR}/include
else
  ifneq ($(wildcard $(MPIDIR)/include/openmpi/mpi.h),)
    FFLAGS += -I${MPIDIR}/include/openmpi
    F90FLAGS += -I${MPIDIR}/include/openmpi
    CPPFLAGS += -I${MPIDIR}/include/openmpi
  else
    $(error MPI distribution not found. Check settings in ~/.modelErc)
  endif
endif
LIBS += -L${MPIDIR}/lib
endif

# try to work around memory leak
CPPFLAGS += -DMPITYPE_LOOKUP_HACK

#LIBS += -lcprts -limf -lm -lcxa -lunwind -lrt -ldl \
#-lfmpi -lmpi -lstdc++ -threads

##LIBS += -lcprts -limf -lm -lcxa -lunwind -lrt -ldl \
##-lmpi_f90 -lmpi_f77 -lmpi -lstdc++ -threads

LIBS += -lmpi_f77 -lmpi -lmpi_cxx -lstdc++
ifneq ($(shell uname),Darwin)
  LIBS += -lrt
endif

#LIBS +=  -lmpi -lmpigc3 -lmpigc4 -lmpigf -lmpigi -lmpiic4 -lmpiic -lmpiif -threads

