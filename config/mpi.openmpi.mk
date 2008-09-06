
ifneq (${MPIDIR},)
FFLAGS += -I${MPIDIR}/include
F90FLAGS += -I${MPIDIR}/include
LIBS += -L${MPIDIR}/lib
endif


#LIBS += -lcprts -limf -lm -lcxa -lunwind -lrt -ldl \
#-lfmpi -lmpi -lstdc++ -threads

LIBS += -lcprts -limf -lm -lcxa -lunwind -lrt -ldl \
-lmpi_f90 -lmpi_f77 -lmpi -lstdc++ -threads


#LIBS +=  -lmpi -lmpigc3 -lmpigc4 -lmpigf -lmpigi -lmpiic4 -lmpiic -lmpiif -threads

