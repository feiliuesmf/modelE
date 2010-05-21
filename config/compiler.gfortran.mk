
F90 = gfortran
FMAKEDEP = $(SCRIPTS_DIR)/sfmakedepend
CPPFLAGS += -DCOMPILER_G95
FFLAGS = -cpp -fconvert=big-endian -O2
F90FLAGS = -cpp -fconvert=big-endian -O2
LFLAGS =

# machine-specific options
ifeq ($(MACHINE),IRIX64)
FFLAGS += -mabi=64 
F90FLAGS += -mabi=64
LFLAGS += -mabi=64
endif

# uncomment next two lines for extensive debugging
# the following switch adds extra debugging
ifeq ($(COMPILE_WITH_TRAPS),YES)
FFLAGS += -trap=INVALID,DIVBYZERO,OVERFLOW -B111 -C
LFLAGS += -lefence
endif
