
F90 = g95
FMAKEDEP = $(SCRIPTS_DIR)/sfmakedepend
CPPFLAGS += -DCOMPILER_G95
FFLAGS = -fno-second-underscore -O0 -fendian=big
F90FLAGS = -fno-second-underscore -O0 -ffree-form -fendian=big
LFLAGS =
# uncomment next two lines for extensive debugging
# the following switch adds extra debugging
ifeq ($(COMPILE_WITH_TRAPS),YES)
FFLAGS += -trap=INVALID,DIVBYZERO,OVERFLOW -B111 -C
LFLAGS += -lefence
endif
