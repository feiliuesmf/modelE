.SUFFIXES:

ifndef GISSCLIM_DIR
       GISSCLIM_DIR := $(shell pwd)/../..
endif

#CASE = $(MAKECMDGOALS)

#DEPENDS_ON = shared Ent

FSRCS = GHY.f SNOW_DRV.f SNOW.f

#ifeq ($(MAKECMDGOALS),main)
#  FSRCS += lsm_standalone.f
#endif

include $(GISSCLIM_DIR)/arch/base.mk
sinclude $(DEPENDFILE)
include $(GISSCLIM_DIR)/arch/Rules.make

#FFLAGS += $(addprefix -I../,$(DEPENDS_ON)) -DUSE_ENT
#INCS += -I$(shell pwd)/../shared
#ESMF_Interface = ..

#hack: compile for Ent by default
FFLAGS += -DUSE_ENT

# = $(patsubst %, SRC_%, $(MAKECMDGOALS))

#$(CASE): Makefile $(FSRC)
#	echo $^

#dependencies (should be done later automatically):
#sinclude $(DEPENDFILE)


#LIBLIST += $(addsuffix /lib.a,$(addprefix ../,$(DEPENDS_ON)))

#main: $(OBJS)
#	$(F90) $(LFLAGS) $(EXTRA_LFLAGS) $(OBJS) $(F90OBJS) $(ESMF_OBJS) \
#	  $(LIBLIST) $(LIBS) -o main  $(LINK_OUTPUT)


