.SUFFIXES:

ifndef GISSCLIM_DIR
       GISSCLIM_DIR := $(shell pwd)/../..
endif
ifndef CONFIG_DIR
       CONFIG_DIR := $(GISSCLIM_DIR)/config
endif

#CASE = $(MAKECMDGOALS)

#DEPENDS_ON = shared Ent

FSRCS = GHY.f SNOW_DRV.f SNOW.f

#ifeq ($(MAKECMDGOALS),main)
#  FSRCS += lsm_standalone.f
#endif

include $(CONFIG_DIR)/base.mk
sinclude $(DEPENDFILE)
include $(CONFIG_DIR)/rules.mk

#FFLAGS += $(addprefix -I../,$(DEPENDS_ON)) -DUSE_ENT
#INCS += -I$(shell pwd)/../shared
#ESMF_Interface = ..

#hack: compile for Ent by default
FFLAGS += -DUSE_ENT


