# makefile for the component giss_LSM
.SUFFIXES:

SUPPORTED_OPTIONS = ESMF

FSRCS = GHY.f SNOW_DRV.f SNOW.f


ifndef GISSCLIM_DIR
       GISSCLIM_DIR := $(shell pwd)/../..
endif
ifndef CONFIG_DIR
       CONFIG_DIR := $(GISSCLIM_DIR)/config
endif

include $(CONFIG_DIR)/base.mk
sinclude $(DEPENDFILE)
include $(CONFIG_DIR)/rules.mk

#hack: compile for Ent by default
CPPFLAGS += -DUSE_ENT


