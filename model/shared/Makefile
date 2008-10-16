.SUFFIXES:

ifndef GISSCLIM_DIR
       GISSCLIM_DIR := $(shell pwd)/../..
endif
ifndef CONFIG_DIR
       CONFIG_DIR := $(GISSCLIM_DIR)/config
endif

#CASE = $(MAKECMDGOALS)

#FSRCS = $(shell ls *.f)
FSRCS = CONST.f  PARAM.f  PARSER.f  SYSTEM.f  TRIDIAG.f  UTILDBL.f

include $(CONFIG_DIR)/base.mk
#dependencies (should be done later automatically):
sinclude $(DEPENDFILE)
include $(CONFIG_DIR)/rules.mk

CPPFLAGS += -I$(GISSCLIM_DIR)/include
