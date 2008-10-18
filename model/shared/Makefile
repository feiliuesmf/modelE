# makefile for the component shared
.SUFFIXES:

SUPPORTED_OPTIONS = ESMF

FSRCS = CONST.f  PARAM.f  PARSER.f  SYSTEM.f  TRIDIAG.f  UTILDBL.f


ifndef GISSCLIM_DIR
       GISSCLIM_DIR := $(shell pwd)/../..
endif
ifndef CONFIG_DIR
       CONFIG_DIR := $(GISSCLIM_DIR)/config
endif

include $(CONFIG_DIR)/base.mk
sinclude $(DEPENDFILE)
include $(CONFIG_DIR)/rules.mk

CPPFLAGS += -I$(GISSCLIM_DIR)/include
