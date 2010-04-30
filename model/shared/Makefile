# makefile for the component shared
.SUFFIXES:

SUPPORTED_OPTIONS = ESMF

FSRCS = CONST.f  PARAM.f  PARSER.f  SYSTEM.f  UTILDBL.f cubic_eq.f \
        array_bundle.f string_funcs.f
F90SRCS = Parser_mod.F90 Dictionary_mod.F90 Geometry_mod.F90 \
	JulianCalendar_mod.F90 stop_model.F90



ifndef GISSCLIM_DIR
       GISSCLIM_DIR := $(shell pwd)/../..
endif
ifndef CONFIG_DIR
       CONFIG_DIR := $(GISSCLIM_DIR)/config
endif

include $(CONFIG_DIR)/base.mk
include $(CONFIG_DIR)/rules.mk

