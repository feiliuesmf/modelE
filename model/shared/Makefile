# makefile for the component shared
.SUFFIXES:

SUPPORTED_OPTIONS = ESMF

FSRCS = SYSTEM.f
F90SRCS  = Constants_mod.F90 StringUtilities_mod.F90 CubicEquation_mod.F90  ArrayBundle_mod.F90
F90SRCS += GenericType_mod.F90 KeyValuePair_mod.F90 Parser_mod.F90 Dictionary_mod.F90
F90SRCS += Geometry_mod.F90 JulianCalendar_mod.F90 stop_model.F90 FileManager_mod.F90
F90SRCS += RunTimeControls_mod.F90 SystemTimers_mod.F90
F90SRCS += Precision_mod.F90 Utilities.F90 Random_mod.F90 GetTime_mod.F90


ifndef GISSCLIM_DIR
       GISSCLIM_DIR := $(shell pwd)/../..
endif
ifndef CONFIG_DIR
       CONFIG_DIR := $(GISSCLIM_DIR)/config
endif

include $(CONFIG_DIR)/base.mk
include $(CONFIG_DIR)/rules.mk

