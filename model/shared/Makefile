# makefile for the component shared
.SUFFIXES:

SUPPORTED_OPTIONS = ESMF

FSRCS =
F90SRCS  = Constants_mod.F90 StringUtilities_mod.F90 CubicEquation_mod.F90  ArrayBundle_mod.F90
F90SRCS += GenericType_mod.F90 KeyValuePair_mod.F90 Parser_mod.F90 Dictionary_mod.F90
F90SRCS += Geometry_mod.F90 stop_model.F90 FileManager_mod.F90
F90SRCS += RunTimeControls_mod.F90 SystemTimers_mod.F90
F90SRCS += Precision_mod.F90 Utilities.F90 Random_mod.F90 GetTime_mod.F90 System.F90
F90SRCS += SpecialFunctions.F90
F90SRCS += TimeConstants.F90 BaseTime.F90 TimeAndCalendar.F90 Calendar.F90 Month.F90 Time.F90
F90SRCS += ModelClock.F90 JulianCalendar_mod.F90

ifndef GISSCLIM_DIR
       GISSCLIM_DIR := $(shell pwd)/../..
endif
ifndef CONFIG_DIR
       CONFIG_DIR := $(GISSCLIM_DIR)/config
endif

include $(CONFIG_DIR)/base.mk
include $(CONFIG_DIR)/rules.mk

ifeq ($(USE_PFUNIT),YES)
$(OBJS): F90FLAGS+=-DUSE_PFUNIT -I$(PFUNIT)/mod
endif

