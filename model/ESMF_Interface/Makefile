ifneq ($(NEW_MAKEFILE_FORMAT),YES)

.SUFFIXES:

MODEL_E_ROOT = ../..
CONFIG_DIR = $(MODEL_E_ROOT)/config

all:  ESMF_MOD.o ESMF_CUSTOM_MOD.o assert.o

# dependencies
ESMF_MOD.o: 		ESMF_MOD.F90
assert.o:		assert.F90 ESMF_MOD.o  assert.h
ESMF_CUSTOM_MOD.o: 	ESMF_CUSTOM_MOD.F90 assert.o ESMF_MOD.o assert.h

clean:
	rm -f *.o *.mod *.sig *.sig1

include $(CONFIG_DIR)/rules.mk

else
###################################################################

.SUFFIXES:

SUPPORTED_OPTIONS =

F90SRCS = ESMF_MOD.F90 assert.F90 ESMF_CUSTOM_MOD.F90
FSRCS = DOMAIN_DECOMP.f

ifndef GISSCLIM_DIR
       GISSCLIM_DIR := $(shell pwd)/../..
endif
ifndef CONFIG_DIR
       CONFIG_DIR := $(GISSCLIM_DIR)/config
endif

include $(CONFIG_DIR)/base.mk
#sinclude $(DEPENDFILE)
include $(CONFIG_DIR)/rules.mk

no_default:
	echo "no default target!!!"
	echo "GISSCLIM_DIR = $(GISSCLIM_DIR)"
	echo "CONFIG_DIR = $(CONFIG_DIR)"
	exit 1

#####################################################################
endif
