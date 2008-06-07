.SUFFIXES:

ifndef GISSCLIM_DIR
       GISSCLIM_DIR := $(shell pwd)/../..
endif

#CASE = $(MAKECMDGOALS)

FSRCS = $(shell ls *.f)

include $(GISSCLIM_DIR)/arch/base.mk
include $(GISSCLIM_DIR)/arch/Rules.make

#FFLAGS += -I../shared

# = $(patsubst %, SRC_%, $(MAKECMDGOALS))

#$(CASE): Makefile $(FSRC)
#	echo $^

#dependencies (should be done later automatically):
