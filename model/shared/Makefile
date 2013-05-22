# makefile for the component shared
.SUFFIXES:

SUPPORTED_OPTIONS = ESMF
F90SRCS := $(wildcard *.F90)

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

