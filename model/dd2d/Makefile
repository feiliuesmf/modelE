# makefile for the component dd2d
.SUFFIXES:

SUPPORTED_OPTIONS = ESMF

FSRCS = dd2d_utils.f pario_nc.f cdl_mod.f


ifndef GISSCLIM_DIR
       GISSCLIM_DIR := $(shell pwd)/../..
endif
ifndef CONFIG_DIR
       CONFIG_DIR := $(GISSCLIM_DIR)/config
endif

include $(CONFIG_DIR)/base.mk
include $(CONFIG_DIR)/rules.mk

