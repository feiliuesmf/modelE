#makefile for the component Ent
.SUFFIXES:

SUPPORTED_OPTIONS = ESMF PFT_MODEL PS_MODEL

FSRCS = 
F90SRCS = 

FSRCS += ent_prescribed_drv.f # needed for lsm_standalone

FSRCS += ent_mod.f ent.f \
cohorts.f patches.f entcells.f \
canopyradiation.f \
growthallometry.f reproduction.f phenology.f \
disturbance.f soilbgc.f ent_const.f \
ent_types.f  ent_prescribed_updates.f

# different versions of the code
ifeq ($(PFT_MODEL),ENT)
FSRCS += ent_pfts_ENT.f  ent_ENTveg.f
else
FSRCS += ent_pfts.f  ent_GISSveg.f
endif

ifeq ($(PS_MODEL),FBB)
FSRCS += FBBphotosynthesis.f  canopyspitters.f 
      ifeq ($(PFT_MODEL),ENT)
      FSRCS += FBBpfts_ENT.f
      else
      FSRCS += FBBpfts.f
      endif
else
FSRCS += biophysics.f
endif


ifndef GISSCLIM_DIR
       GISSCLIM_DIR := $(shell pwd)/../..
endif
ifndef CONFIG_DIR
       CONFIG_DIR := $(GISSCLIM_DIR)/config
endif

include $(CONFIG_DIR)/base.mk
#dependencies (should be done later automatically):
sinclude $(DEPENDFILE)
include $(CONFIG_DIR)/rules.mk

ifeq ($(PFT_MODEL),ENT)
CPPFLAGS += -DPFT_MODEL_ENT
endif

