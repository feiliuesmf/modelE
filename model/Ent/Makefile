#makefile for the component Ent
.SUFFIXES:

SUPPORTED_OPTIONS = ESMF PFT_MODEL PS_MODEL RAD_MODEL

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
  FSRCS += ent_ENTveg.f
  ifeq ($(FLUXNET),YES)
    FSRCS += ent_pfts_ENT_FLUXNET.f
  else
    FSRCS += ent_pfts_ENT.f
  endif
else
    FSRCS += ent_GISSveg.f
  ifeq ($(FLUXNET),YES)
    FSRCS += ent_pfts_FLUXNET.f
  else
    FSRCS += ent_pfts.f
  endif
endif

ifeq ($(PS_MODEL),FBB)
  FSRCS += FBBphotosynthesis.f
  ifeq ($(RAD_MODEL),GORT)
    FSRCS += canopygort.f
  else
    FSRCS += canopyspitters.f
  endif
  ifeq ($(PFT_MODEL),ENT)
    ifeq ($(FLUXNET),YES)
      FSRCS += FBBpfts_ENT_FLUXNET.f
    else
      FSRCS += FBBpfts_ENT.f
    endif
  else
    ifeq ($(FLUXNET),YES)
      FSRCS += FBBpfts_FLUXNET.f
    else
      FSRCS += FBBpfts.f
    endif
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
#sinclude $(DEPENDFILE)
include $(CONFIG_DIR)/rules.mk

ifeq ($(PFT_MODEL),ENT)
CPPFLAGS += -DPFT_MODEL_ENT
endif

ifeq ($(FLUXNET),YES)
CPPFLAGS += -DSOILCARB_SITE
endif

ifeq ($(ENT_STANDALONE_DIAG),YES)
CPPFLAGS += -DENT_STANDALONE_DIAG
endif

