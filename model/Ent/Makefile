.SUFFIXES:

ifndef GISSCLIM_DIR
       GISSCLIM_DIR := $(shell pwd)/../..
endif

#DEPENDS_ON = shared

#CASE = $(MAKECMDGOALS)
# online/offline source code
ifeq ($(ONLINE),YES)
FSRCS = 
F90SRCS = 
else
FSRCS = ent_prog.f gcm_utils.f gcm_system.f \
ent_forcings.f ent_prescribed_drv.f
F90SRCS = 
endif

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



include $(GISSCLIM_DIR)/arch/base.mk
#dependencies (should be done later automatically):
sinclude $(DEPENDFILE)
include $(GISSCLIM_DIR)/arch/Rules.make

#FFLAGS += -I$(shell pwd)/../shared
#INCS += -I$(shell pwd)/../shared
#ESMF_Interface = ..

ifeq ($(PFT_MODEL),ENT)
CPPFLAGS += -DPFT_MODEL_ENT
endif


# = $(patsubst %, SRC_%, $(MAKECMDGOALS))

#$(CASE): Makefile $(FSRC)
#	echo $^


#deplist:
#	echo $(THIS)_dir: $(DEPENDS_ON:=_dir)
