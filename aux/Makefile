# 
# this Makefile compiles programs in "aux" directory
#

.SUFFIXES:

.PHONY: clean vclean

MODEL_E_ROOT = ..
MODEL_DIR = $(MODEL_E_ROOT)/model
DECKS_DIR = $(MODEL_E_ROOT)/decks
SCRIPTS_DIR = $(MODEL_E_ROOT)/exec

RUN     =
RUNDECK = $(DECKS_DIR)/$(RUN).R

OBJ_DIR = .
BIN_DIR = $(DECKS_DIR)/$(RUN)_bin

GCMLIB = $(BIN_DIR)/$(RUN).a

# getting options from configuration files
# default modelE configuration file:
include $(SCRIPTS_DIR)/modelErc
# configuration file with user-defined options:
sinclude ~/.modelErc         


# list of targets to be built by default
all: $(BIN_DIR)/qc $(BIN_DIR)/CMPE001 $(BIN_DIR)/conv_rsf $(BIN_DIR)/conv_gic $(BIN_DIR)/rsf2gic $(BIN_DIR)/conv_strat_rsf  $(BIN_DIR)/RMS $(BIN_DIR)/acc2nc $(BIN_DIR)/CMPE002 $(BIN_DIR)/ext_gic
qflux: $(BIN_DIR)/vertflux.exe $(BIN_DIR)/ocnIC.exe $(BIN_DIR)/otspec.E001.exe 
# $(BIN_DIR)/rsfE001.addO.exe obsolete
odeep: $(BIN_DIR)/mkdeep.exe

# include file with implicit rules
include $(MODEL_DIR)/Rules.make

# look for .mod files in MODEL_DIR
ifeq ($(COMPILER),Absoft)
FFLAGS += -p $(MODEL_DIR) -I$(MODEL_DIR)
FFLAGSF += -p $(MODEL_DIR) -I$(MODEL_DIR)
else
FFLAGS += -I$(MODEL_DIR)
FFLAGSF += -I$(MODEL_DIR)
endif

#
#  targets
#

$(BIN_DIR)/CMPE002: CMPE002.F90 $(GCMLIB)
	$(F90) $(FFLAGSF) $(EXTRA_FFLAGS) $(CPPFLAGS) $(LFLAGS) -o $@ $^ $(LIBS)
CMPE001.o: $(GCMLIB)

$(BIN_DIR)/CMPE001: CMPE001.o $(GCMLIB)
	$(F90) $(LFLAGS) -o $@ $^ $(LIBS)

RMS.o: $(GCMLIB)

$(BIN_DIR)/RMS: RMS.o $(GCMLIB)
	$(F90) $(LFLAGS) -o $@ $^ $(LIBS)

conv_rsf.o: $(GCMLIB)

$(BIN_DIR)/conv_rsf: conv_rsf.o $(GCMLIB)
	$(F90) $(LFLAGS) -o $@ $^ $(LIBS)

conv_strat_rsf.o: $(GCMLIB)

$(BIN_DIR)/conv_strat_rsf: conv_strat_rsf.o $(GCMLIB)
	$(F90) $(LFLAGS) -o $@ $^ $(LIBS)

conv_gic.o: $(GCMLIB)

$(BIN_DIR)/conv_gic: conv_gic.o $(GCMLIB)
	$(F90) $(LFLAGS) -o $@ $^ $(LIBS)

rsf2gic.o: $(GCMLIB)

$(BIN_DIR)/rsf2gic: rsf2gic.o $(GCMLIB)
	$(F90) $(LFLAGS) -o $@ $^ $(LIBS)

ext_gic.o: $(GCMLIB)

$(BIN_DIR)/ext_gic: ext_gic.o $(GCMLIB)
	$(F90) $(LFLAGS) -o $@ $^ $(LIBS)

qc.o: $(GCMLIB)

$(BIN_DIR)/qc: qc.o $(GCMLIB)
	$(F90) $(LFLAGS) -o $@ $^ $(LIBS)

vertflux.E001.o: $(GCMLIB) MAPS.o

$(BIN_DIR)/vertflux.exe: vertflux.E001.o MAPS.o $(GCMLIB)
	$(F90) $(LFLAGS) -o $@ $^ $(LIBS)

ocnIC.E001.o: $(GCMLIB) MAPS.o

$(BIN_DIR)/ocnIC.exe: ocnIC.E001.o MAPS.o $(GCMLIB)
	$(F90) $(LFLAGS) -o $@ $^ $(LIBS)

otspec.E001.o: $(GCMLIB) MAPS.o

$(BIN_DIR)/otspec.E001.exe: otspec.E001.o MAPS.o $(GCMLIB)
	$(F90) $(LFLAGS) -o $@ $^ $(LIBS)

#$(BIN_DIR)/rsfE001.addO.exe: rsfE001.addO.o $(GCMLIB)
#	$(F90) $(FFLAGS) $(LFLAGS) -o $@ $^
#
#rsfE001.addO.o: $(GCMLIB)

mkdeep.o: $(GCMLIB)

$(BIN_DIR)/mkdeep.exe: mkdeep.o $(GCMLIB)
	$(F90) $(LFLAGS) -o $@ $^ $(LIBS)

$(GCMLIB): FORCE
	$(MAKE) -C $(MODEL_DIR) $@ RUN=$(RUN)

ifdef NETCDFHOME
acc2nc.o: $(GCMLIB)

$(BIN_DIR)/acc2nc: acc2nc.o $(GCMLIB)
	$(F90) $(LFLAGS) -o $@ $^ $(LIBS)
else
$(BIN_DIR)/acc2nc:
	@echo "NetCDF library is not specified"
	@echo "will skip $(BIN_DIR)/acc2nc"
endif

clean:
	rm -f *.LST *~ *.ERR *.CHK *.L

vclean:
	rm -f *.o *.mod *.vo *.inc .depend* *.sig

