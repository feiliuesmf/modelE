# 
# this Makefile compiles programs in "aux" directory
#

.SUFFIXES:

.PHONY: clean vclean

MODEL_E_ROOT = ..
MODEL_DIR = $(MODEL_E_ROOT)/model
DECKS_DIR = .
SCRIPTDIR = $(MODEL_E_ROOT)/exec

RUN     =
RUNDECK = $(DECKS_DIR)/$(RUN).R

OBJ_DIR = .
BIN_DIR = $(DECKS_DIR)/$(RUN)_bin

GCMLIB = $(BIN_DIR)/$(RUN).a

# netcdf library
NCLIB = /usr/local/netcdf-3.4/lib64/libnetcdf.a

# list of targets to be built by default
all: $(BIN_DIR)/qc $(BIN_DIR)/CMPE001 $(BIN_DIR)/acc2nc $(BIN_DIR)/conv_rsf $(BIN_DIR)/conv_gic $(BIN_DIR)/conv_strat_rsf
qflux: $(BIN_DIR)/vertflux.exe $(BIN_DIR)/dec31.exe $(BIN_DIR)/otspec.E001.exe $(BIN_DIR)/rsfE001.addO.exe

# include file with implicit rules
include $(MODEL_DIR)/Rules.make

# look for .mod files in MODEL_DIR
FFLAGS += -I $(MODEL_DIR)


#
#  targets
#

$(BIN_DIR)/CMPE001:    CMPE001.f
	$(F90) $(FFLAGS) -o $@ $^

$(BIN_DIR)/conv_rsf: conv_rsf.o $(GCMLIB)
	$(F90) $(FFLAGS) -o $@ $^

conv_rsf.o: $(GCMLIB)

$(BIN_DIR)/conv_strat_rsf: conv_strat_rsf.o $(GCMLIB)
	$(F90) $(FFLAGS) -o $@ $^

conv_strat_rsf.o: $(GCMLIB)

$(BIN_DIR)/conv_gic: conv_gic.o $(GCMLIB)
	$(F90) $(FFLAGS) -o $@ $^

conv_gic.o: $(GCMLIB)

$(BIN_DIR)/qc: qc.o $(GCMLIB)
	$(F90) $(FFLAGS) -o $@ $^

qc.o: $(GCMLIB)

$(BIN_DIR)/vertflux.exe: vertflux.E001.o MAPS.o $(GCMLIB)
	$(F90) $(FFLAGS) -o $@ $^

vertflux.E001.o: $(GCMLIB)

$(BIN_DIR)/dec31.exe: dec31.E001.o MAPS.o $(GCMLIB)
	$(F90) $(FFLAGS) -o $@ $^

dec31.E001.o: $(GCMLIB)

$(BIN_DIR)/otspec.E001.exe: otspec.E001.o MAPS.o $(GCMLIB)
	$(F90) $(FFLAGS) -o $@ $^

otspec.E001.o: $(GCMLIB)

$(BIN_DIR)/rsfE001.addO.exe: rsfE001.addO.o $(GCMLIB)
	$(F90) $(FFLAGS) -o $@ $^

rsfE001.addO.o: $(GCMLIB)

$(GCMLIB): FORCE
	$(MAKE) -C $(MODEL_DIR) $@ RUN=$(RUN)

$(BIN_DIR)/acc2nc: acc2nc.o $(GCMLIB) $(NCLIB)
	$(F90) $(FFLAGS) -o $@ $^

acc2nc.o: $(GCMLIB)

clean:
	rm -f *.LST *~ *.ERR *.CHK *.L

vclean:
	rm -f *.o *.mod *.vo *.inc .depend* *.sig

