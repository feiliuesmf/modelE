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
BIN_DIR = $(DECKS_DIR)/$(RUN).bin

GCMLIB = $(BIN_DIR)/$(RUN).a

# netcdf library
NCLIB	= /usr/local/netcdf-3.4/lib64/libnetcdf.a

# list of targets to be built by default
all: $(BIN_DIR)/qc $(BIN_DIR)/CMPE001 $(BIN_DIR)/conv_rsf
all: $(BIN_DIR)/acc2nc

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

$(BIN_DIR)/qc: qc.o $(GCMLIB)
	$(F90) $(FFLAGS) -o $@ $^

qc.o: $(GCMLIB)

$(GCMLIB): FORCE
	$(MAKE) -C $(MODEL_DIR) $@ RUN=$(RUN)

$(BIN_DIR)/acc2nc: acc2nc.o $(GCMLIB) $(NCLIB)
	$(F90) $(FFLAGS) -o $@ $^

acc2nc.o: $(GCMLIB)

clean:
	rm -f *.LST *~ *.ERR *.CHK *.L

vclean:
	rm -f *.o *.mod *.vo *.inc .depend* *.sig

