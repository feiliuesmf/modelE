#
# this Makefile compiles programs in "aux" directory
#

.SUFFIXES:

.PHONY: clean vclean

all: qc cmp conv_rsf acc2nc

MODEL_E_ROOT = ..
MODEL_DIR = $(MODEL_E_ROOT)/model
DECKS_DIR = .
SCRIPTDIR = $(MODEL_E_ROOT)/exec

RUN     =
RUNDECK = $(DECKS_DIR)/$(RUN).R

OBJ_DIR = .
BIN_DIR = $(DECKS_DIR)/$(RUN).bin

include $(MODEL_DIR)/Rules.make

# look for .mod files in MODEL_DIR
FFLAGS += -I $(MODEL_DIR)

# netcdf library
NCLIB	= /usr/local/netcdf-3.4/lib64/libnetcdf.a

#
#  targets
#

cmp:    CMPE001.f 
	$(F90) $(FFLAGS) -o CMPE001 CMPE001.f

conv_rsf: conv_rsf.o $(MODEL_DIR)/gcmlib.a
	$(F90) $(FFLAGS) -o $(BIN_DIR)/conv_rsf $^

conv_rsf.o: $(MODEL_DIR)/gcmlib.a

qc: qc.o $(MODEL_DIR)/gcmlib.a
	$(F90) $(FFLAGS) -o $(BIN_DIR)/qc $^

qc.o: $(MODEL_DIR)/gcmlib.a

$(MODEL_DIR)/gcmlib.a: FORCE
	$(MAKE) -C $(MODEL_DIR) gcmlib RUN=$(RUN)

acc2nc: acc2nc.o $(MODEL_DIR)/NCACC.o $(MODEL_DIR)/gcmlib.a $(NCLIB)
	$(F90) $(FFLAGS) -o $(BIN_DIR)/$@ $^

$(MODEL_DIR)/NCACC.o: $(MODEL_DIR)/gcmlib.a # recomp if any changes
	$(MAKE) -C $(MODEL_DIR) NCACC.o

clean:
	rm -f *.LST *~ *.ERR *.CHK *.L

vclean:
	rm -f *.o *.mod *.vo *.inc .depend* *.sig

