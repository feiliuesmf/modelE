#
# this Makefile compiles programs in "aux" directory
#

.SUFFIXES:

.PHONY: clean vclean

all: qc cmp conv_rsf

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

#
#  targets
#

cmp:    CMPE001.f 
	$(F90) $(FFLAGS) -o $(BIN_DIR)/CMPE001 CMPE001.f

conv_rsf: conv_rsf.o $(MODEL_DIR)/gcmlib.a
	$(F90) $(FFLAGS) -o $(BIN_DIR)/conv_rsf $^

conv_rsf.o: $(MODEL_DIR)/gcmlib.a

qc: qc.o $(MODEL_DIR)/gcmlib.a
	$(F90) $(FFLAGS) -o $(BIN_DIR)/qc $^

qc.o: $(MODEL_DIR)/gcmlib.a

$(MODEL_DIR)/gcmlib.a: FORCE
	$(MAKE) -C $(MODEL_DIR) gcmlib RUN=$(RUN)

clean:
	rm -f *.LST *~ *.ERR *.CHK *.L

vclean:
	rm -f *.o *.mod *.vo *.inc .depend* *.sig

