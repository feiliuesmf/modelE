.SUFFIXES:
.PHONY: mkrun depend gcm setup clean vclean newstart mkexe cmp .R 

MODEL_E_ROOT = ..
MODEL_DIR = $(MODEL_E_ROOT)/model
AUX_DIR = $(MODEL_E_ROOT)/aux
SCRIPTS_DIR = $(MODEL_E_ROOT)/exec
DECKS_DIR = $(shell pwd)
BIN_DIR = $(DECKS_DIR)/$(RUN).bin

SETUP = $(SCRIPTS_DIR)/setup_e

RUNSRC  = E001

help:
	$(MAKE) -C $(MODEL_DIR) help

mkrun:
	$(MAKE) -C $(MODEL_DIR) $@ RUN=$(RUN) DECKS_DIR=$(DECKS_DIR) RUNSRC=$(RUNSRC)

gcm $(BIN_DIR)/$(RUN).exe:
	$(MAKE) -C $(MODEL_DIR) $@ RUN=$(RUN) DECKS_DIR=$(DECKS_DIR)

setup: gcm
	@echo '--------- Starting setup for $(RUN) --------------------------'
	@$(SETUP) $(RUN)

aux:
	$(MAKE) -C $(AUX_DIR) RUN=$(RUN) DECKS_DIR=$(DECKS_DIR)

clean vclean:
	$(MAKE) -C $(MODEL_DIR) $@
	$(MAKE) -C $(AUX_DIR) $@

mkexe:  $(BIN_DIR)/$(RUN).exe
	@if [ ! -s /u/cmrun/$(RUN)/I ] ; then \
	echo 'directory /u/cmrun/$(RUN) or input files inside it are missing';\
	echo 'you should make "setup" first'; \
	exit 1; fi
	cp -f $(BIN_DIR)/$(RUN).exe /u/cmrun/$(RUN)/

newstart:
	@if [ "$(RUN)" = "" ]; then \
	echo 'You should specify run name on a command line (i.e. RUN=...)'; \
	exit 1; fi
	rm -i /u/cmrun/$(RUN)/*

