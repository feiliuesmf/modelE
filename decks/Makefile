.SUFFIXES:
.PHONY: rundeck depend gcm setup clean vclean newstart exe cmp .R 

MODEL_E_ROOT = ..
MODEL_DIR = $(MODEL_E_ROOT)/model
AUX_DIR = $(MODEL_E_ROOT)/aux
SCRIPTS_DIR = $(MODEL_E_ROOT)/exec
DECKS_DIR = $(shell pwd)
BIN_DIR = $(DECKS_DIR)/$(RUN)_bin

SETUP = $(SCRIPTS_DIR)/setup_e

RUNSRC = E001

ALIAS_LIST = CMPE001 qc

help:
	$(MAKE) -C $(MODEL_DIR) help

rundeck:
	$(MAKE) -C $(MODEL_DIR) $@ RUN=$(RUN) DECKS_DIR=$(DECKS_DIR) RUNSRC=$(RUNSRC)

# alias for backwards compatibility
mkrun: rundeck     

gcm $(BIN_DIR)/$(RUN).exe:
	$(MAKE) -C $(MODEL_DIR) $@ RUN=$(RUN) DECKS_DIR=$(DECKS_DIR)

setup: gcm
	@echo '--------- Starting setup for $(RUN) --------------------------'
	@$(SETUP) $(RUN)

aux:
	$(MAKE) -C $(AUX_DIR) RUN=$(RUN) DECKS_DIR=$(DECKS_DIR)

auxqflux:
	$(MAKE) qflux -C $(AUX_DIR) RUN=$(RUN) DECKS_DIR=$(DECKS_DIR)

clean vclean:
	$(MAKE) -C $(MODEL_DIR) $@
	$(MAKE) -C $(AUX_DIR) $@

exe:  $(BIN_DIR)/$(RUN).exe
	@if [ ! -s /u/cmrun/$(RUN)/I ] ; then \
	echo 'directory /u/cmrun/$(RUN) or input files inside it are missing';\
	echo 'you should make "setup" first'; \
	exit 1; fi
	mv -f $(BIN_DIR)/$(RUN).exe /u/cmrun/$(RUN)/

# alias for backwards compatibility
mkexe:  exe

newstart:
	@if [ "$(RUN)" = "" ]; then \
	echo 'You should specify run name on a command line (i.e. RUN=...)'; \
	exit 1; fi
	rm -i /u/cmrun/$(RUN)/*

alias: aux
	for i in $(ALIAS_LIST); do ln -sf $(BIN_DIR)/$$i .; done

unalias:
	for i in $(ALIAS_LIST); do rm -f $$i; done

