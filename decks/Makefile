.SUFFIXES:
.PHONY: rundeck depend gcm setup clean vclean newstart exe cmp .R htmldoc

# DECKS_REPOSITORY - a directory for permanenet storage of run info 
DECKS_REPOSITORY = /u/cmrun/modelE/decks
# OVERWRITE - can "gmake rundeck" overwrite files already in repository
OVERWRITE = NO

MODEL_E_ROOT = ..
MODEL_DIR = $(MODEL_E_ROOT)/model
AUX_DIR = $(MODEL_E_ROOT)/aux
SCRIPTS_DIR = $(MODEL_E_ROOT)/exec
DECKS_DIR = $(shell pwd)
BIN_DIR = $(DECKS_DIR)/$(RUN)_bin

SETUP = $(SCRIPTS_DIR)/setup_e

RUNSRC = E001

ALIAS_LIST = CMPE001 qc

sinclude $(MODEL_DIR)/Rules.make

help:
	$(MAKE) -C $(MODEL_DIR) help

rundeck:
	@if [ -s $(DECKS_REPOSITORY)/$(RUN).R -a $(OVERWRITE) != YES ];\
	then echo '$(RUN).R already exists in repository'; \
	echo 'use OVERWRITE=YES to overwrite'; \
	exit 1; fi
	$(MAKE) -C $(MODEL_DIR) $@ RUN=$(RUN) DECKS_DIR=$(DECKS_DIR) RUNSRC=$(RUNSRC)
	cp -f $(RUN).R $(DECKS_REPOSITORY)/$(RUN).R

# alias for backwards compatibility
mkrun: rundeck     

gcm $(BIN_DIR)/$(RUN).exe:
	$(MAKE) -C $(MODEL_DIR) $@ RUN=$(RUN) DECKS_DIR=$(DECKS_DIR)

setup: gcm
	@echo '------- Loooks like it was compiled OK -----------------------'
	@echo '------- Saving Rundeck and other info to global repository ---'
	@cp -f $(RUN).R $(DECKS_REPOSITORY)/$(RUN).R
	@echo '### Information below describes your run. Do not delete! ###'\
	>> $(DECKS_REPOSITORY)/$(RUN).R
	@date >> $(DECKS_REPOSITORY)/$(RUN).R
	@$(F90) -version >> $(DECKS_REPOSITORY)/$(RUN).R 2>&1
	@if [ -f $(MODEL_DIR)/CVS/Tag ]; then \
	echo "CVS Repository Tag: `cat $(MODEL_DIR)/CVS/Tag`" \
	>> $(DECKS_REPOSITORY)/$(RUN).R; else \
	echo 'CVS Repository: MAIN Branch' \
	>> $(DECKS_REPOSITORY)/$(RUN).R; fi
	@echo '--------- Starting setup for $(RUN) --------------------------'
	@$(SETUP) $(RUN)

aux:
	$(MAKE) -C $(AUX_DIR) RUN=$(RUN) DECKS_DIR=$(DECKS_DIR)

auxqflux:
	$(MAKE) qflux -C $(AUX_DIR) RUN=$(RUN) DECKS_DIR=$(DECKS_DIR)

auxdeep:
	$(MAKE) odeep -C $(AUX_DIR) RUN=$(RUN) DECKS_DIR=$(DECKS_DIR)

clean vclean:
	$(MAKE) -C $(MODEL_DIR) $@
	$(MAKE) -C $(AUX_DIR) $@

exe:  $(BIN_DIR)/$(RUN).exe
	@if [ ! -s /u/cmrun/$(RUN)/I ] ; then \
	echo 'directory /u/cmrun/$(RUN) or input files inside it are missing';\
	echo 'you should make "setup" first'; \
	exit 1; fi
	mv -f $(BIN_DIR)/$(RUN).exe /u/cmrun/$(RUN)/
	@echo 'adding corresponding info to rundeck repository'
	@perl -e "while(<>){ if(/^### Information/){exit;} print; }" \
	<$(DECKS_REPOSITORY)/$(RUN).R  >$(DECKS_REPOSITORY)/$(RUN).R.tmp
	@echo '### Executable was recompiled ' >> $(DECKS_REPOSITORY)/$(RUN).R
	@date >> $(DECKS_REPOSITORY)/$(RUN).R
	@$(F90) -version >> $(DECKS_REPOSITORY)/$(RUN).R 2>&1
	@if [ -f $(MODEL_DIR)/CVS/Tag ]; then \
	echo "CVS Repository Tag: `cat $(MODEL_DIR)/CVS/Tag`" \
	>> $(DECKS_REPOSITORY)/$(RUN).R; else \
	echo 'CVS Repository: MAIN Branch' \
	>> $(DECKS_REPOSITORY)/$(RUN).R; fi
	@echo 'diff old.R new.R :' >> $(DECKS_REPOSITORY)/$(RUN).R
	@if ! diff $(DECKS_REPOSITORY)/$(RUN).R.tmp $(RUN).R \
	>> $(DECKS_REPOSITORY)/$(RUN).R; then echo 'rundeck changed'; fi
	@rm -f $(DECKS_REPOSITORY)/$(RUN).R.tmp
	@echo done
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

htmldoc:
	$(MAKE) -C $(MODEL_DIR) $@ RUN=$(RUN) DECKS_DIR=$(DECKS_DIR)

