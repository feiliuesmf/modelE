# Makefile for Model E Atmosphere GCM B399+ 
.SUFFIXES:
.SUFFIXES: .o .U .S .GCM .OCN .f .F
.SUFFIXES:

.PHONY: rundeck depend gcm setup clean vclean newstart exe cmp .R htmldoc

#
# Print Usage if no targets were given
#
print_usage help all:
	@echo
	@echo '               GCM Model E Makefile'
	@echo 
	@echo 'This Makefile works with GNU make ("gmake" on most systems).'
	@echo 'It WILL NOT work with standard "make" which comes with the '
	@echo 'system (unless you are on Linux). If you see error messages'
	@echo 'above or below this message, that means that you are using'
	@echo 'wrong (not GNU) "make" command.'
	@echo 'This Makefile automatically recognizes and supports three'
	@echo 'types of systems: SGI, IBM workstation and Linux PC. Though'
	@echo '"setup" currently works only on SGIs.'
	@echo
	@echo The following is a list of supported targets:
	@echo
	@echo ' rundeck - create new rundeck'
	@echo ' depend  - create dependencies for specified rundeck'
	@echo ' gcm     - compile object files and build executable for specified rundeck'
	@echo ' aux     - compile standard auxiliary programs'
	@echo ' auxqflux- compile auxiliary programs for computing qflux'
	@echo ' setup   - do setup for specified rundeck'
	@echo ' clean   - remove listings and error files'
	@echo ' vclean  - remove object files, .mod files and dependencies'
	@echo ' newstart- remove all files in the run directory'
	@echo ' exe     - compile gcm and put executable into RUN directory'
	@echo
	@echo '          Recommended usage:'
	@echo 'The targets are made in order: "rundeck", "depend", "gcm", "setup",'
	@echo 'but since Makefile automatically creates necessary targets, you'
	@echo 'do not have to make all of them manually.'
	@echo 'You should always start with "rundeck", it will NOT be created automatically'
	@echo 'as it was in previous version. You may (and should) always skip "depend".'
	@echo 'So if you want to compile gcm and set up a run with a name "my_run"'
	@echo 'typically you will have to execute the commands:'
	@echo
	@echo '       gmake rundeck RUN=my_run'
	@echo '       gmake gcm RUN=my_run'
	@echo '       gmake setup RUN=my_run'
	@echo
	@echo 'You can actually skip "gmake gcm..." either and gcm will be'
	@echo 'compiled automatically, but you may not want to do it in most'
	@echo 'cases, since you want to look at the warning massages before'
	@echo 'doing "setup".'
	@echo


#
# Directory structure and file names
#

MODEL_E_ROOT = ..
MODEL_DIR = .
DECKS_DIR = .
SCRIPTDIR = $(MODEL_E_ROOT)/exec

RUNSRC  = E001
RUN     =
RUNDECKSRC = $(MODEL_DIR)/$(RUNSRC).R
RUNDECK = $(DECKS_DIR)/$(RUN).R

OBJ_DIR = .
BIN_DIR = $(DECKS_DIR)/$(RUN)_bin
DOC_DIR = $(DECKS_DIR)/$(RUN)_doc

DEPENDFILE = $(OBJ_DIR)/.depend.$(RUN)

RUN_H = rundeck_opts.h
RUN_H_TMP = tmp_cpp.h

#
# Lists of objects and dependencies
# 

#ifneq ($(RUN).R,.R)
ifeq ($(shell [ -s $(RUNDECK) ] && echo yes),yes)
OBJ_LIST = $(shell [ -s $(RUNDECK) ] && perl \
	-e 'while(<> !~ /^Object +modules/ ){};' \
	-e 'while(($$_=<>) !~ /(^Data +input)|(^ *$$)/)' \
	-e '{chop; s/ *(!.*)*$$//; s/~/$(MACHINE)/g;print " $$_";}' \
	< $(RUNDECK))

sinclude $(DEPENDFILE)
FSRCS = $(OBJ_LIST:=.f)
OBJS = $(OBJ_LIST:=.o)

MODS = $(sort $(shell perl \
	-e 'while(<>){ tr/a-z/A-Z/;'\
	-e 'if(/^\s*module\s+(\w+)\s*$$/i) {print "$$1.mod ";} };' \
	$(FSRCS)))
endif


#
#include rules here
#

sinclude $(MODEL_DIR)/Rules.make



#########             TARGETS               ########

# Make a current rundeck
rundeck:
	@if [ ! -s $(RUNDECKSRC) ]; then \
	echo 'Source rundeck $(RUNDECKSRC) does not exist'; exit 1; fi
	@if [ "$(RUN)" = "" ]; then \
	echo 'You should specify run name on a command line (i.e. RUN=...)'; \
	exit 1; fi
	@echo '---- Regenerating Rundeck $(RUNDECK) from default ($(RUNDECKSRC)) ----'
	@if [ $(RUNDECK) != $(RUNDECKSRC) ]; then  \
	sed 's/^'$(RUNSRC)'\([:\. ]\)/'$(RUN)'\1/g' $(RUNDECKSRC) > $(RUNDECK); \
	fi;
	@echo '---------         Rundeck  created                   ---------' 

depend $(DEPENDFILE): $(RUNDECK) $(FSRCS)
	@echo '-------- Dependecies are out of date or do not exist ---------'
	@echo '---------         regenerating them                  ---------'
	[ -f $(RUN_H) ] || touch $(RUN_H)
	$(FMAKEDEP) -f $(DEPENDFILE) $^
	@echo '---------         done with dependencies             ---------'

$(RUN_H): $(RUNDECK) FORCE
	-@rm -f $(RUN_H_TMP)
	@perl -e 'while(($$_=<>) && ! /^ *Preprocessor +Options/i ){};' \
	-e 'while(($$_=<>) && ! /^ *End +Preprocessor +Options/i )' \
	-e '{chop; s/ *(!.*)*$$//; $$_ && print "$$_\n";}' \
	< $(RUNDECK) > $(RUN_H_TMP)
	@[ -f $(RUN_H) ] && diff $(RUN_H_TMP) $(RUN_H) >/dev/null 2>&1  \
	|| ( mv $(RUN_H_TMP) $(RUN_H); echo 'Updated CPP Options' )

$(BIN_DIR):
	mkdir $(BIN_DIR)

# Create an executable
gcm $(BIN_DIR)/$(RUN).exe: $(RUNDECK) $(OBJS)
	$(F90) $(LFLAGS) $(OBJS) $(LIBS) -o $(RUN).exe  $(LINK_OUTPUT)
ifdef LINK_OUTPUT
	@if [ -s $(RUN).ERR ] ; then cat $(RUN).ERR; else rm -f $(RUN).ERR; fi
endif
	@if [ ! -d $(BIN_DIR) ]; then mkdir $(BIN_DIR); fi
	mv $(RUN).exe $(BIN_DIR)
	@echo '--------------------------------------------------------------'
	@echo '---------    GCM successfuly compiled                ---------'
	@echo '---------    executable $(RUN).exe was created       ---------'
	@echo '--------------------------------------------------------------'

gcmlib $(BIN_DIR)/$(RUN).a: $(RUNDECK) $(MODS) $(OBJS) NCACC.o
	ar cr $(RUN).a $(OBJS) NCACC.o
	@if [ ! -d $(BIN_DIR) ]; then mkdir $(BIN_DIR); fi
	mv $(RUN).a $(BIN_DIR)


# Setup the gcm
setup:  $(BIN_DIR)/$(RUN).exe
	@echo '--------- Starting setup for $(RUN) --------------------------'
	@$(SETUP) $(RUN)

# Mkexe the gcm
exe:	$(BIN_DIR)/$(RUN).exe
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

htmldoc:
	[ -d $(DOC_DIR) ] || mkdir $(DOC_DIR)
	$(SCRIPTDIR)/gcmdoc.pl -D $(DOC_DIR) $(FSRCS)

clean:
	rm -f *.LST *~ *.ERR *.CHK *.L

vclean:
	rm -f *.o *.mod *.vo *.inc .depend* *.sig

