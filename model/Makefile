# Makefile for Model E Atmosphere GCM B399+ 
.SUFFIXES:
.SUFFIXES: .o .U .S .GCM .OCN .f .F

.PHONY: mkrun depend gcm setup clean vclean newstart mkexe cmp .R depend_depend


######  Some user customizable settings:   ########

# EXTRA_FFLAGS specifies some extra flags you want to pass 
# to Fortarn compiler, like
# -g        - include debugging information
# -listing  - create listings (.L)
EXTRA_FFLAGS = -listing

# COMP_OUTPUT, LINK_OUTPUT specify if you want to send error output 
# to the screen or to a file ( no arg = screen )
#COMP_OUTPUT =   
#LINK_OUTPUT =
COMP_OUTPUT = > $*.ERR 2>&1
LINK_OUTPUT = > $(RUN).ERR 2>&1

#######       End of user settings       #########

# Print Usage if no targets were given
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
	@echo '"setup" currently works only on SGI's.'
	@echo
	@echo The following is a list of supported targets:
	@echo
	@echo ' mkrun   - create new rundeck'
	@echo ' depend  - create dependencies for specified rundeck'
	@echo ' gcm     - compile object files and build executable for specified rundeck'
	@echo ' setup   - do setup for specified rundeck'
	@echo ' clean   - remove listings and error files'
	@echo ' vclean  - remove object files, .mod files and dependencies'
	@echo ' newstart - remove all files in the run directory'
	@echo ' mkexe   - compile gcm and put executable into RUN directory'
	@echo ' cmp     - compile CMPE001'
	@echo
	@echo '          Recommended usage:'
	@echo 'The targets are made in order: "mkrun", "depend", "gcm", "setup",'
	@echo 'but since Makefile automatically creates necessary targets, you'
	@echo 'do not have to make all of them manually.'
	@echo 'You can skip "mkrun" step unless you want to edit the rundeck'
	@echo 'before compiling gcm. You may (and should) always skip "depend".'
	@echo 'So if you want to compile gcm and set up a run with a name "my_run"'
	@echo 'typically you will have to execute just two commands:'
	@echo
	@echo '       gmake gcm RUN=my_run'
	@echo '       gmake setup RUN=my_run'
	@echo
	@echo 'You can actually skip the first line either and gcm will be'
	@echo 'compiled automatically, but you may not want to do it in most
	@echo 'cases, since you want to look at the warning massages before
	@echo 'doing "setup".'
	@echo 'Rundeck is generated automatically only if corresponding '
	@echo 'file (*.R) does not exist. If it does exist, it will not be'
	@echo 'regenerated even if the source rundeck is newer. This prevents'
	@echo 'accidental overwriting of a rundeck you created or edited'
	@echo 'manually. Use "mkrun" to overwrite it. '
	@echo
 

# General Definitions

RUNSRC  = E001
RUN     =
DEPENDFILE = .depend.$(RUN)

SCRIPTDIR = ./SCRIPTS

NO_COMMAND = @echo Requested target is not supported on $(UNAME); exit 1;

F90 = $(NO_COMMAND)
FMAKEDEP = $(NO_COMMAND)
F =  $(NO_COMMAND)
U = $(NO_COMMAND)
SETUP = $(NO_COMMAND)
MACHINE = not_specified

UNAME = $(shell uname)

# SGI - specific options here
ifeq ($(UNAME),IRIX64)
MACHINE = SGI
F90 = f90
FMAKEDEP = $(SCRIPTDIR)/sfmakedepend -H
F       = $(SCRIPTDIR)/fco2_90
U	= $(SCRIPTDIR)/uco2_f90
SETUP	= $(SCRIPTDIR)/setup_e
LIBS	= -L/u/cmrun -lGCM -lgP 
FFLAGS = -O2 -64 -mips4 -static -OPT:reorg_comm=off -w2 
LFLAGS = -64 -O2 -mips4 -static -lfastm -mp -OPT:reorg_common=OFF -Wl,-woff,134 -Wl,-woff,15
endif

# Linux - specific options here
ifeq ($(UNAME),Linux)
MACHINE = Linux
F90 = f90
FMAKEDEP = $(SCRIPTDIR)/sfmakedepend -m vo
FFLAGS = -O2
LFLAGS = -O2
endif

# IBM - specific options here
ifeq ($(UNAME),AIX)
MACHINE = IBM
F90 = xlf90_r
FMAKEDEP = perl $(SCRIPTDIR)/sfmakedepend
FFLAGS = -O2 -qfixed
LFLAGS = -O2
endif

# end of machine - specific options


# Lists of objects and dependencies
# .depend is created with "make depend" command

ifneq ($(RUN).R,.R)
sinclude $(DEPENDFILE)
endif
FSRCS = $(OBJ_LIST:=.f)
OBJS = $(OBJ_LIST:=.o)


# Implicit rules

%.mod:
	@echo 'hack to make dependencies work' >/dev/null

# Standard fortran
# .timestemp is a hack to set proper times on .o and .mod
.f.o:
	@touch .timestamp
	$(F90) -c $(FFLAGS) $(EXTRA_FFLAGS) $<  $(COMP_OUTPUT)
	@touch -r .timestamp $@
ifdef COMP_OUTPUT
	@if [ -s $*.ERR ] ; then cat $*.ERR; else rm -f $*.ERR; fi
endif

.F.o:
	@touch .timestamp
	$(F90) -cpp -c $(FFLAGS) $(EXTRA_FFLAGS) $<  $(COMP_OUTPUT)
	@touch -r .timestamp $@
ifdef COMP_OUTPUT
	@if [ -s $*.ERR ] ; then cat $*.ERR; else rm -f $*.ERR; fi
endif

# Update files
.U.o:
	$(U)  $*

# GISS-Fortran source files (with line numbers)
.S.o:
	$(F) $<

# MAPS.GCM
.GCM.o:
	$(F) $<

# FUNTABLE.OCN
.OCN.o:
	$(F) $<


#########             TARGETS               ########

# Make a current rundeck
mkrun $(RUN).R:
	@if [ ! -s $(RUNSRC).R ]; then \
	echo 'Source rundeck $(RUNSRC).R does not exist'; exit 1; fi
	@if [ "$(RUN)" = "" ]; then \
	echo 'You should specify run name on a command line (i.e. RUN=...)'; \
	exit 1; fi
	@echo '---- Regenerating Rundeck $(RUN).R from default ($(RUNSRC).R) ----'
	@if [ $(RUN) != $(RUNSRC) ]; then  \
	sed 's/^'$(RUNSRC)'\([:\. ]\)/'$(RUN)'\1/g' $(RUNSRC).R > $(RUN).R ; \
	fi;
	@echo '---------         Rundeck  created                   ---------' 

# Extract list of object files from the rundeck $(RUN).R
# and create dependencies
depend $(DEPENDFILE):  $(RUN).R $(FSRCS)
	@echo '-------- Dependecies are out of date or do not exist ---------'
	@echo '---------         regenerating them                  ---------'
	@perl -e 'print "OBJ_LIST = "; while(<> !~ /^Object +modules/ ){};' \
	-e 'while(($$_=<>) !~ /(^Data +input)|(^ *$$)/)' \
	-e '{chop; s/ *(!.*)*$$//; s/~/$(MACHINE)/g;print "\\\n$$_";}print"\n\n";' \
	< $(RUN).R > $(DEPENDFILE)
	@$(MAKE) -s depend_depend RUN=$(RUN)

depend_depend:
	$(FMAKEDEP) -f $(DEPENDFILE) $(FSRCS)
	@echo '---------         done with dependencies             ---------'

# Create an executable
gcm $(RUN).exe: $(RUN).R $(OBJS)
	$(F90) $(LFLAGS) $(OBJS) -o $(RUN).exe  $(LINK_OUTPUT)
ifdef LINK_OUTPUT
	@if [ -s $(RUN).ERR ] ; then cat $(RUN).ERR; else rm -f $(RUN).ERR; fi
endif
	@echo '--------------------------------------------------------------'
	@echo '---------    GCM successfuly compiled                ---------'
	@echo '---------    executable $(RUN).exe was created       ---------'
	@echo '--------------------------------------------------------------'

# Setup the gcm
setup:  $(RUN).exe
	@echo '--------- Starting setup for $(RUN) --------------------------'
	@$(SETUP) $(RUN)

# Mkexe the gcm
mkexe:	$(RUN).exe
	@if [ ! -s /u/cmrun/$(RUN)/I ] ; then \
	echo 'directory /u/cmrun/$(RUN) or input files inside it are missing';\
	echo 'you should make "setup" first'; \
	exit 1; fi
	mv -f $(RUN).exe /u/cmrun/$(RUN)/

newstart:
	@if [ "$(RUN)" = "" ]; then \
        echo 'You should specify run name on a command line (i.e. RUN=...)'; \
        exit 1; fi
	rm -i /u/cmrun/$(RUN)/*

clean:
	rm -f *.LST *~ *.ERR *.CHK *.L

vclean:
	rm -f *.o *.mod *.vo *.inc .depend*

cmp:    CMPE001.f
	$(F90) $(FFLAGS) -o CMPE001 CMPE001.f


