#
# this file contains rules shared by Makefiles in "model" and "aux"
#

.PHONY: FORCE


######  Some user customizable settings:   ########

# EXTRA_FFLAGS specifies some extra flags you want to pass 
# to Fortarn compiler, like
# -g        - include debugging information
# -listing  - create listings (.L)
EXTRA_FFLAGS =

# CPPFLAGS specifies some pre-processor directives you want to pass 
# to Fortarn compiler, like
# -DCOMPILER=Absoft    - compile code that is specific to Absoft
CPPFLAGS =

ifeq ($(OUTPUT_TO_FILES),YES)
  COMP_OUTPUT = > $*.ERR 2>&1 || { r=$$? ; cat $*.ERR ; exit $$r ; }
  LINK_OUTPUT = > $(RUN).ERR 2>&1 || { r=$$? ; cat $(RUN).ERR ; exit $$r ; }
else
  COMP_OUTPUT =
  LINK_OUTPUT =
endif

# if -s specified enable some extra messages
ifeq ($(findstring s,$(MFLAGS)),s)
MSG = 
else
MSG = > /dev/null
endif

#
# starting machine - specific options
#

NO_COMMAND = @echo Requested target is not supported on $(UNAME); exit 1;
F90 = $(NO_COMMAND)
FMAKEDEP = $(NO_COMMAND)
F =  $(NO_COMMAND)
U = $(NO_COMMAND)
SETUP = $(SCRIPTS_DIR)/setup_e
CPP = $(NO_COMMAND)
MACHINE = not_specified
COMPILER = not_specified
LIBS =

UNAME = $(shell uname)

# SGI - specific options here
ifeq ($(UNAME),IRIX64)
MACHINE = SGI
F90 = f90
CPP = /lib/cpp -P
FMAKEDEP = $(SCRIPTS_DIR)/sfmakedepend -H
F       = $(SCRIPTS_DIR)/fco2_90
U	= $(SCRIPTS_DIR)/uco2_f90
FFLAGS = -cpp -O2 -64 -mips4 -OPT:reorg_comm=off -w2 -OPT:Olimit=5745
LFLAGS = -64 -O2 -mips4 -lfastm -mp -OPT:reorg_common=OFF -Wl,-woff,134 -Wl,-woff,15
endif

# Linux - specific options here
ifeq ($(UNAME),Linux)
MACHINE = Linux
#COMPILER = Vast
#F90 = f90
#FMAKEDEP = $(SCRIPTS_DIR)/sfmakedepend -m vo
# this is for Intel Compiler
#COMPILER = Intel
#F90 = ifc
#FMAKEDEP = $(SCRIPTS_DIR)/sfmakedepend -m d -int
#FFLAGS = -O2
#LFLAGS = -O2
COMPILER = Absoft
F90 = f90
CPP = /usr/bin/cpp -P -traditional
FMAKEDEP = $(SCRIPTS_DIR)/sfmakedepend -h
CPPFLAGS= -DCONVERT_BIGENDIAN
FFLAGS = -O -B100 
LFLAGS = -lf90math -lV77 -lU77
endif

# IBM - specific options here
ifeq ($(UNAME),AIX)
MACHINE = IBM
F90 = xlf90_r
CPP = /lib/cpp -P
FMAKEDEP = perl $(SCRIPTS_DIR)/sfmakedepend
FFLAGS = -O2 -qfixed
LFLAGS = -O2
endif

# DEC Alpha - specific options here
# This option has not been tested yet. It may need some adjustment.
ifeq ($(UNAME),OSF1)
MACHINE = DEC
F90 = f90
CPP = /lib/cpp -P
FMAKEDEP = $(SCRIPTS_DIR)/sfmakedepend
CPPFLAGS= -DCONVERT_BIGENDIAN
FFLAGS = -O2
LFLAGS = -O2
endif


CPPFLAGS += -DMACHINE_$(MACHINE)

# end of machine - specific options

#
# Check for extra options specified in modelErc
#

ifdef NETCDFHOME
  LIBS += -L$(NETCDFHOME) -lnetcdf
endif

#
# Pattern  rules
#

FORCE:

*.mod: FORCE

%.mod: 
	@echo -n checking $@:
#	@echo 'called rule for $@, depends on $^ built from: '`cat $@.sig`
	@if [ ! -s $@.sig ] || [ "`cat $@.sig`" != "$<" ]; then \
	echo "  dependencies for $@ have changed - recompiling"; \
	rm -f $<; $(MAKE) $< RUN=$(RUN); else echo '  ok'; fi

# Standard fortran
# .timestemp is a hack to set proper times on .o and .mod
# For the Absoft compiler, we need to force a cpp run through
%.o: %.f
	@echo -n compiling $< ... $(MSG)
	@touch .timestamp
ifeq ($(COMPILER),Absoft)
	cp $*.f $*.F
	$(F90) -c $(FFLAGS) $(EXTRA_FFLAGS) $(CPPFLAGS) $*.F  $(COMP_OUTPUT)
	rm -f $*.F
else
	$(F90) -c $(FFLAGS) $(EXTRA_FFLAGS) $(CPPFLAGS) $*.f  $(COMP_OUTPUT)
endif
	-@if [ `ls | grep ".mod" | tail -1` ] ; then for i in *.mod; \
	  do if [ ! -s $$i.sig ] || [ $$i -nt $$i.sig ] ; then \
	  echo $@ > $$i.sig; fi; done; fi 
	@touch -r .timestamp $@
	@if [ -s $*.ERR ] ; then echo $(MSG); else echo Done $(MSG); fi
ifdef COMP_OUTPUT
	@if [ -s $*.ERR ] ; then cat $*.ERR; else  rm -f $*.ERR; fi
endif

# cpp preprocessing
%.cpp: %.f
	@echo preprocessing $<  $(MSG)
	$(CPP) $(CPPFLAGS) $*.f > $*.cpp

# Update files
%.o: %.U
	$(U)  $*

# GISS-Fortran source files (with line numbers) 
%.o: %.S
	$(F) $<

# MAPS.GCM
%.o: %.GCM
	$(F) $<

# FUNTABLE.OCN
%.o: %.OCN
	$(F) $<


# end of Pattern  rules

