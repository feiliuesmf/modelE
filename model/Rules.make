#
# this files contains rules shared by Makefiles in "model" and "aux"
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

# COMP_OUTPUT, LINK_OUTPUT specify if you want to send error output 
# to the screen or to a file ( no arg = screen )
#COMP_OUTPUT =   
#LINK_OUTPUT =
COMP_OUTPUT = > $*.ERR 2>&1 || { r=$$? ; cat $*.ERR ; exit $$r ; }
LINK_OUTPUT = > $(RUN).ERR 2>&1 || { r=$$? ; cat $(RUN).ERR ; exit $$r ; }

# overwriting above options if environment var MODELE_MAKE_OUTPUT=SCREEN
ifeq ($(MODELE_MAKE_OUTPUT),SCREEN)
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
MACHINE = not_specified
COMPILER = not_specified

UNAME = $(shell uname)

# SGI - specific options here
ifeq ($(UNAME),IRIX64)
MACHINE = SGI
F90 = f90
FMAKEDEP = $(SCRIPTS_DIR)/sfmakedepend -H
F       = $(SCRIPTS_DIR)/fco2_90
U	= $(SCRIPTS_DIR)/uco2_f90
LIBS	= -L/usr/local/netcdf-3.4/lib64 -lnetcdf 
FFLAGS = -cpp -O2 -64 -mips4 -OPT:reorg_comm=off -w2 
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
FMAKEDEP = $(SCRIPTS_DIR)/sfmakedepend -h
CPPFLAGS= -DCONVERT_BIGENDIAN
FFLAGS = -O -B100 
LFLAGS = -lf90math -lV77 -lU77
endif

# IBM - specific options here
ifeq ($(UNAME),AIX)
MACHINE = IBM
F90 = xlf90_r
FMAKEDEP = perl $(SCRIPTS_DIR)/sfmakedepend
FFLAGS = -O2 -qfixed
LFLAGS = -O2
endif

# end of machine - specific options

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
	@echo compiling $@ $(MSG)
	@touch .timestamp
	@echo -n Compiling $<...
	@if [ $(COMPILER) = Absoft ] ; then cpp -traditional -E $(CPPFLAGS) $< $*.F ; $(F90) -c $(FFLAGS) $(EXTRA_FFLAGS) $*.F  $(COMP_OUTPUT) ; rm $*.F ; else $(F90) -c $(FFLAGS) $(EXTRA_FFLAGS) $(CPPFLAGS) $*.f  $(COMP_OUTPUT) ; fi
	-@if [ -s `ls *.mod | tail -1 ` ] ; then for i in *.mod; do if [ ! -s $$i.sig ] || [ $$i -nt $$i.sig ] ; then echo $@ > $$i.sig; fi; done; fi
	@touch -r .timestamp $@
ifdef COMP_OUTPUT
	@if [ -s $*.ERR ] ; then echo ; cat $*.ERR; else echo Done ; rm -f $*.ERR; fi
endif

%.o: %.F
	@echo compiling $@ $(MSG)
	@touch .timestamp
	$(F90) -c $(FFLAGS) $(EXTRA_FFLAGS) $<  $(COMP_OUTPUT)
	-@if [ -s *.mod ] ; then for i in *.mod; do if [ ! -s $$i.sig ] || [ $$i -nt $$i.sig ] ; then echo $@ > $$i.sig; fi; done; fi
	@touch -r .timestamp $@
ifdef COMP_OUTPUT
	@if [ -s $*.ERR ] ; then cat $*.ERR; else rm -f $*.ERR; fi
endif

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

# RADIATION.f
#RADIATION.o: RADIATION.f
#	@touch .timestamp
#	$(F90) -c -static $(FFLAGS) $(EXTRA_FFLAGS) $<  $(COMP_OUTPUT)
#	-@if [ -s *.mod ] ; then for i in *.mod; do if [ ! -s $$i.sig ] || [ $$i -nt $$i.sig ] ; then echo $@ > $$i.sig; fi; done; fi
#	@touch -r .timestamp $@
#ifdef COMP_OUTPUT
#	@if [ -s $*.ERR ] ; then cat $*.ERR; else rm -f $*.ERR; fi
#endif



# end of Pattern  rules

