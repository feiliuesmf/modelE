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

# COMP_OUTPUT, LINK_OUTPUT specify if you want to send error output 
# to the screen or to a file ( no arg = screen )
#COMP_OUTPUT =   
#LINK_OUTPUT =
COMP_OUTPUT = > $*.ERR 2>&1
LINK_OUTPUT = > $(RUN).ERR 2>&1

#
# starting machine - specific options
#

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
FFLAGS = -cpp -O2 -64 -mips4 -OPT:reorg_comm=off -w2 
LFLAGS = -64 -O2 -mips4 -lfastm -mp -OPT:reorg_common=OFF -Wl,-woff,134 -Wl,-woff,15
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

#
# Pattern  rules
#

FORCE:

*.mod: FORCE

%.mod: 
	@echo -n checking $@:
#	@echo 'called rule for $@, depends on $^ built from: '`cat $@.sig`
	@if [ "`cat $@.sig`" != "$<" ]; then \
	echo "  dependencies for $@ have changed - recompiling"; \
	rm -f $<; $(MAKE) $< RUN=$(RUN); else echo '  ok'; fi

# Standard fortran
# .timestemp is a hack to set proper times on .o and .mod
%.o: %.f
	@touch .timestamp
	$(F90) -c $(FFLAGS) $(EXTRA_FFLAGS) $<  $(COMP_OUTPUT)
	-@if [ -s *.mod ] ; then for i in *.mod; do if [ ! -s $$i.sig ] || [ $$i -nt $$i.sig ] ; then echo $@ > $$i.sig; fi; done; fi
	@touch -r .timestamp $@
ifdef COMP_OUTPUT
	@if [ -s $*.ERR ] ; then cat $*.ERR; else rm -f $*.ERR; fi
endif

%.o: %.F
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

# RE001.f
RE001.o: RE001.f
	@touch .timestamp
	$(F90) -c -static $(FFLAGS) $(EXTRA_FFLAGS) $<  $(COMP_OUTPUT)
	-@if [ -s *.mod ] ; then for i in *.mod; do if [ ! -s $$i.sig ] || [ $$i -nt $$i.sig ] ; then echo $@ > $$i.sig; fi; done; fi
	@touch -r .timestamp $@
ifdef COMP_OUTPUT
	@if [ -s $*.ERR ] ; then cat $*.ERR; else rm -f $*.ERR; fi
endif



# end of Pattern  rules
