# Makefile for Model E Atmosphere GCM B399+ 
.SUFFIXES:
.SUFFIXES: .o .U .S .GCM .OCN .f .F .mod

######  Some user customizable settings:   ########

# EXTRA_FFLAGS specifies some extra flags you want to pass 
#to Fortarn compiler, like
# -g        - include debugging information
# -listing  - create listings (.L)
EXTRA_FFLAGS = -listing

# PROC_OUTPUT specifies if you want to send error output 
# to the screen or to a file ( no arg = screen )
#PROC_OUTPUT =   
PROC_OUTPUT =  > $*.ERR 2>&1 

#######       End of user settings       #########


# General Definitions

F90 = f90
F       = fco2_90
U	= uco2_f90
SETUP	= setup_f90
MKEXE	= mkexe_f90
RUNSRC	= E001
RUN	= E001dummy       #$(RUNSRC)
EXECDIR = /u/cmrun/modelE/exec
LIBS	= -L/u/cmrun -lGCM -lgP 
DEPENDFILE = .depend.$(RUN)

# Compiler flags

FFLAGS = -O2 -64 -mips4 -static -OPT:reorg_comm=off -w2 $(EXTRA_FFLAGS)

# Lists of objects and dependencies
# .depend is created with "make depend" command

sinclude $(DEPENDFILE)
FSRCS = $(OBJ_LIST:=.f)
OBJS = $(OBJ_LIST:=.o)


# Implicit rules

# Standard fortran
# .timestemp is a hack to set proper times on .o and .mod
.f.o:
	@touch .timestamp
	$(F90) -c $(FFLAGS) $<  $(PROC_OUTPUT)
	@touch -r .timestamp $@

.F.o:
	@touch .timestamp
	$(F90) -cpp -c $(FFLAGS) $<  $(PROC_OUTPUT)
	@touch -r .timestamp $@

# Update files
.U.o:
	$(EXECDIR)/$(U)  $*

# GISS-Fortran source files (with line numbers)
.S.o:
	$(EXECDIR)/$(F) $<

# MAPS.GCM
.GCM.o:
	$(EXECDIR)/$(F) $<

# FUNTABLE.OCN
.OCN.o:
	$(EXECDIR)/$(F) $<


#########             TARGETS               ########

# Extract list of object files from the rundeck $(RUN).R
# and create dependencies
depend:
	perl -e 'print "OBJ_LIST = "; while(<> !~ /^Object +modules/ ){};' \
	-e 'while(($$_=<>) !~ /(^Data +input)|(^ *$$)/)' \
	-e '{chop; s/ *(!.*)*$$//;print "\\\n$$_";}print"\n\n";' \
	< $(RUN).R > $(DEPENDFILE)
	make depend_depend RUN=$(RUN)

depend_depend:
	sfmakedepend -H -f $(DEPENDFILE) $(FSRCS)

# Compile the gcm
gcm:	$(OBJS)

# Make a current rundeck
mkrun: $(RUNSRC).R
	if [ $(RUN) != $(RUNSRC) ]; then  \
	sed 's/^'$(RUNSRC)'\([:\. ]\)/'$(RUN)'\1/g' $(RUNSRC).R > $(RUN).R ; \
	fi;

# Setup the gcm
setup:	gcm
	$(EXECDIR)/$(SETUP) $(RUN)

# Mkexe the gcm
mkexe:	gcm
	$(EXECDIR)/$(MKEXE) $(RUN)

newstart:
	rm -i /u/cmrun/$(RUN)/*

clean:
	rm -f *.LST *~ *.ERR *.CHK *.L

vclean:
	rm -f *.o *.mod

test:
	echo $(RUN) $(EXECDIR) $(F) $(U) $(SETUP)

cmp:    CMPE001.f
	f90 -o CMPE001 CMPE001.f

