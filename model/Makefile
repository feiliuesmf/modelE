# Makefile for Model E Atmosphere GCM B399+ 
# Begin macro definitions

F	= fco2_f90
U	= uco2_f90
SETUP	= setup_f90
MKEXE	= mkexe_f90
RUN	= E001
EXECDIR = /u/cmrun/modelE/exec
LIBS	= -L/u/cmrun -lGCM -lgP 
# Model common blocks 
COM	= E001M12.COM SOMTQ.COM
SCOM    = soils45.COM 
RCOM    = BR00B.COM

# Model object files
MOBJ	= ME001M12.o
QOBJ    = QUSEM12.o SUBSIDEM12.o
DOBJ	= DE001M12.o
PBLOBJ	= PBLE001.o
SLOBJ	= SLE001.o 
POBJ	= PE001M12.o 
SOBJ	= EE001M12.o SE001M12.o 
COBJ    = CE001M12.o CE001M12.o
ROBJ	= RE001.o
UOBJ    = FFT72.o  UTILDBL.o  RANDSGI.o   
OBJ	= $(MOBJ) $(QOBJ) $(DOBJ) $(PBLOBJ) $(SLOBJ) $(POBJ) $(SOBJ) $(COBJ) $(ROBJ) $(UOBJ)

# Make rules for generating object files from GISS-Fortran files
.SUFFIXES:
.SUFFIXES: .o .U .S .GCM .OCN .f .F

# Standard fortran
.f.o:
	$(EXECDIR)/$(F) $<

# Update files
.U.o:
	$(EXECDIR)/$(U)  $*

# GISS-Fortran source files (with line numbers)
.S.o:
	$(EXECDIR)/$(F) $<

.F.o:
	$(EXECDIR)/$(F) $<

# MAPS.GCM
.GCM.o:
	$(EXECDIR)/$(F) $<


# FUNTABLE.OCN
.OCN.o:
	$(EXECDIR)/$(F) $<

# Make object files dependent on relevant common block and *.S source files
$(MOBJ):	$(COM) 
$(POBJ):	$(COM) $(SCOM) 
$(PBLOBJ):	$(COM) 
$(SOBJ):	$(COM) $(SCOM) 
$(DOBJ):	$(COM) 
$(QOBJ):	$(COM) 
$(COBJ):	$(COM) 
$(ROBJ):	$(RCOM)

# Compile the gcm
gcm:	$(OBJ)

# Setup the gcm
setup:	gcm
	$(EXECDIR)/$(SETUP) $(RUN)

# Mkexe the gcm
mkexe:	gcm
	$(EXECDIR)/$(MKEXE) $(RUN)

newstart:
	rm -i /u/cmrun/$(RUN)/*

clean:
	rm -f *.LST *~ *.ERR *.CHK

vclean:
	rm -f $(OBJ) 

test:
	echo $(RUN) $(EXECDIR) $(F) $(U) $(SETUP)

cmp:    CMPE001.f
	f90 -o CMPE001 CMPE001.f
