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
COM	= E001M12_COM.o SOMTQ_COM.o GEOM_B.o
SCOM    = soils45.COM 
RCOM    = BR00B.COM

# Model object files
MOBJ	= ME001M12.o
QOBJ    = QUSEM12.o SUBSIDEM12.o
DOBJ	= DE001M12.o
PBLOBJ	= PBLE001.o
SLOBJ	= SLE001.o 
POBJ	= PE001M12.o 
EOBJ	= EE001M12.o
SOBJ	= SE001M12.o 
COBJ    = CE001M12.o
ROBJ	= RE001.o
UOBJ    = UTILDBL.o
FFTOBJ  = FFT72.o
RNDOBJ  = RANDSGI.o
OBJ	= $(MOBJ) $(QOBJ) $(DOBJ) $(PBLOBJ) $(SLOBJ) $(POBJ) $(SOBJ) $(COBJ) $(ROBJ) $(UOBJ) \
          $(FFTOBJ) $(RNDOBJ) $(COM) $(EOBJ)

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

# Make object files dependent on relevant modules
$(MOBJ):	$(COM) $(SLOBJ) $(RNDOBJ) $(COBJ) $(PBLOBJ)
$(POBJ):	$(COM) $(ROBJ)  $(RNDOBJ) $(COBJ)
$(PBLOBJ):	$(COM) 
$(SOBJ):	$(COM) $(SLOBJ) $(PBLOBJ)
$(EOBJ):	$(COM) $(SLOBJ) $(COBJ) $(PBLOBJ)
$(SLOBJ):	$(COM) $(SCOM)
$(DOBJ):	$(COM)
$(QOBJ):	$(COM) 
$(COBJ):	$(COM) $(RNDOBJ)
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

# need to alter vclean to remove relevant .mod files
vclean:
	rm -f $(OBJ) 

test:
	echo $(RUN) $(EXECDIR) $(F) $(U) $(SETUP)

cmp:    CMPE001.f
	f90 -o CMPE001 CMPE001.f
