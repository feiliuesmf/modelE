# Makefile for Model E Atmosphere GCM B399+ 
# Begin macro definitions

F	= fco2_f90
U	= uco2_f90
SETUP	= setup_f90
MKEXE	= mkexe_f90
RUNSRC	= E001
RUN	= E001dummy       #$(RUNSRC)
EXECDIR = /u/cmrun/modelE/exec
LIBS	= -L/u/cmrun -lGCM -lgP 
# Model common blocks 
MCOM	= E001M12_COM.o const.o
GEOCOM	= GEOM_B.o
QUSCOM	= SOMTQ_COM.o
GHYCOM  = GHYCOM.o
DAGCOM  = DAGCOM.o
RADCOM  = RADNCB.o
PBLCOM  = PBLCOM.o
CLDCOM  = CLD01_COM_E001.o
LKCOM   = LAKES_COM.o
COM	= $(MCOM) $(QUSCOM) $(GEOCOM)
RCOM    = BR00B.COM # will eventually be contained in RE001.f,RADNCB.f

# Model object files
MOBJ	= ME001M12.o
DYOBJ   = DYNE001.o
QOBJ    = QUSEM12.o SUBSIDEM12.o
DOBJ	= DE001M12.o
PBLOBJ	= PBLE001.o
PBLDRV	= PBLDRV.o
SLOBJ	= SLE001.o 
POBJ	= PE001M12.o 
EOBJ	= EE001M12.o
SOBJ	= SE001M12.o 
COBJ1   = CLD01.o
COBJ2   = CLD01_DRV_E001.o
ROBJ	= RE001.o
OOBJ	= OCNE001.o
UOBJ    = UTILDBL.o
FFTOBJ  = FFT72.o
RNDOBJ  = RANDSGI.o
OUTOBJ  = POUT.o
LOBJ    = LAKES.o
LIOBJ   = LANDICE.o
OBJ	= $(MOBJ) $(DYOBJ) $(QOBJ) $(DOBJ) $(PBLOBJ) $(SLOBJ) $(POBJ) \
	  $(SOBJ) $(COBJ1) $(COBJ2) $(CLDCOM) $(ROBJ) $(UOBJ) $(OOBJ)  \
	  $(COM) $(EOBJ) $(GHYCOM) $(DAGCOM) $(RADCOM) $(PBLCOM) $(PBLDRV) \
          $(FFTOBJ) $(RNDOBJ) $(OUTOBJ) $(LKCOM) $(LOBJ) $(LIOBJ)

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
$(QUSCOM):	$(MCOM)
$(GEOCOM):	$(MCOM)
$(GHYCOM):	$(MCOM) $(SLOBJ)
$(DAGCOM):	$(MCOM) $(RADCOM)
$(RADCOM):	$(MCOM)
$(CLDCOM):	$(MCOM) $(COBJ1)
$(PBLCOM):	$(MCOM) $(PBLOBJ)
$(LKCOM):	$(MCOM)
$(DYOBJ):	$(COM) $(PBLCOM) $(OOBJ) $(DAGCOM)
$(MOBJ):	$(COM) $(GHYCOM) $(DAGCOM) $(RNDOBJ) $(CLDCOM) $(PBLOBJ) $(DYOBJ) $(RADCOM) $(PBLCOM) $(OOBJ) $(UOBJ) $(LOBJ)
$(POBJ):	$(COM) $(DAGCOM) $(ROBJ)  $(RNDOBJ) $(CLDCOM) $(RADCOM) $(PBLCOM) $(DYOBJ) $(OOBJ) $(UOBJ) $(LOBJ) $(LIOBJ)
$(PBLDRV):	$(COM) $(PBLOBJ) $(PBLCOM) $(DYOBJ) $(OOBJ) $(UOBJ)
$(SOBJ):	$(COM) $(DAGCOM) $(PBLOBJ) $(RADCOM) $(PBLCOM) $(DYOBJ) $(OOBJ)
$(EOBJ):	$(COM) $(DAGCOM) $(GHYCOM) $(CLDCOM) $(PBLOBJ) $(RADCOM) $(PBLCOM) $(DYOBJ) $(UOBJ)
$(DOBJ):	$(COM) $(DYOBJ) $(DAGCOM) $(RADCOM) $(PBLCOM) $(OOBJ)
$(QOBJ):	$(COM) 
$(ROBJ):	$(RCOM)
$(OOBJ):	$(COM) $(UOBJ)
$(COBJ1):	$(COM)
$(COBJ2):	$(COM) $(COBJ1) $(CLDCOM) $(OOBJ) $(DAGCOM) $(RNDOBJ) $(DYOBJ) $(PBLCOM)
$(LOBJ):	$(COM) $(LKCOM) $(DAGCOM)

# Compile the gcm
gcm:	$(OBJ)

# Make a current rundeck
$(RUN).R:	$(RUNSRC).R
mkrun:	$(RUN).R
	mkrun $(RUNSRC) $(RUN) 

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
