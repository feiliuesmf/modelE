all:  ESMF_MOD.o ESMF_CUSTOM_MOD.o FieldComm.o assert.o

# dependencies
ESMF_MOD.o: 		ESMF_MOD.F90
assert.o:		assert.F90 ESMF_MOD.o  assert.h
FieldComm.o: 		FieldComm.F90 ESMF_MOD.o assert.o assert.h
ESMF_CUSTOM_MOD.o: 	ESMF_CUSTOM_MOD.F90 assert.o FieldComm.o ESMF_MOD.o assert.h

clean:
	rm *.o *.mod *.sig

# Free source form fortran
# .timestemp is a hack to set proper times on .o and .mod
# For the Absoft/Lahey/PGI compilers, we need to force a cpp run through
%.o: %.F90
	@echo $(ECHO_FLAGS)  compiling $< ... $(MSG) $(MACHINE) $(F90)
	@touch .timestamp
ifeq ($(MACHINE),MAC)
	$(CPP) $(CPPFLAGS) $< | sed -n '/^#pragma/!p' > $*_cpp.F90
	$(F90) -c $(FFLAGS) $(EXTRA_FFLAGS) $(RFLAGS) $*_cpp.F90 \
	  -o $@ $(COMP_OUTPUT)
	rm -f $*_cpp.F90
else
ifeq ($(COMPILER),Absoft)
	$(CPP) $(CPPFLAGS) $< $*_cpp.f
	$(F90) -c $(FFLAGS) $(EXTRA_FFLAGS) $(RFLAGS) $*_cpp.f -o $@ $(COMP_OUTPUT)
	rm -f $*_cpp.f
else
	@$(F90) -c $(FFLAGS) $(EXTRA_FFLAGS) $(CPPFLAGS) $(RFLAGS) $<  $(COMP_OUTPUT)
endif
endif
	-@if [ `ls | grep ".mod" | tail -1` ] ; then for i in *.mod; \
	  do if [ ! -s $$i.sig ] || [ `find $$i -newer $$i.sig` ] ; then \
	  echo $@ > $$i.sig; fi; done; fi 
	@touch -r .timestamp $@
	@if [ -s $*.ERR ] ; then echo $(MSG); else echo Done $(MSG); fi
ifdef COMP_OUTPUT
	@if [ -s $*.ERR ] ; then cat $*.ERR; else  rm -f $*.ERR; fi
endif



