.SUFFIXES:
all:  ESMF_MOD.o ESMF_CUSTOM_MOD.o FieldComm.o assert.o

# dependencies
ESMF_MOD.o: 		ESMF_MOD.F90
assert.o:		assert.F90 ESMF_MOD.o  assert.h
ESMF_CUSTOM_MOD.o: 	ESMF_CUSTOM_MOD.F90 assert.o ESMF_MOD.o assert.h

clean:
	rm *.o *.mod *.sig

include ../Rules.make

