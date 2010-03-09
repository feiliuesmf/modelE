SHELL=/bin/sh

#FC=    ifort -assume byterecl -convert big_endian
#LD=    ifort -assume byterecl -convert big_endian
FC=     ifort 
LD=     ifort

#FFLAGS = -traceback
FFLAGS = 
LDFLAGS =  
CPPFLAGS =
LIBS =

avg_obj=const_proc.o hycom_dimen.o hycom_arrays.o avg.o gethyb.o extrct.o \
	prtmsk.o bigrid.o gtdpth.o meshsz.o restep.o findmx.o thruflow.o
avg_src = $(avg_obj:.o=.f)

ov2_obj=const_proc.o hycom_dimen.o hycom_arrays.o pltov_mon.o gethyb.o \
	extrct.o prtmsk.o bigrid.o gtdpth.o meshsz.o restep.o \
	findmx.o thruflow.o

latlonz_obj=const_proc.o hycom_dimen.o hycom_arrays.o cpler.o \
	latlonz.o bigrid.o gethyb.o extrct.o gtdpth.o restep.o grdtrns.o 
latlonz_src = $(latlonz_obj:.o=.f)

get2d_obj=const_proc.o hycom_dimen.o get2d.o 
get2d_src = $(get2d_obj:.o=.f)

avg:    $(avg_obj)
	$(LD) -o avg $(avg_obj) $(LIBS) $(LDFLAGS)

pltov_mon:  $(ov2_obj)
	$(LD) -o pltov_mon $(ov2_obj) $(LIBS) $(LDFLAGS)

latlonz:    $(latlonz_obj)
	$(LD) -o latlonz $(latlonz_obj) $(LIBS) $(LDFLAGS)

get2d:  $(get2d_obj)
	$(LD) -o get2d $(get2d_obj) $(LIBS) $(LDFLAGS)

avg.o: avg.f const_proc.f hycom_dimen.f
	$(LD) $(LDFLAGS) -c avg.f
bigrid.o: bigrid.f hycom_dimen.f
	$(LD) $(LDFLAGS) -c bigrid.f
const_proc.o: const_proc.f
	$(LD) $(LDFLAGS) -c const_proc.f
extrct.o: extrct.f const_proc.f hycom_dimen.f
	$(LD) $(LDFLAGS) -c extrct.f
get2d.o: get2d.f const_proc.f hycom_dimen.f
	$(LD) $(LDFLAGS) -c get2d.f
gethyb.o: gethyb.f const_proc.f hycom_dimen.f hycom_arrays.f
	$(LD) $(LDFLAGS) -c gethyb.f
gtdpth.o: gtdpth.f const_proc.f hycom_dimen.f
	$(LD) $(LDFLAGS) -c gtdpth.f
grdtrns.o: grdtrns.f
	$(LD) $(LDFLAGS) -c grdtrns.f
hycom_dimen.o: hycom_dimen.f
	$(LD) $(LDFLAGS) -c hycom_dimen.f
hycom_arrays.o: hycom_arrays.f
	$(LD) $(LDFLAGS) -c hycom_arrays.f
latlonz.o: latlonz.f const_proc.f hycom_dimen.f hycom_arrays.f
	$(LD) $(LDFLAGS) -c latlonz.f
restep.o: restep.f const_proc.f hycom_dimen.f
	$(LD) $(LDFLAGS) -c restep.f

tar: makefile input readme $(latlonz_src) $(get2d_src) $(avg_src)
	tar -cvf proc.tar makefile input readme $(latlonz_src) $(get2d_src) $(avg_src)