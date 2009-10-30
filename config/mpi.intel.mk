
#LIBS += -lcprts -limf -lm -lcxa -lunwind -lrt -ldl \
#-lfmpi -lmpi -lstdc++ -threads

ifeq ($(IFORT_RELEASE),11.1)
LIBS += -limf -lm -lrt -ldl \
-lmpiif -lmpi -lstdc++ -threads
else
LIBS += -lcprts -limf -lm -lcxa -lunwind -lrt -ldl \
-lmpiif -lmpi -lstdc++ -threads
endif

#-lmpigf -lmpi -lstdc++ -threads


#LIBS +=  -lmpi -lmpigc3 -lmpigc4 -lmpigf -lmpigi -lmpiic4 -lmpiic -lmpiif -threads

