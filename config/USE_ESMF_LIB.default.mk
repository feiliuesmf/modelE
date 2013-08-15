# these are the options needed to compile the code with ESMF library
# (MPI support should be included separately)

# conditional compilation and linking with ESMF
ifdef ESMFMKFILE
  ifneq ($(ESMFMKFILE),)
    include $(ESMFMKFILE)
    CPPFLAGS += -DUSE_ESMF_LIB
    FFLAGS   += $(ESMF_F90COMPILEPATHS)
    F90FLAGS += $(ESMF_F90COMPILEPATHS)
    VPATH    += ${ESMFINCLUDEDIR}
    CFLAGS   += $(ESMF_CXXCOMPILEPATHS)
    LIBS     += $(ESMF_F90LINKPATHS) $(ESMF_F90LINKRPATHS) $(ESMF_F90ESMFLINKLIBS)
  else
    $(error environment variable ESMFMKFILE cannot be empty and must point to a esmf.mk file)
  endif
else
  $(error environment variable ESMFMKFILE must be set when USE_ESMF_LIB is YES)
endif
