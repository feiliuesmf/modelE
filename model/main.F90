program main
#ifdef USE_ESMF_LIB
  use main_nuopc_mod
  call modelE_NUOPC_mainDriver()
#else
  use MODELE
  call modelE_mainDriver()
#endif
end program main
