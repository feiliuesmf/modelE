#ifdef USE_ESMF_LIB
module main_nuopc_mod
  !
  ! If USE_ESMF_LIB is turned on this main program is compiled.
  ! This is the NUOPC main program for modele
  ! Fei.Liu@noaa.gov 8/13
  !
  !-----------------------------------------------------------------------------
  ! Generic ESMF Main
  !-----------------------------------------------------------------------------

  use ESMF

  use modele_driver, only: &
    driver_SS => SetServices

  implicit none
  private
  public modelE_NUOPC_mainDriver

  contains

    subroutine modelE_NUOPC_mainDriver()
    
    integer                       :: rc, userRc
    type(ESMF_GridComp)           :: drvComp

    ! Initialize ESMF
    call ESMF_Initialize(defaultCalkind=ESMF_CALKIND_GREGORIAN, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      call ESMF_Finalize(endflag=ESMF_END_ABORT)
      
    call ESMF_LogWrite("ModelE STARTING", ESMF_LOGMSG_INFO, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      call ESMF_Finalize(endflag=ESMF_END_ABORT)
      
    !-----------------------------------------------------------------------------
    
    ! -> CREATE THE DRIVER
    drvComp = ESMF_GridCompCreate(name="driver", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      call ESMF_Finalize(endflag=ESMF_END_ABORT)
      
    ! -> SET DRIVER SERVICES
    call ESMF_GridCompSetServices(drvComp, driver_SS, userRc=userRc, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      call ESMF_Finalize(endflag=ESMF_END_ABORT)
    if (ESMF_LogFoundError(rcToCheck=userRc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      call ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! INITIALIZE THE DRIVER
    call ESMF_GridCompInitialize(drvComp, userRc=userRc, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      call ESMF_Finalize(endflag=ESMF_END_ABORT)
    if (ESMF_LogFoundError(rcToCheck=userRc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      call ESMF_Finalize(endflag=ESMF_END_ABORT)
        
    ! RUN THE DRIVER
    call ESMF_GridCompRun(drvComp, userRc=userRc, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      call ESMF_Finalize(endflag=ESMF_END_ABORT)
    if (ESMF_LogFoundError(rcToCheck=userRc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      call ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    ! FINALIZE THE DRIVER
    call ESMF_GridCompFinalize(drvComp, userRc=userRc, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      call ESMF_Finalize(endflag=ESMF_END_ABORT)
    if (ESMF_LogFoundError(rcToCheck=userRc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      call ESMF_Finalize(endflag=ESMF_END_ABORT)

    !-----------------------------------------------------------------------------
    
    call ESMF_LogWrite("ModelE FINISHED", ESMF_LOGMSG_INFO, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      call ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! Finalize ESMF
    call ESMF_Finalize()
  end subroutine modelE_NUOPC_mainDriver

end module main_nuopc_mod
#endif  
