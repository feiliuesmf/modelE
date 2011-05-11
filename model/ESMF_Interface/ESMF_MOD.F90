!***********************************************************
!* WARNING                                                 *
!*                                                         *
!* This file has been automatically generated via a rather *
!* perverse set of preprocessors.  Do _NOT_ directly edit  *
!* this file.  Instead, on should obtain the $?.nw (noweb) *
!* and make changes there.                                 *
!*                                                         *
!* For more information, contact:                          *
!*        Tom Clune <Thomas.L.Clune@nasa.gov>              *
!*                                                         *
!***********************************************************
Module ESMF_MOD_private
  Implicit None
  Private

  ! Derived Types
  Public :: Clock

  ! Constants
  Public :: ESMF_FAILURE
  Public :: ESMF_SUCCESS
  Integer, Parameter :: ESMF_FAILURE = -1
  Integer, Parameter :: ESMF_SUCCESS = 0

  type Clock
     Private
     Integer :: placeholder
  end type Clock

End Module ESMF_MOD_private

#ifndef USE_ESMF
Module ESMF_MOD
  Use ESMF_MOD_private, Only: ESMF_FAILURE, ESMF_SUCCESS
  use ESMF_Mod_private, only: ESMF_Clock => Clock

  Implicit None
  Private

  Public :: ESMF_Clock

  Public :: ESMF_Initialize, ESMF_Finalize

Contains

  Subroutine ESMF_Initialize(rc)
    Integer, Optional, Intent(Out) :: rc
    If (Present(rc)) rc = ESMF_SUCCESS
  End Subroutine ESMF_Initialize

  Subroutine ESMF_Finalize(rc)
    Integer, Optional, Intent(Out) :: rc
    If (Present(rc)) rc = ESMF_SUCCESS
  End Subroutine ESMF_Finalize

End Module ESMF_MOD
#endif
