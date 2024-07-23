!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module SP_ModProc
  implicit none
  SAVE

  ! MPI information
  !----------------------------------------------------------------------------
  integer:: iComm  = -1 ! For MPI communicators
  integer:: iProc  = -1 ! Current processor index
  integer:: nProc  = -1 ! Total processor number
  integer:: iError = -1 ! Error message
end module SP_ModProc
!==============================================================================
