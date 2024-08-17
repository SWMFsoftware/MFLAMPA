!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module SP_ModProc

  implicit none
  SAVE

  ! Public MPI information
  !----------------------------------------------------------------------------
  ! For MPI communicators
  integer :: iComm  = -1
  ! Total processor number
  integer :: nProc  = -1
  ! Current processor index
  integer :: iProc  = -1
  ! Error message
  integer :: iError = -1

  ! Processors needed to work on the same field line, if nProc > nLine
  ! For communicators on the same field line
  integer :: iCommSameLine = -1
  ! Total processor number in iCommSameLine, on the same field line
  integer :: nProcSameLine = -1
  ! Processor index (from 0) on the same field line
  integer :: iProcSameLine0 = -1
  ! Processor index (Not from 0, but just iProc) on field line
  integer, allocatable :: iProcSameLine_I(:)

end module SP_ModProc
!==============================================================================
